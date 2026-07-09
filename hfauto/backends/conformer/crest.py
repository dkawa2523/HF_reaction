from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math
import re

from hfauto.backends.conformer.rdkit import RDKitConformerBackend
from hfauto.chemistry.hf_builder import XYZ, read_xyz, write_xyz
from hfauto.core.executables import resolve_executable, run_command
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.units import HARTREE_TO_KCAL_MOL, R_KCAL_MOL_K


@dataclass
class ParsedConformer:
    xyz: XYZ
    energy: float | None = None
    relative_energy_kcal_mol: float = 0.0


def read_multixyz(path: str | Path) -> list[XYZ]:
    lines = Path(path).read_text(encoding="utf-8", errors="ignore").splitlines()
    out: list[XYZ] = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        try:
            n = int(lines[i].strip().split()[0])
        except Exception:
            break
        comment = lines[i + 1] if i + 1 < len(lines) else ""
        block = "\n".join([str(n), comment, *lines[i + 2 : i + 2 + n]]) + "\n"
        tmp = Path(path).with_suffix(".tmp.xyz")
        tmp.write_text(block, encoding="utf-8")
        out.append(read_xyz(tmp))
        try:
            tmp.unlink()
        except Exception:
            pass
        i += n + 2
    return out



def parse_multi_xyz(path: str | Path) -> list[XYZ]:
    """Backward-compatible alias for read_multixyz."""
    return read_multixyz(path)


def _floats(line: str) -> list[float]:
    return [float(x) for x in re.findall(r"[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?|[-+]?\d+(?:[Ee][-+]?\d+)", line)]


def parse_crest_energies(workdir: str | Path, conformers: list[XYZ] | None = None, n_expected: int | None = None) -> list[float | None]:
    """Parse CREST energies conservatively.

    CREST output variants differ by version. We first read crest.energies, then
    comments in crest_conformers.xyz. If values look like Hartree totals, they are
    converted to relative kcal/mol later; if they look like already relative kcal,
    they are kept as-is by the heuristic in ``assign_relative_energies``.
    """
    wd = Path(workdir)
    if wd.is_file():
        # Compatibility mode: directly parse a crest.energies-like file.
        raw_file: list[float | None] = []
        for line in wd.read_text(encoding="utf-8", errors="ignore").splitlines():
            vals = _floats(line)
            if vals:
                raw_file.append(vals[1] if len(vals) >= 2 and abs(vals[0] - round(vals[0])) < 1e-6 else vals[0])
        if n_expected is not None:
            return assign_relative_energies(raw_file, int(n_expected))
        return raw_file
    conformers = conformers or []
    candidates = [wd / "crest.energies", wd / "crest_conformers.xyz"]
    raw: list[float | None] = []
    if candidates[0].exists():
        for line in candidates[0].read_text(encoding="utf-8", errors="ignore").splitlines():
            vals = _floats(line)
            if not vals:
                continue
            # Common formats have rank/index first and energy/relative energy later.
            if len(vals) >= 2 and abs(vals[0] - round(vals[0])) < 1e-6:
                raw.append(vals[1])
            else:
                raw.append(vals[0])
    if len(raw) >= len(conformers):
        return raw[: len(conformers)]
    # Fall back to comment-line parsing.
    raw = []
    for xyz in conformers:
        vals = _floats(xyz.comment)
        raw.append(vals[0] if vals else None)
    return raw[: len(conformers)]


def assign_relative_energies(raw: list[float | None], n: int) -> list[float]:
    vals = [x for x in raw if x is not None]
    if not vals:
        return [0.0] * n
    # Negative large values are likely Hartree total energies.
    if min(vals) < -1.0:
        min_e = min(vals)
        rel = [0.0 if x is None else (x - min_e) * HARTREE_TO_KCAL_MOL for x in raw]
    else:
        # Positive small numbers are likely relative energies already in kcal/mol.
        min_e = min(vals)
        rel = [0.0 if x is None else max(0.0, x - min_e) for x in raw]
    if len(rel) < n:
        rel.extend([max(rel or [0.0]) + 999.0] * (n - len(rel)))
    return rel[:n]


class CRESTConformerBackend:
    """CREST/GFN-xTB conformer search adapter with RDKit fallback.

    Production execution is enabled with ``allow_subprocess: true``. Otherwise
    this backend returns RDKit conformers marked as fallback, so the whole package
    remains runnable in review and CI environments without CREST installed.
    """

    name = "crest"

    def __init__(self, **kwargs):
        self.defaults = kwargs

    def generate(self, molecule: Artifact, config: dict, workdir: str) -> list[Artifact]:
        cfg = {**self.defaults, **(config or {})}
        exe = resolve_executable("crest", cfg.get("executable"))
        allow = bool(cfg.get("allow_subprocess", cfg.get("mode", "fallback") in {"auto", "subprocess"}))
        if not allow or exe is None:
            artifacts = RDKitConformerBackend().generate(molecule, {**cfg, "source_fallback": "crest_unavailable_rdkit"}, workdir)
            for a in artifacts:
                if a.status.status == "success":
                    a.data["source"] = "rdkit_fallback_for_crest"
                    a.method = {**(a.method or {}), "requested_engine": self.name, "intended_engine": self.name, "fallback_engine": "rdkit"}
                    a.qc["fallback_reason"] = "crest_not_enabled_or_not_found"
                    a.qc["intended_engine"] = self.name
                    a.qc["crest_executable"] = exe
            return artifacts

        wd = Path(workdir) / molecule.data["mol_id"] / "crest"
        wd.mkdir(parents=True, exist_ok=True)
        # Use RDKit to create a robust initial geometry for CREST.
        seed_artifacts = RDKitConformerBackend().generate(molecule, {**cfg, "selected_max": 1}, str(wd / "seed"))
        seed = next((a for a in seed_artifacts if a.status.status == "success"), None)
        if seed is None:
            return seed_artifacts
        seed_xyz = wd / "input.xyz"
        seed_xyz.write_text(Path(seed.paths["xyz"]).read_text(encoding="utf-8"), encoding="utf-8")

        cmd = [exe, str(seed_xyz)]
        gfn = cfg.get("gfn", 2)
        if gfn is not None:
            cmd += [f"--gfn{int(gfn)}"]
        charge = int(molecule.data.get("formal_charge", 0) or 0)
        mult = int(molecule.data.get("multiplicity", 1) or 1)
        if charge:
            cmd += ["--chrg", str(charge)]
        if mult > 1:
            cmd += ["--uhf", str(mult - 1)]
        if cfg.get("quick", True):
            cmd.append("--quick")
        if "ewin" in cfg:
            cmd += ["--ewin", str(cfg["ewin"])]
        if "nci" in cfg and cfg["nci"]:
            cmd.append("--nci")
        threads = cfg.get("threads")
        env = {"OMP_NUM_THREADS": str(threads)} if threads else None
        result = run_command(cmd, wd, timeout_s=int(cfg.get("timeout_s", 7200)), env=env, stdout_name="crest.stdout", stderr_name="crest.stderr")
        if not result.ok:
            fail = Artifact.failure(
                f"crest_failed_{molecule.artifact_id}",
                "conformer",
                f"CREST failed returncode={result.returncode} timed_out={result.timed_out}",
                category="crest_failed",
                parents=[molecule.artifact_id],
                recoverable=True,
                recommended_fallback="rdkit_conformer_backend",
                command_result=result.to_dict(),
            )
            return [fail]

        conf_path = wd / "crest_conformers.xyz"
        if not conf_path.exists():
            return [
                Artifact.failure(
                    f"crest_failed_{molecule.artifact_id}",
                    "conformer",
                    "CREST completed but crest_conformers.xyz was not produced",
                    category="crest_missing_conformers",
                    parents=[molecule.artifact_id],
                    recommended_fallback="rdkit_conformer_backend",
                    command_result=result.to_dict(),
                )
            ]
        conformers = read_multixyz(conf_path)
        raw_e = parse_crest_energies(wd, conformers)
        rel_e = assign_relative_energies(raw_e, len(conformers))
        order = sorted(range(len(conformers)), key=lambda i: rel_e[i])
        max_selected = int(cfg.get("selected_max", cfg.get("max_conformers", 20)))
        energy_window = float(cfg.get("energy_window_kcal_mol", cfg.get("ewin", 8.0)))
        selected = [i for i in order if rel_e[i] <= energy_window][:max_selected]
        if not selected and order:
            selected = order[:1]
        temperature = float(cfg.get("boltzmann_temperature_K", 298.15))
        weights_raw = [math.exp(-rel_e[i] / (R_KCAL_MOL_K * temperature)) for i in selected]
        z = sum(weights_raw) or 1.0
        artifacts: list[Artifact] = []
        out_base = Path(workdir) / molecule.data["mol_id"]
        out_base.mkdir(parents=True, exist_ok=True)
        for rank, (i, wraw) in enumerate(zip(selected, weights_raw)):
            conf_id = f"{molecule.data['mol_id']}_conf{rank:04d}"
            xyz_path = out_base / f"conf{rank:04d}.xyz"
            xyz = conformers[i]
            xyz.comment = f"generated_by=hfauto source=crest crest_index={i} relE_kcal_mol={rel_e[i]:.6f}"
            write_xyz(xyz, xyz_path)
            rec = {
                "conformer_id": conf_id,
                "mol_id": molecule.data["mol_id"],
                "source": "crest",
                "relative_energy_kcal_mol": float(rel_e[i]),
                "boltzmann_weight_298K": float(wraw / z),
                "xyz_path": str(xyz_path),
                "selected_for_hf_build": True,
                "rmsd_cluster_id": f"crest_rank_{rank:04d}",
                "extras": {
                    "crest_index": int(i),
                    "raw_energy": raw_e[i] if i < len(raw_e) else None,
                    "command_result": result.to_dict(),
                    "seed_source": seed.artifact_id,
                },
            }
            artifacts.append(
                Artifact(
                    artifact_id=conf_id,
                    artifact_type="conformer",
                    parents=[molecule.artifact_id, seed.artifact_id],
                    paths={"xyz": str(xyz_path), "crest_conformers": str(conf_path), "stdout": result.stdout_path, "stderr": result.stderr_path},
                    data=rec,
                    method={"engine": self.name, "gfn": gfn, "command": result.command},
                    qc={
                        "generated": True,
                        "relative_energy_kcal_mol": float(rel_e[i]),
                        "boltzmann_weight_298K": float(wraw / z),
                        "conformer_engine": self.name,
                        "fallback_reason": None,
                    },
                )
            )
        return artifacts

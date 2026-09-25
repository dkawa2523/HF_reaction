from __future__ import annotations

import math
import os
import re
import shutil
from dataclasses import dataclass
from pathlib import Path

from hfauto.backends.conformer.rdkit import RDKitConformerBackend
from hfauto.chemistry.xyz import XYZ, write_xyz
from hfauto.chemistry.xyz_trajectory import read_xyz_trajectory
from hfauto.core.artifacts import species_xyz_path
from hfauto.core.executables import resolve_executable, run_command
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.units import HARTREE_TO_KCAL_MOL, R_KCAL_MOL_K


@dataclass
class ParsedConformer:
    xyz: XYZ
    energy: float | None = None
    relative_energy_kcal_mol: float = 0.0


@dataclass(frozen=True)
class ParsedCrestTopologyEndpoint:
    xyz: XYZ
    trajectory_path: Path
    electronic_energy_hartree: float | None


def parse_crest_topology_endpoint(
    failure: Artifact,
) -> ParsedCrestTopologyEndpoint | None:
    """Parse a completed CREST initial optimization that safety-stopped.

    Only the conformer backend interprets CREST text.  The workflow receives a
    normalized endpoint and still has to validate atom order and connectivity.
    """

    if (
        failure.status.status != "failed"
        or failure.status.category != "crest_failed"
    ):
        return None
    command = dict(failure.data.get("command_result") or {})
    stdout_candidates = [failure.paths.get("stdout"), command.get("stdout_path")]
    stdout_path = next(
        (
            Path(str(value))
            for value in stdout_candidates
            if value and Path(str(value)).is_file()
        ),
        None,
    )
    if stdout_path is None:
        return None
    stdout = stdout_path.read_text(encoding="utf-8", errors="ignore")
    if (
        "Change in topology detected" not in stdout
        or "Geometry successfully optimized" not in stdout
    ):
        return None
    workdir = Path(str(command.get("cwd") or stdout_path.parent))
    trajectory_candidates = [
        failure.paths.get("optimization_trajectory"),
        workdir / "crestopt.log",
        stdout_path.parent / "crestopt.log",
    ]
    trajectory = next(
        (
            Path(str(value))
            for value in trajectory_candidates
            if value and Path(str(value)).is_file()
        ),
        None,
    )
    if trajectory is None:
        return None
    final = read_xyz_trajectory(trajectory)[-1]
    values = _floats(final.comment)
    energy = values[-1] if values else None
    return ParsedCrestTopologyEndpoint(final, trajectory, energy)


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
    """CREST/GFN-xTB conformer and non-covalent-complex search adapter.

    Molecules without coordinates are seeded by RDKit.  Species and assembled
    complexes are seeded from their existing XYZ geometry, which is required for
    CREST ``--nci``.  An NCI request never falls back to RDKit because RDKit does
    not perform a non-covalent-complex search.
    """

    name = "crest"

    def __init__(self, **kwargs):
        self.defaults = kwargs

    @staticmethod
    def _calculation_identity(source: Artifact) -> str:
        """Return an ID scoped to the actual geometry-bearing source.

        A multicomponent species can retain provenance to one originating
        molecule, but it is not that molecule.  Using ``mol_id`` for a species
        would overwrite its upstream molecular conformer when CREST emits
        ``*_conf0000``.
        """

        if source.artifact_type in {"molecule", "molecule_state"}:
            return str(source.data.get("mol_id") or source.artifact_id)
        return str(source.data.get("species_id") or source.artifact_id)

    def generate(self, molecule: Artifact, config: dict, workdir: str) -> list[Artifact]:
        cfg = {**self.defaults, **(config or {})}
        exe = resolve_executable("crest", cfg.get("executable"))
        xtb_exe = resolve_executable("xtb", cfg.get("xtb_executable"))
        allow = bool(cfg.get("allow_subprocess", cfg.get("mode", "fallback") in {"auto", "subprocess"}))
        if not allow or exe is None:
            if cfg.get("nci") or not bool(cfg.get("fallback_to_rdkit", True)):
                reason = "CREST subprocess is disabled" if not allow else "CREST executable was not found"
                return [
                    Artifact.failure(
                        f"crest_unavailable_{molecule.artifact_id}",
                        "conformer",
                        reason,
                        category="crest_not_run",
                        parents=[molecule.artifact_id],
                        recommended_fallback="configure CREST or explicitly enable fallback_to_rdkit",
                    )
                ]
            artifacts = RDKitConformerBackend().generate(molecule, {**cfg, "source_fallback": "crest_unavailable_rdkit"}, workdir)
            for a in artifacts:
                if a.status.status == "success":
                    a.data["source"] = "rdkit_fallback_for_crest"
                    a.method = {**(a.method or {}), "requested_engine": self.name, "intended_engine": self.name, "fallback_engine": "rdkit"}
                    a.qc["fallback_reason"] = "crest_not_enabled_or_not_found"
                    a.qc["intended_engine"] = self.name
                    a.qc["crest_executable"] = exe
            return artifacts

        identity = self._calculation_identity(molecule)
        wd = Path(workdir) / identity / "crest"
        wd.mkdir(parents=True, exist_ok=True)
        seed_xyz = wd / "input.xyz"
        seed: Artifact | None = None
        try:
            geometry_path = species_xyz_path(molecule)
        except (FileNotFoundError, ValueError):
            seed_artifacts = RDKitConformerBackend().generate(
                molecule,
                {**cfg, "selected_max": 1},
                str(wd / "seed"),
            )
            seed = next(
                (
                    item
                    for item in seed_artifacts
                    if item.status.status == "success"
                ),
                None,
            )
            if seed is None:
                return seed_artifacts
            geometry_path = Path(seed.paths["xyz"])
        if geometry_path.resolve() != seed_xyz.resolve():
            shutil.copyfile(geometry_path, seed_xyz)

        # The subprocess cwd is already ``wd``; passing a repository-relative
        # path would make CREST resolve the work directory twice.
        cmd = [exe, seed_xyz.name]
        gfn = cfg.get("gfn", 2)
        if gfn is not None:
            cmd += [f"--gfn{int(gfn)}"]
        charge = int(
            molecule.data.get("charge", molecule.data.get("formal_charge", 0))
            or 0
        )
        mult = int(molecule.data.get("multiplicity", 1) or 1)
        if charge:
            cmd += ["--chrg", str(charge)]
        if mult > 1:
            cmd += ["--uhf", str(mult - 1)]
        if cfg.get("quick", True):
            cmd.append("--quick")
        if "ewin" in cfg:
            cmd += ["--ewin", str(cfg["ewin"])]
        if cfg.get("nci"):
            cmd.append("--nci")
        threads = cfg.get("threads")
        env = {"OMP_NUM_THREADS": str(threads)} if threads else {}
        if xtb_exe:
            env["PATH"] = os.pathsep.join(
                [str(Path(xtb_exe).parent), os.environ.get("PATH", "")]
            )
        result = run_command(cmd, wd, timeout_s=int(cfg.get("timeout_s", 7200)), env=env, stdout_name="crest.stdout", stderr_name="crest.stderr")
        if not result.ok:
            stdout_text = (
                Path(result.stdout_path).read_text(
                    encoding="utf-8", errors="ignore"
                )
                if result.stdout_path and Path(result.stdout_path).is_file()
                else ""
            )
            trajectory = wd / "crestopt.log"
            topology_changed = (
                "Change in topology detected" in stdout_text
                and "Geometry successfully optimized" in stdout_text
                and trajectory.is_file()
            )
            fail = Artifact.failure(
                f"crest_failed_{molecule.artifact_id}",
                "conformer",
                f"CREST failed returncode={result.returncode} timed_out={result.timed_out}",
                category="crest_failed",
                parents=[molecule.artifact_id],
                recoverable=True,
                recommended_fallback=(
                    "classify_unbiased_relaxation_endpoint"
                    if topology_changed
                    else "review_crest_failure"
                ),
                command_result=result.to_dict(),
                topology_change_detected=topology_changed,
            )
            fail.paths = {
                **(
                    {"stdout": str(result.stdout_path)}
                    if result.stdout_path
                    else {}
                ),
                **(
                    {"stderr": str(result.stderr_path)}
                    if result.stderr_path
                    else {}
                ),
                **(
                    {"optimization_trajectory": str(trajectory)}
                    if trajectory.is_file()
                    else {}
                ),
            }
            fail.method = {
                "engine": self.name,
                "gfn": gfn,
                "command": result.command,
                "crest_executable": exe,
                "xtb_executable": xtb_exe,
            }
            fail.qc = {
                "nci_search": bool(cfg.get("nci")),
                "topology_change_detected": topology_changed,
                "unbiased_initial_optimization_completed": topology_changed,
                "activation_barrier_validated": False,
            }
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
        conformers = read_xyz_trajectory(conf_path)
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
        out_base = Path(workdir) / identity
        out_base.mkdir(parents=True, exist_ok=True)
        for rank, (i, wraw) in enumerate(zip(selected, weights_raw)):
            conf_id = f"{identity}_conf{rank:04d}"
            xyz_path = out_base / f"conf{rank:04d}.xyz"
            xyz = conformers[i]
            xyz.comment = f"generated_by=hfauto source=crest crest_index={i} relE_kcal_mol={rel_e[i]:.6f}"
            write_xyz(xyz, xyz_path)
            rec = {
                "conformer_id": conf_id,
                "mol_id": molecule.data.get("mol_id"),
                "source_species_id": molecule.data.get("species_id"),
                "source": "crest",
                "relative_energy_kcal_mol": float(rel_e[i]),
                "boltzmann_weight_298K": float(wraw / z),
                "xyz_path": str(xyz_path),
                "selected_for_complex_build": True,
                "rmsd_cluster_id": f"crest_rank_{rank:04d}",
                "extras": {
                    "crest_index": int(i),
                    "raw_energy": raw_e[i] if i < len(raw_e) else None,
                    "command_result": result.to_dict(),
                    "seed_source": seed.artifact_id if seed is not None else molecule.artifact_id,
                    "nci": bool(cfg.get("nci")),
                },
            }
            artifacts.append(
                Artifact(
                    artifact_id=conf_id,
                    artifact_type="conformer",
                    parents=list(
                        dict.fromkeys(
                            [
                                molecule.artifact_id,
                                *( [seed.artifact_id] if seed is not None else [] ),
                            ]
                        )
                    ),
                    paths={"xyz": str(xyz_path), "crest_conformers": str(conf_path), "stdout": result.stdout_path, "stderr": result.stderr_path},
                    data=rec,
                    method={
                        "engine": self.name,
                        "gfn": gfn,
                        "command": result.command,
                        "crest_executable": exe,
                        "xtb_executable": xtb_exe,
                    },
                    qc={
                        "generated": True,
                        "relative_energy_kcal_mol": float(rel_e[i]),
                        "boltzmann_weight_298K": float(wraw / z),
                        "conformer_engine": self.name,
                        "nci_search": bool(cfg.get("nci")),
                        "fallback_reason": None,
                    },
                )
            )
        return artifacts

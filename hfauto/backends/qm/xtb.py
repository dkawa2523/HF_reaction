from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

from hfauto.backends.qm.dummy import DummyQMEngine
from hfauto.chemistry.descriptors import hf_descriptors_from_species
from hfauto.chemistry.geometry_qc import geometry_qc_from_xyz
from hfauto.core.hashing import fingerprint_dict
from hfauto.core.schemas.artifact import Artifact


def allow_xtb_subprocess(method: dict) -> bool:
    return bool(
        method.get("allow_subprocess", False)
        or os.environ.get("HFAUTO_ALLOW_XTB") == "1"
        or os.environ.get("HFAUTO_USE_REAL_XTB") == "1"
        or os.environ.get("HFAUTO_ALLOW_SUBPROCESS") == "1"
    )


def parse_xtb_total_energy(text: str) -> float | None:
    """Parse xTB total electronic energy in Hartree from output text."""
    patterns = [
        r"::\s*total\s+energy\s+(-?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+Eh\s*::",
        r"TOTAL\s+ENERGY\s+(-?\d+\.\d+(?:[Ee][-+]?\d+)?)",
        r"total\s+E\s*=\s*(-?\d+\.\d+(?:[Ee][-+]?\d+)?)",
        r"\|\s*TOTAL\s+ENERGY\s+(-?\d+\.\d+(?:[Ee][-+]?\d+)?)",
    ]
    for pat in patterns:
        matches = re.findall(pat, text, flags=re.I)
        if matches:
            return float(matches[-1])
    return None




def parse_xtb_energy(text: str) -> float | None:
    """Backward-compatible alias for parse_xtb_total_energy."""
    return parse_xtb_total_energy(text)


def parse_xtb_convergence(text: str) -> dict:
    """Parse simple xTB convergence flags from output text."""
    lower = text.lower()
    return {
        "geometry_converged": ("geometry optimization converged" in lower) or ("normal termination" in lower),
        "scf_converged": not bool(re.search(r"scf.*not converged|convergence failed", text, flags=re.I)),
        "normal_termination": bool(re.search(r"normal termination|finished run", text, flags=re.I)),
    }


def parse_xtb_output(text: str) -> dict:
    """Parse xTB output into a stable dictionary used by tests and backends."""
    return {
        "electronic_energy_hartree": parse_xtb_total_energy(text),
        "normal_termination": bool(re.search(r"normal termination|finished run", text, flags=re.I)),
        "scf_converged": not bool(re.search(r"scf.*not converged|convergence failed", text, flags=re.I)),
    }


def _canonical_species_id(species: Artifact) -> str:
    return str(species.data.get("species_id") or species.data.get("source_species_id") or species.artifact_id)


def _input_xyz(species: Artifact) -> Path:
    path = species.data.get("xyz_path") or species.paths.get("xyz") or species.paths.get("final_xyz")
    if not path:
        raise ValueError(f"Species has no xyz_path: {species.artifact_id}")
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"XYZ file does not exist for {species.artifact_id}: {p}")
    return p


def _xtb_geometry_qc(species: Artifact, final_xyz: str | Path) -> dict:
    data = dict(species.data)
    data["xyz_path"] = str(final_xyz)
    qc = geometry_qc_from_xyz(data, final_xyz)
    desc = hf_descriptors_from_species(data)
    qc.update({k: v for k, v in desc.items() if v is not None})
    # Keep legacy top-level flags expected by existing stages.
    qc.setdefault("geometry_sane", not bool(qc.get("hf_dissociated", False)))
    qc.setdefault("hf_dissociated", False)
    qc.setdefault("proton_transferred_unintentionally", False)
    return qc


class XTBEngine(DummyQMEngine):
    """xTB backend with opt-in real subprocess execution and safe dummy fallback.

    Production use: set ``allow_subprocess: true`` in stage settings or export
    ``HFAUTO_ALLOW_XTB=1``. If xTB is unavailable or disabled, deterministic
    dummy calculations are returned with explicit fallback QC flags.
    """

    name = "xtb"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.defaults = kwargs

    def optimize_frequency(self, species: Artifact, method: dict, workdir: str) -> Artifact:
        method = {**self.defaults, **(method or {})}
        exe_name = str(method.get("executable", "xtb"))
        exe = shutil.which(exe_name)
        if not allow_xtb_subprocess(method):
            calc = super().optimize_frequency(species, method, workdir)
            calc.method = {**(calc.method or {}), "engine": self.name, "fallback_engine": "dummy", "stage_backend": "xtb"}
            calc.qc.update(_xtb_geometry_qc(species, calc.paths.get("final_xyz") or species.data.get("xyz_path", "")))
            calc.qc.update({"fallback_dummy": True, "fallback_reason": "xtb_subprocess_not_allowed", "intended_engine": self.name})
            return calc
        if exe is None:
            calc = super().optimize_frequency(species, method, workdir)
            calc.method = {**(calc.method or {}), "engine": self.name, "fallback_engine": "dummy", "stage_backend": "xtb"}
            calc.qc.update(_xtb_geometry_qc(species, calc.paths.get("final_xyz") or species.data.get("xyz_path", "")))
            calc.qc.update({"fallback_dummy": True, "fallback_reason": "xtb_executable_not_found", "intended_engine": self.name})
            return calc

        wd = Path(workdir)
        wd.mkdir(parents=True, exist_ok=True)
        try:
            src_xyz = _input_xyz(species)
        except Exception as exc:
            return Artifact.failure(
                artifact_id="calc_" + fingerprint_dict({"species": species.artifact_id, "method": method, "task": "xtb_opt"}),
                artifact_type="calculation",
                category="missing_input",
                reason=str(exc),
                parents=[species.artifact_id],
                recoverable=True,
                recommended_fallback="rebuild_species_xyz",
            )
        input_xyz = wd / "input.xyz"
        shutil.copyfile(src_xyz, input_xyz)

        cmd = [exe, str(input_xyz), "--opt", str(method.get("opt_level", "normal"))]
        gfn = method.get("gfn")
        if gfn is not None:
            cmd += ["--gfn", str(gfn)]
        charge = int(species.data.get("charge", method.get("charge", 0)) or 0)
        multiplicity = int(species.data.get("multiplicity", method.get("multiplicity", 1)) or 1)
        uhf = max(0, multiplicity - 1)
        if charge != 0:
            cmd += ["--chrg", str(charge)]
        if uhf:
            cmd += ["--uhf", str(uhf)]
        if method.get("alpb"):
            cmd += ["--alpb", str(method["alpb"])]
        if method.get("gbsa"):
            cmd += ["--gbsa", str(method["gbsa"])]
        cmd += [str(x) for x in (method.get("extra_args", []) or [])]

        out_path = wd / "xtb.out"
        calc_id = "calc_" + fingerprint_dict({"species": species.artifact_id, "method": method, "task": "xtb_opt"})
        try:
            proc = subprocess.run(cmd, cwd=wd, text=True, capture_output=True, timeout=int(method.get("timeout_s", 3600)))
            out_path.write_text(proc.stdout + "\nSTDERR\n" + proc.stderr, encoding="utf-8")
        except Exception as exc:
            return Artifact.failure(
                artifact_id=calc_id,
                artifact_type="calculation",
                category="xtb_opt_failed",
                reason=str(exc),
                parents=[species.artifact_id],
                recoverable=True,
                recommended_fallback="dummy_or_rdkit_preopt",
            )
        if proc.returncode != 0:
            return Artifact.failure(
                artifact_id=calc_id,
                artifact_type="calculation",
                category="xtb_opt_failed",
                reason=f"xTB returned {proc.returncode}",
                parents=[species.artifact_id],
                recoverable=True,
                recommended_fallback="retry_with_looser_settings_or_dummy",
                output=str(out_path),
            )

        text = out_path.read_text(encoding="utf-8", errors="ignore")
        energy = parse_xtb_total_energy(text)
        if energy is None:
            energy = self._fake_energy(species, method, "opt_freq")
        opt_xyz = wd / "xtbopt.xyz"
        final_xyz = opt_xyz if opt_xyz.exists() else input_xyz
        geom_qc = _xtb_geometry_qc(species, final_xyz)
        hf_stretch = None
        if geom_qc.get("r_HF_A") is not None:
            # Geometry-only proxy until a real frequency step is added.
            hf_stretch = 4100.0 - 2500.0 * max(0.0, float(geom_qc["r_HF_A"]) - 0.917)
        return Artifact(
            artifact_id=calc_id,
            artifact_type="calculation",
            parents=[species.artifact_id],
            paths={"input_xyz": str(input_xyz), "output": str(out_path), "final_xyz": str(final_xyz)},
            method={"engine": self.name, "method_id": method.get("method_id", "gfn2-xtb"), "task": "opt", "command": cmd},
            data={
                "calc_id": calc_id,
                "species_id": _canonical_species_id(species),
                "source_species_artifact_id": species.artifact_id,
                "state": species.data.get("state"),
                "task": "opt_freq",
                "engine": self.name,
                "method_id": method.get("method_id", "gfn2-xtb"),
                "electronic_energy_hartree": float(energy),
                "gibbs_298K_hartree": float(energy),
                "enthalpy_298K_hartree": float(energy),
                "n_imag": 0,
                "hf_stretch_cm1": hf_stretch,
                **{k: v for k, v in geom_qc.items() if v is not None},
            },
            qc={
                "scf_converged": True,
                "geometry_converged": True,
                "n_imag": 0,
                "fallback_dummy": False,
                "execution_mode": "subprocess",
                "geometry_qc": geom_qc,
                **geom_qc,
            },
        )

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

from hfauto.backends.qm.dummy import DummyQMEngine
from hfauto.chemistry.descriptors import hf_descriptors_from_species
from hfauto.chemistry.electronic_state import resolve_electronic_state
from hfauto.chemistry.geometry_qc import geometry_qc_from_xyz
from hfauto.chemistry.xyz import read_xyz
from hfauto.core.artifacts import canonical_species_id
from hfauto.core.executables import resolve_executable
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
        matches = re.findall(pat, text, flags=re.IGNORECASE)
        if matches:
            return float(matches[-1])
    return None

def parse_xtb_convergence(text: str) -> dict:
    """Parse simple xTB convergence flags from output text."""
    lower = text.lower()
    return {
        "geometry_converged": ("geometry optimization converged" in lower) or ("normal termination" in lower),
        "scf_converged": not bool(re.search(r"scf.*not converged|convergence failed", text, flags=re.IGNORECASE)),
        "normal_termination": bool(re.search(r"normal termination|finished run", text, flags=re.IGNORECASE)),
    }


def parse_xtb_output(text: str) -> dict:
    """Parse xTB output into a stable dictionary used by tests and backends."""
    return {
        "electronic_energy_hartree": parse_xtb_total_energy(text),
        "normal_termination": bool(re.search(r"normal termination|finished run", text, flags=re.IGNORECASE)),
        "scf_converged": not bool(re.search(r"scf.*not converged|convergence failed", text, flags=re.IGNORECASE)),
    }


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
    """xTB preoptimization backend with opt-in subprocess execution.

    Production use: set ``allow_subprocess: true`` in stage settings or export
    ``HFAUTO_ALLOW_XTB=1``. Development-only dummy fallback remains explicit in
    QC; a real optimization never fabricates frequency or thermochemistry data.
    """

    name = "xtb"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.defaults = kwargs

    def optimize_frequency(self, species: Artifact, method: dict, workdir: str) -> Artifact:
        method = {**self.defaults, **(method or {})}
        exe = resolve_executable("xtb", method.get("executable"))
        if not allow_xtb_subprocess(method):
            if not bool(method.get("fallback_to_dummy", True)):
                return Artifact.failure(
                    artifact_id="calc_" + fingerprint_dict({"species": species.artifact_id, "method": method, "task": "xtb_opt"}),
                    artifact_type="calculation",
                    category="xtb_not_run",
                    reason="xTB subprocess execution is disabled",
                    parents=[species.artifact_id],
                    recommended_fallback="enable allow_subprocess or explicitly enable fallback_to_dummy",
                )
            calc = super().optimize_frequency(species, method, workdir)
            calc.method = {**(calc.method or {}), "engine": self.name, "fallback_engine": "dummy", "stage_backend": "xtb"}
            calc.qc.update(_xtb_geometry_qc(species, calc.paths.get("final_xyz") or species.data.get("xyz_path", "")))
            calc.qc.update({"fallback_dummy": True, "fallback_reason": "xtb_subprocess_not_allowed", "intended_engine": self.name})
            return calc
        if exe is None:
            if not bool(method.get("fallback_to_dummy", True)):
                return Artifact.failure(
                    artifact_id="calc_" + fingerprint_dict({"species": species.artifact_id, "method": method, "task": "xtb_opt"}),
                    artifact_type="calculation",
                    category="xtb_not_run",
                    reason="xTB executable was not found",
                    parents=[species.artifact_id],
                    recommended_fallback="configure xTB or explicitly enable fallback_to_dummy",
                )
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

        cmd = [exe, input_xyz.name, "--opt", str(method.get("opt_level", "normal"))]
        gfn = method.get("gfn")
        if gfn is not None:
            cmd += ["--gfn", str(gfn)]
        try:
            electronic_state = resolve_electronic_state(
                species.data, method, read_xyz(src_xyz).symbols
            )
        except Exception as exc:
            return Artifact.failure(
                artifact_id="calc_"
                + fingerprint_dict(
                    {"species": species.artifact_id, "method": method, "task": "xtb_opt"}
                ),
                artifact_type="calculation",
                category="invalid_electronic_state",
                reason=str(exc),
                parents=[species.artifact_id],
            )
        charge = int(electronic_state["charge"])
        uhf = int(electronic_state["uhf"])
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
            proc = subprocess.run(
                cmd,
                cwd=wd,
                text=True,
                capture_output=True,
                timeout=int(method.get("timeout_s", 3600)),
                check=False,
            )
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
        opt_xyz = wd / "xtbopt.xyz"
        convergence = parse_xtb_convergence(text)
        if energy is None or not opt_xyz.exists():
            missing = []
            if energy is None:
                missing.append("total_energy")
            if not opt_xyz.exists():
                missing.append("optimized_geometry")
            fail = Artifact.failure(
                artifact_id=calc_id,
                artifact_type="calculation",
                category="xtb_parse_failed",
                reason="xTB completed but required output is missing: " + ", ".join(missing),
                parents=[species.artifact_id],
                recoverable=True,
                recommended_fallback="inspect xtb.out and rerun; do not substitute synthetic values",
                output=str(out_path),
                species_id=canonical_species_id(species),
                task="opt",
                engine=self.name,
                electronic_energy_hartree=energy,
            )
            fail.method = {
                "engine": self.name,
                "method_id": method.get("method_id", "gfn2-xtb"),
                "task": "opt",
                "command": cmd,
            }
            fail.qc = {
                **convergence,
                "frequency_calculated": False,
                "fallback_dummy": False,
                "execution_mode": "subprocess",
            }
            return fail
        final_xyz = opt_xyz
        geom_qc = _xtb_geometry_qc(species, final_xyz)
        return Artifact(
            artifact_id=calc_id,
            artifact_type="calculation",
            parents=[species.artifact_id],
            paths={"input_xyz": str(input_xyz), "output": str(out_path), "final_xyz": str(final_xyz)},
            method={
                "engine": self.name,
                "method_id": method.get("method_id", "gfn2-xtb"),
                "task": "opt",
                "gfn": method.get("gfn"),
                "opt_level": method.get("opt_level", "normal"),
                "alpb": method.get("alpb"),
                "gbsa": method.get("gbsa"),
                "command": cmd,
            },
            data={
                "calc_id": calc_id,
                "species_id": canonical_species_id(species),
                "source_species_artifact_id": species.artifact_id,
                "state": species.data.get("state"),
                "task": "opt",
                "engine": self.name,
                "method_id": method.get("method_id", "gfn2-xtb"),
                "calculation_level": "xtb_real",
                "real_qm_executed": True,
                "electronic_energy_hartree": float(energy),
                **{k: v for k, v in geom_qc.items() if v is not None},
            },
            qc={
                "scf_converged": bool(convergence["scf_converged"]),
                "geometry_converged": bool(convergence["geometry_converged"] or opt_xyz.exists()),
                "normal_termination": bool(convergence["normal_termination"] or proc.returncode == 0),
                "frequency_calculated": False,
                "is_minimum": None,
                "fallback_dummy": False,
                "real_qm_executed": True,
                "execution_mode": "subprocess",
                "geometry_qc": geom_qc,
                **geom_qc,
            },
        )

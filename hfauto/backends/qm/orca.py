from __future__ import annotations

"""ORCA QM backend for DFT opt/freq and single-point calculations.

The backend is intentionally auditable and conservative:
- ORCA subprocess execution is opt-in via config or environment variables;
- input files are always written, even when execution is disabled;
- output parsing is isolated in small functions with parser tests;
- missing ORCA can either produce a recoverable failure Artifact or an explicit
  dummy fallback Artifact, depending on configuration.

This keeps Phase 4 useful both on developer laptops without ORCA and on HPC
nodes where ORCA is installed.
"""

from dataclasses import asdict
import os
import re
from pathlib import Path
from typing import Any

import numpy as np

from hfauto.backends.qm.dummy import DummyQMEngine
from hfauto.chemistry.descriptors import hf_descriptors_from_species
from hfauto.chemistry.geometry_qc import geometry_qc_from_xyz
from hfauto.chemistry.hf_builder import XYZ, read_xyz, write_xyz
from hfauto.core.executables import resolve_executable, run_command
from hfauto.core.hashing import fingerprint_dict
from hfauto.core.schemas.artifact import Artifact


def allow_orca_subprocess(method: dict[str, Any]) -> bool:
    """Return True only when ORCA execution is explicitly enabled."""
    return bool(
        method.get("allow_subprocess", False)
        or os.environ.get("HFAUTO_ALLOW_ORCA") == "1"
        or os.environ.get("HFAUTO_USE_REAL_ORCA") == "1"
        or os.environ.get("HFAUTO_ALLOW_SUBPROCESS") == "1"
    )


def _last_float(patterns: list[str], text: str) -> float | None:
    for pat in patterns:
        matches = re.findall(pat, text, flags=re.I | re.M)
        if matches:
            val = matches[-1]
            if isinstance(val, tuple):
                val = next((x for x in val if x not in {"", None}), val[-1])
            try:
                return float(val)
            except Exception:
                continue
    return None


def parse_orca_frequencies(text: str) -> list[float]:
    """Parse vibrational frequencies in cm^-1 from common ORCA text formats."""
    freqs: list[float] = []
    # Typical ORCA frequency table: "   6:       123.45 cm**-1"
    for m in re.findall(r"^\s*\d+\s*:\s*(-?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)\s+cm\*\*-1", text, flags=re.M):
        try:
            freqs.append(float(m))
        except ValueError:
            pass
    if freqs:
        return freqs

    # Fallback for compact/sample outputs containing "Frequencies -- ...".
    for line in text.splitlines():
        if re.search(r"frequenc", line, flags=re.I):
            vals = re.findall(r"-?\d+\.\d+|-?\d+", line)
            for v in vals:
                try:
                    x = float(v)
                    if abs(x) > 1.0:  # avoid atom indices or tiny formatting numbers
                        freqs.append(x)
                except ValueError:
                    pass
    return freqs


def parse_orca_cartesian_blocks(text: str) -> list[XYZ]:
    """Extract CARTESIAN COORDINATES blocks printed by ORCA.

    ORCA may not write a final ``.xyz`` depending on input keywords. This parser
    provides a fallback final geometry for downstream descriptor/QC stages.
    """
    lines = text.splitlines()
    blocks: list[XYZ] = []
    i = 0
    header_re = re.compile(r"CARTESIAN COORDINATES\s*\(ANGSTROEM\)", re.I)
    coord_re = re.compile(
        r"^\s*([A-Z][a-z]?)\s+(-?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)\s+(-?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)\s+(-?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)"
    )
    while i < len(lines):
        if not header_re.search(lines[i]):
            i += 1
            continue
        i += 1
        # Skip separator/blank lines.
        while i < len(lines) and (not lines[i].strip() or set(lines[i].strip()) <= {"-"}):
            i += 1
        symbols: list[str] = []
        coords: list[list[float]] = []
        while i < len(lines):
            line = lines[i]
            if not line.strip() or set(line.strip()) <= {"-"}:
                break
            m = coord_re.match(line)
            if not m:
                # End block on the first non-coordinate line after data starts.
                if symbols:
                    break
                i += 1
                continue
            symbols.append(m.group(1))
            coords.append([float(m.group(2)), float(m.group(3)), float(m.group(4))])
            i += 1
        if symbols:
            blocks.append(XYZ(symbols=symbols, coords=np.array(coords, dtype=float), comment="extracted_from_orca_output"))
        i += 1
    return blocks


def parse_orca_output(text: str | Path) -> dict[str, Any]:
    """Parse stable calculation fields from ORCA output text or an output path.

    The parser is deliberately tolerant; it returns ``None`` for absent fields
    rather than raising. Backend execution/QC decides whether missing fields are
    fatal for a specific task.
    """
    if isinstance(text, Path):
        text = text.read_text(encoding="utf-8", errors="ignore")
    energy = _last_float([r"FINAL\s+SINGLE\s+POINT\s+ENERGY\s+(-?\d+\.\d+(?:[Ee][-+]?\d+)?)"], text)
    zpe = _last_float(
        [
            r"Zero\s+point\s+energy\s+\.\.\.\s+(-?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+Eh",
            r"Zero\s+point\s+energy\s*[:=]\s*(-?\d+\.\d+(?:[Ee][-+]?\d+)?)",
            r"ZPE\s*(?:correction)?\s*[:=]\s*(-?\d+\.\d+(?:[Ee][-+]?\d+)?)",
        ],
        text,
    )
    enthalpy = _last_float(
        [
            r"Total\s+Enthalpy\s+\.\.\.\s+(-?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+Eh",
            r"Final\s+enthalpy\s+\.\.\.\s+(-?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+Eh",
            r"Total\s+Enthalpy\s*[:=]\s*(-?\d+\.\d+(?:[Ee][-+]?\d+)?)",
        ],
        text,
    )
    gibbs = _last_float(
        [
            r"Final\s+Gibbs\s+free\s+energy\s+\.\.\.\s+(-?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+Eh",
            r"Final\s+Gibbs\s+free\s+energy\s*[:=]\s*(-?\d+\.\d+(?:[Ee][-+]?\d+)?)",
            r"Gibbs\s+free\s+energy\s+\.\.\.\s+(-?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+Eh",
        ],
        text,
    )
    g_corr = _last_float(
        [
            r"G-E\(el\)\s+\.\.\.\s+(-?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+Eh",
            r"Total\s+Gibbs\s+correction\s*[:=]\s*(-?\d+\.\d+(?:[Ee][-+]?\d+)?)",
        ],
        text,
    )
    h_corr = _last_float(
        [
            r"H-E\(el\)\s+\.\.\.\s+(-?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+Eh",
            r"Thermal\s+Enthalpy\s+correction\s*[:=]\s*(-?\d+\.\d+(?:[Ee][-+]?\d+)?)",
        ],
        text,
    )
    if gibbs is None and energy is not None and g_corr is not None:
        gibbs = energy + g_corr
    if enthalpy is None and energy is not None and h_corr is not None:
        enthalpy = energy + h_corr

    freqs = parse_orca_frequencies(text)
    imag = [f for f in freqs if f < 0.0]
    positive = [f for f in freqs if f > 0.0]
    hf_stretch = max(positive) if positive else None
    version = None
    m = re.search(r"Program\s+Version\s+([0-9][^\s]*)", text, flags=re.I) or re.search(r"ORCA\s+VERSION\s+([0-9][^\s]*)", text, flags=re.I)
    if m:
        version = m.group(1)
    dipole = _last_float(
        [
            r"Magnitude\s+\(Debye\)\s*:\s*(-?\d+\.\d+(?:[Ee][-+]?\d+)?)",
            r"Total\s+Dipole\s+Moment\s*[:=].*?(-?\d+\.\d+(?:[Ee][-+]?\d+)?)\s*Debye",
        ],
        text,
    )
    normal = bool(re.search(r"ORCA\s+TERMINATED\s+NORMALLY", text, flags=re.I))
    scf_bad = bool(re.search(r"SCF\s+NOT\s+CONVERGED|SCF.*FAILED|convergence\s+failure", text, flags=re.I))
    scf_good = bool(re.search(r"SCF\s+CONVERGED|SCF.*CONVERGED\s+AFTER", text, flags=re.I)) or (energy is not None and not scf_bad)
    opt_good = bool(re.search(r"THE\s+OPTIMIZATION\s+HAS\s+CONVERGED|OPTIMIZATION\s+CONVERGED|OPTIMIZATION\s+RUN\s+DONE", text, flags=re.I))
    opt_bad = bool(re.search(r"OPTIMIZATION\s+DID\s+NOT\s+CONVERGE|GEOMETRY\s+OPTIMIZATION\s+FAILED", text, flags=re.I))
    return {
        "program_version": version,
        "electronic_energy_hartree": energy,
        "zpe_hartree": zpe,
        "enthalpy_298K_hartree": enthalpy,
        "gibbs_298K_hartree": gibbs,
        "thermal_correction_gibbs_hartree": None if gibbs is None or energy is None else gibbs - energy,
        "thermal_correction_enthalpy_hartree": None if enthalpy is None or energy is None else enthalpy - energy,
        "frequencies_cm1": freqs,
        "n_imag": len(imag),
        "lowest_freq_cm1": min(freqs) if freqs else None,
        "imag_freq_cm1": min(imag) if imag else None,
        "hf_stretch_cm1": hf_stretch,
        "dipole_D": dipole,
        "normal_termination": normal,
        "scf_converged": scf_good,
        "geometry_converged": opt_good and not opt_bad,
    }


def _canonical_species_id(species: Artifact) -> str:
    return str(species.data.get("species_id") or species.data.get("source_species_id") or species.artifact_id)


class OrcaInputRenderer:
    """Render ORCA input from a Species Artifact and method dictionary."""

    @staticmethod
    def keywords(method: dict[str, Any], task: str) -> list[str]:
        kw = list(method.get("orca_keywords", []) or [])
        if not kw:
            method_id = str(method.get("method_id", "r2SCAN-3c"))
            if method_id.lower() in {"r2scan3c", "r2scan-3c", "r2scan3c_orca"}:
                kw = ["r2SCAN-3c"]
            elif method_id.lower() in {"wb97xd4_def2_tzvppd_orca", "wb97xd4-def2tzvppd", "wb97x-d4/def2-tzvppd"}:
                kw = ["wB97X-D4", "def2-TZVPPD"]
            else:
                kw = [method_id]
        lower_tokens = [str(x).lower() for x in kw]
        lower_joined = " ".join(lower_tokens)
        if task == "opt_freq":
            # TightOpt/LooseOpt are convergence settings, not a substitute for
            # the Opt job keyword; append Opt unless an actual optimization or
            # path-search keyword is present.
            has_opt_job = any(k in {"opt", "optts"} or "neb" in k for k in lower_tokens)
            if not has_opt_job:
                kw.append("Opt")
                lower_tokens.append("opt")
            if not any(k in {"freq", "numfreq"} for k in lower_tokens):
                kw.append("Freq")
                lower_tokens.append("freq")
        if "tightscf" not in lower_tokens and "tightscf" not in lower_joined:
            kw.append("TightSCF")
        return kw

    @staticmethod
    def render(species: Artifact, method: dict[str, Any], task: str) -> str:
        xyz_path = Path(str(species.data.get("xyz_path") or species.paths.get("xyz") or species.paths.get("final_xyz"))).resolve()
        charge = int(species.data.get("charge", 0) or 0)
        mult = int(species.data.get("multiplicity", 1) or 1)
        nprocs = int(method.get("nprocs", method.get("ncores", 1)) or 1)
        mem_mb = int(method.get("memory_mb", method.get("maxcore_mb", 2000)) or 2000)
        lines = [
            "# generated_by=hfauto",
            f"# species_id={_canonical_species_id(species)}",
            f"# task={task}",
            "! " + " ".join(OrcaInputRenderer.keywords(method, task)),
            "",
            f"%pal nprocs {nprocs} end",
            f"%maxcore {mem_mb}",
        ]
        for block in method.get("extra_blocks", []) or []:
            lines += ["", str(block).rstrip()]
        lines += ["", f"* xyzfile {charge} {mult} {xyz_path}", ""]
        return "\n".join(lines)


class OrcaEngine:
    """ORCA backend wrapper used by DFTMinimaStage and SinglePointStage."""

    name = "orca"

    def __init__(self, **kwargs: Any):
        self.config = kwargs

    def render_input(self, species: Artifact, method: dict[str, Any], workdir: str | Path, task: str) -> Path:
        """Public renderer used by tests and review tooling.

        It writes the same ORCA input that production execution uses, without
        launching ORCA.
        """
        return self._write_input(species, self._method(method), Path(workdir), task)

    def _method(self, method: dict[str, Any]) -> dict[str, Any]:
        merged = {**self.config, **(method or {})}
        return merged

    def _write_input(self, species: Artifact, method: dict[str, Any], workdir: Path, task: str) -> Path:
        workdir.mkdir(parents=True, exist_ok=True)
        inp = workdir / "orca.inp"
        inp.write_text(OrcaInputRenderer.render(species, method, task), encoding="utf-8")
        return inp

    def _find_final_xyz(self, workdir: Path, output_text: str, species: Artifact) -> Path | None:
        candidates = [workdir / "orca.xyz"] + sorted(workdir.glob("*.xyz"), key=lambda p: p.stat().st_mtime if p.exists() else 0.0, reverse=True)
        for p in candidates:
            if p.exists() and p.is_file():
                try:
                    read_xyz(p)
                    return p
                except Exception:
                    continue
        blocks = parse_orca_cartesian_blocks(output_text)
        if blocks:
            final = workdir / "final_from_output.xyz"
            write_xyz(blocks[-1], final)
            return final
        # For single point calculations, input geometry is the final geometry.
        if species.data.get("xyz_path") and Path(str(species.data["xyz_path"])).exists():
            return Path(str(species.data["xyz_path"]))
        return None

    def _dummy_fallback(self, species: Artifact, method: dict[str, Any], workdir: Path, task: str, reason: str, input_path: Path) -> Artifact:
        dummy = DummyQMEngine().single_point(species, {**method, "task": task}, str(workdir / "dummy_fallback")) if task == "single_point" else DummyQMEngine().optimize_frequency(species, {**method, "task": task}, str(workdir / "dummy_fallback"))
        dummy.artifact_id = "calc_" + fingerprint_dict({"species": species.artifact_id, "method": method, "task": task, "engine": self.name, "fallback": "dummy"})
        dummy.paths = {**dummy.paths, "orca_input": str(input_path)}
        dummy.method = {
            "engine": self.name,
            "method_id": method.get("method_id", "orca"),
            "task": task,
            "backend_fallback": "dummy",
            "intended_engine": "orca",
        }
        dummy.qc = {
            **dummy.qc,
            "fallback_dummy": True,
            "engine_is_dummy": True,
            "real_orca_executed": False,
            "fallback_reason": reason,
            "scientific_use": "software_test_only_not_dft",
        }
        dummy.data = {
            **dummy.data,
            "engine": self.name,
            "requested_engine": "orca",
            "method_id": method.get("method_id", "orca"),
            "task": task,
            "calculation_level": "dummy_fallback_for_orca",
            "real_orca_executed": False,
        }
        return dummy

    def _run(self, species: Artifact, method_in: dict[str, Any], workdir: str, task: str) -> Artifact:
        method = self._method(method_in)
        wd = Path(workdir)
        inp = self._write_input(species, method, wd, task)
        calc_id = "calc_" + fingerprint_dict({"species": species.artifact_id, "method": method, "task": task, "engine": self.name})
        exe = resolve_executable("orca", method.get("executable"))
        fallback_to_dummy = bool(method.get("fallback_to_dummy", False))
        if not allow_orca_subprocess(method) or exe is None or method.get("dry_run", False):
            reason = "ORCA execution disabled" if not allow_orca_subprocess(method) else "ORCA executable not found"
            if method.get("dry_run", False):
                reason = "dry_run=True"
            if fallback_to_dummy:
                return self._dummy_fallback(species, method, wd, task, reason, inp)
            return Artifact.failure(
                artifact_id=calc_id,
                artifact_type="calculation",
                category="orca_not_run",
                reason=f"{reason}. Input was written.",
                parents=[species.artifact_id],
                recoverable=True,
                recommended_fallback="set allow_subprocess=true, configure executable, or enable fallback_to_dummy",
                input=str(inp),
                data={"species_id": _canonical_species_id(species), "task": task, "engine": self.name, "method_id": method.get("method_id", "orca")},
            )

        result = run_command(
            [str(exe), str(inp.name)],
            cwd=wd,
            timeout_s=int(method.get("timeout_s", 86400)),
            env=method.get("env"),
            stdout_name="orca.out",
            stderr_name="orca.err",
        )
        stdout_path = Path(result.stdout_path)
        stderr_path = Path(result.stderr_path)
        output_text = stdout_path.read_text(encoding="utf-8", errors="ignore")
        if stderr_path.exists():
            err = stderr_path.read_text(encoding="utf-8", errors="ignore")
            if err.strip():
                output_text += "\nSTDERR\n" + err
        parsed = parse_orca_output(output_text)
        if (not result.ok) or parsed.get("electronic_energy_hartree") is None:
            fail = Artifact.failure(
                artifact_id=calc_id,
                artifact_type="calculation",
                category="orca_failed",
                reason=f"ORCA returncode={result.returncode}, timed_out={result.timed_out}, parsed_energy={parsed.get('electronic_energy_hartree')}",
                parents=[species.artifact_id],
                recoverable=True,
                recommended_fallback="inspect_orca_output, retry with safer SCF/Opt settings, or lower method",
                input=str(inp),
                output=str(stdout_path),
                stderr=str(stderr_path),
                command_result=asdict(result),
                data={"species_id": _canonical_species_id(species), "task": task, "engine": self.name, "method_id": method.get("method_id", "orca"), **parsed},
            )
            fail.method = {"engine": self.name, "method_id": method.get("method_id", "orca"), "task": task}
            fail.qc = {"real_orca_executed": True, "orca_output_parsed": bool(parsed), "fallback_dummy": False, **{k: parsed.get(k) for k in ["scf_converged", "geometry_converged", "normal_termination", "n_imag"]}}
            return fail

        final_xyz = self._find_final_xyz(wd, output_text, species)
        species_for_qc = dict(species.data)
        if final_xyz is not None:
            species_for_qc["xyz_path"] = str(final_xyz)
        geom_qc = geometry_qc_from_xyz(species_for_qc, final_xyz) if final_xyz is not None else {"geometry_sane": None, "geometry_qc_status": "missing_final_xyz"}
        desc = hf_descriptors_from_species(species_for_qc, parsed)
        data = {
            "calc_id": calc_id,
            "species_id": _canonical_species_id(species),
            "source_species_artifact_id": species.artifact_id,
            "state": species.data.get("state"),
            "task": task,
            "engine": self.name,
            "method_id": method.get("method_id", "orca"),
            "calculation_level": "orca_real",
            "real_orca_executed": True,
            "program_version": parsed.get("program_version"),
            "electronic_energy_hartree": parsed.get("electronic_energy_hartree"),
            "zpe_hartree": parsed.get("zpe_hartree"),
            "enthalpy_298K_hartree": parsed.get("enthalpy_298K_hartree"),
            "gibbs_298K_hartree": parsed.get("gibbs_298K_hartree") or parsed.get("electronic_energy_hartree"),
            "thermal_correction_gibbs_hartree": parsed.get("thermal_correction_gibbs_hartree"),
            "thermal_correction_enthalpy_hartree": parsed.get("thermal_correction_enthalpy_hartree"),
            "n_imag": parsed.get("n_imag"),
            "lowest_freq_cm1": parsed.get("lowest_freq_cm1"),
            "imag_freq_cm1": parsed.get("imag_freq_cm1"),
            "hf_stretch_cm1": parsed.get("hf_stretch_cm1"),
            "dipole_D": parsed.get("dipole_D"),
            **{k: v for k, v in geom_qc.items() if v is not None},
            **{k: v for k, v in desc.items() if v is not None},
        }
        qc = {
            "scf_converged": bool(parsed.get("scf_converged")),
            "geometry_converged": bool(parsed.get("geometry_converged")) if task == "opt_freq" else True,
            "normal_termination": bool(parsed.get("normal_termination")),
            "n_imag": parsed.get("n_imag"),
            "is_minimum": parsed.get("n_imag") == 0 if task == "opt_freq" else None,
            "fallback_dummy": False,
            "engine_is_dummy": False,
            "real_orca_executed": True,
            "orca_output_parsed": True,
            "geometry_sane": geom_qc.get("geometry_sane"),
            "geometry_qc": geom_qc,
        }
        return Artifact(
            artifact_id=calc_id,
            artifact_type="calculation",
            parents=[species.artifact_id],
            paths={
                "input": str(inp),
                "output": str(stdout_path),
                "stderr": str(stderr_path),
                "command_result": str(wd / "command_result.json"),
                "final_xyz": str(final_xyz) if final_xyz is not None else "",
            },
            method={"engine": self.name, "method_id": method.get("method_id", "orca"), "task": task, "keywords": OrcaInputRenderer.keywords(method, task)},
            data=data,
            qc=qc,
            provenance={"command": asdict(result), "created_by": "OrcaEngine"},
        )

    def optimize_frequency(self, species: Artifact, method: dict, workdir: str) -> Artifact:
        return self._run(species, method, workdir, "opt_freq")

    def single_point(self, species: Artifact, method: dict, workdir: str) -> Artifact:
        return self._run(species, method, workdir, "single_point")

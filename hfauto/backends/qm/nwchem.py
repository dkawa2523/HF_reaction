"""Small NWChem backend for minima, frequencies, saddle points, and energies.

The module intentionally owns only three responsibilities: render a readable
NWChem input, execute it, and translate the result to hfauto's calculation
Artifact.  Workflow and ranking policy remain in the stage layer.
"""

from __future__ import annotations

import math
import os
import re
import shutil
from dataclasses import asdict
from pathlib import Path
from typing import Any

from hfauto.backends.qm.dummy import DummyQMEngine
from hfauto.chemistry.descriptors import hf_descriptors_from_species
from hfauto.chemistry.electronic_state import resolve_electronic_state
from hfauto.chemistry.geometry_qc import geometry_qc_from_xyz
from hfauto.chemistry.nwchem_evidence import (
    nwchem_input_geometry_symbols,
    source_file_hashes,
)
from hfauto.chemistry.xyz import XYZ, read_xyz
from hfauto.core.artifacts import canonical_species_id, species_xyz_path
from hfauto.core.executables import CommandResult, resolve_executable, run_command
from hfauto.core.frequency_qc import (
    DEFAULT_IMAGINARY_FREQUENCY_CUTOFF_CM1,
    imaginary_frequency_cutoff_from_method,
    resolve_imaginary_frequency_cutoff,
    significant_imaginary_frequencies,
)
from hfauto.core.hashing import fingerprint_dict, sha256_file
from hfauto.core.schemas.artifact import Artifact

CAL_PER_HARTREE = 627_509.474
THERMOCHEMISTRY_FREQUENCY_CUTOFF_CM1 = 1.0
# Input-renderer/parser semantics; independent of the NWChem program version.
NWCHEM_BACKEND_SCHEMA_VERSION = 1
_BACKEND_SCHEMA_FIELD = "nwchem_backend_schema_version"


def _stamp_backend_schema(artifact: Artifact) -> Artifact:
    """Record the backend semantics used to create a calculation artifact."""

    artifact.data[_BACKEND_SCHEMA_FIELD] = NWCHEM_BACKEND_SCHEMA_VERSION
    artifact.method = {
        **(artifact.method or {}),
        _BACKEND_SCHEMA_FIELD: NWCHEM_BACKEND_SCHEMA_VERSION,
    }
    artifact.provenance[_BACKEND_SCHEMA_FIELD] = NWCHEM_BACKEND_SCHEMA_VERSION
    return artifact


def _execution_provenance(
    result: CommandResult,
    *,
    input_xyz_sha256: str,
    electronic_state: dict[str, int | None],
) -> dict[str, Any]:
    """Return the identical execution evidence for success and failure."""

    return {
        "command": asdict(result),
        "created_by": "NWChemEngine",
        "input_xyz_sha256": input_xyz_sha256,
        "electronic_state": dict(electronic_state),
    }


def allow_nwchem_subprocess(method: dict[str, Any]) -> bool:
    return bool(
        method.get("allow_subprocess", False)
        or os.environ.get("HFAUTO_ALLOW_NWCHEM") == "1"
        or os.environ.get("HFAUTO_ALLOW_SUBPROCESS") == "1"
    )


def _last_float(pattern: str, text: str) -> float | None:
    matches = re.findall(pattern, text, flags=re.IGNORECASE | re.MULTILINE)
    return float(matches[-1]) if matches else None


def _last_population_analysis(text: str, scheme: str) -> list[dict[str, Any]]:
    """Parse the last compact Mulliken/Lowdin atom-population table."""

    header = rf"Total Density - {re.escape(scheme)} Population Analysis"
    starts = list(re.finditer(header, text, flags=re.IGNORECASE))
    if not starts:
        return []
    body = text[starts[-1].end() :]
    boundary = re.search(
        r"Total Density - (?:Mulliken|Lowdin) Population Analysis|"
        r"^\s*center of mass\s*$|\Z",
        body,
        flags=re.IGNORECASE | re.MULTILINE,
    )
    if boundary:
        body = body[: boundary.start()]
    rows: list[dict[str, Any]] = []
    for match in re.finditer(
        r"^\s*(\d+)\s+([A-Za-z]{1,3})\s+(\d+)\s+"
        r"(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s+",
        body,
        flags=re.MULTILINE,
    ):
        atom_index, element, nuclear_charge, electron_population = match.groups()
        z = int(nuclear_charge)
        population = float(electron_population)
        rows.append(
            {
                "atom_index": int(atom_index) - 1,
                "element": element,
                "nuclear_charge": z,
                "electron_population": population,
                "partial_charge_e": float(z - population),
            }
        )
    return rows


def _last_property_mulliken_analysis(text: str) -> list[dict[str, Any]]:
    """Parse Mulliken gross atom populations printed by ``task dft property``."""

    starts = list(
        re.finditer(
            r"-+\s*Total\s+gross population on atoms\s*-+",
            text,
            flags=re.IGNORECASE,
        )
    )
    if not starts:
        return []
    body = text[starts[-1].end() :]
    boundary = re.search(r"^\s*-+\s*Bond indices\s*-+", body, flags=re.MULTILINE)
    if boundary:
        body = body[: boundary.start()]
    rows: list[dict[str, Any]] = []
    for match in re.finditer(
        r"^\s*(\d+)\s+([A-Za-z]{1,3})\s+"
        r"(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s*$",
        body,
        flags=re.MULTILINE,
    ):
        atom_index, element, nuclear_charge, electron_population = match.groups()
        z = round(float(nuclear_charge))
        population = float(electron_population)
        rows.append(
            {
                "atom_index": int(atom_index) - 1,
                "element": element,
                "nuclear_charge": z,
                "electron_population": population,
                "partial_charge_e": float(z - population),
            }
        )
    return rows


def _finite_float(token: str) -> float | None:
    """Parse one NWChem number, rejecting non-finite values."""

    try:
        value = float(token.replace("D", "E").replace("d", "e"))
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def _last_projected_frequency_section(text: str) -> str:
    """Limit frequency decisions to the final geometry's projected analysis."""

    markers = list(
        re.finditer(
            r"^\s*\(Projected Frequencies expressed in cm-1\)\s*$",
            text,
            flags=re.IGNORECASE | re.MULTILINE,
        )
    )
    return text[markers[-1].end() :] if markers else text


def _last_projected_cartesian_modes(
    text: str,
    raw_frequencies: list[float],
    atom_count: int | None,
 ) -> list[list[list[float]]] | None:
    """Read every Cartesian mode from the last projected-frequency section.

    NWChem prints at most six mode columns per block, while every block repeats
    all ``3N`` Cartesian rows.  Any incomplete, duplicate, non-sequential, or
    non-finite block makes the result unusable rather than partially trusted.
    """

    if atom_count is None:
        if not raw_frequencies or len(raw_frequencies) % 3:
            return None
        expected_dimension = len(raw_frequencies)
    else:
        if atom_count <= 0:
            return None
        expected_dimension = 3 * atom_count
    if len(raw_frequencies) != expected_dimension:
        return None

    section_markers = list(
        re.finditer(
            r"^\s*\(Projected Frequencies expressed in cm-1\)\s*$",
            text,
            flags=re.IGNORECASE | re.MULTILINE,
        )
    )
    if not section_markers:
        return None
    lines = text[section_markers[-1].end() :].splitlines()

    frequencies_by_mode: dict[int, float] = {}
    vectors_by_mode: dict[int, list[float]] = {}
    cursor = 0
    while cursor < len(lines):
        stripped = lines[cursor].strip()
        if not stripped:
            cursor += 1
            continue
        if stripped.startswith("-----"):
            break

        header_tokens = stripped.split()
        if not all(token.isdigit() for token in header_tokens):
            return None
        mode_numbers = [int(token) for token in header_tokens]
        next_mode = len(frequencies_by_mode) + 1
        expected_modes = list(range(next_mode, next_mode + len(mode_numbers)))
        if (
            not mode_numbers
            or len(mode_numbers) > 6
            or mode_numbers != expected_modes
            or mode_numbers[-1] > expected_dimension
        ):
            return None

        cursor += 1
        while cursor < len(lines) and not lines[cursor].strip():
            cursor += 1
        if cursor >= len(lines):
            return None
        frequency_tokens = lines[cursor].split()
        if (
            len(frequency_tokens) != len(mode_numbers) + 1
            or frequency_tokens[0].lower() != "p.frequency"
        ):
            return None
        block_frequencies: list[float] = []
        for token in frequency_tokens[1:]:
            value = _finite_float(token)
            if value is None:
                return None
            block_frequencies.append(value)

        cursor += 1
        while cursor < len(lines) and not lines[cursor].strip():
            cursor += 1
        block_vectors = {mode_number: [] for mode_number in mode_numbers}
        for expected_row in range(1, expected_dimension + 1):
            if cursor >= len(lines):
                return None
            row_tokens = lines[cursor].split()
            if len(row_tokens) != len(mode_numbers) + 1:
                return None
            try:
                row_number = int(row_tokens[0])
            except ValueError:
                return None
            if row_number != expected_row:
                return None
            for mode_number, token in zip(mode_numbers, row_tokens[1:]):
                value = _finite_float(token)
                if value is None:
                    return None
                block_vectors[mode_number].append(value)
            cursor += 1

        for mode_number, frequency in zip(mode_numbers, block_frequencies):
            if mode_number in frequencies_by_mode or mode_number in vectors_by_mode:
                return None
            frequencies_by_mode[mode_number] = frequency
            vectors_by_mode[mode_number] = block_vectors[mode_number]

    expected_modes = list(range(1, expected_dimension + 1))
    if list(frequencies_by_mode) != expected_modes:
        return None
    if [frequencies_by_mode[index] for index in expected_modes] != raw_frequencies:
        return None

    modes: list[list[list[float]]] = []
    for mode_number in expected_modes:
        flat_mode = vectors_by_mode[mode_number]
        if len(flat_mode) != expected_dimension:
            return None
        modes.append(
            [
                flat_mode[index : index + 3]
                for index in range(0, expected_dimension, 3)
            ]
        )
    return modes


def _projected_imaginary_modes_from_frequencies(
    text: str,
    raw_frequencies: list[float],
    atom_count: int | None,
    imaginary_frequency_cutoff_cm1: float = (
        DEFAULT_IMAGINARY_FREQUENCY_CUTOFF_CM1
    ),
) -> list[dict[str, Any]] | None:
    """Return every significant mode from one complete projected 3N table."""

    cutoff = resolve_imaginary_frequency_cutoff(
        imaginary_frequency_cutoff_cm1
    )
    imaginary_indices = [
        index
        for index, frequency in enumerate(raw_frequencies)
        if frequency < cutoff
    ]
    if not imaginary_indices:
        return []
    modes = _last_projected_cartesian_modes(text, raw_frequencies, atom_count)
    if modes is None:
        return None
    return [
        {
            "mode_number": mode_index + 1,
            "frequency_cm1": raw_frequencies[mode_index],
            "cartesian_displacements": modes[mode_index],
            "coordinate_convention": (
                "NWChem projected normal-mode eigenvector in Cartesian coordinates"
            ),
            "component_units": "amu^-1/2",
            "mass_weighting_handling": (
                "already transformed by NWChem from the mass-weighted Hessian to "
                "Cartesian coordinates; do not divide by sqrt(atomic_mass) again"
            ),
            "imaginary_frequency_cutoff_cm1": cutoff,
            "raw_frequency_count": len(raw_frequencies),
        }
        for mode_index in imaginary_indices
    ]


def _projected_imaginary_mode_from_frequencies(
    text: str,
    raw_frequencies: list[float],
    atom_count: int | None,
    imaginary_frequency_cutoff_cm1: float = (
        DEFAULT_IMAGINARY_FREQUENCY_CUTOFF_CM1
    ),
) -> dict[str, Any] | None:
    """Return a mode only when the final spectrum has exactly one negative."""

    modes = _projected_imaginary_modes_from_frequencies(
        text,
        raw_frequencies,
        atom_count,
        imaginary_frequency_cutoff_cm1,
    )
    return modes[0] if modes is not None and len(modes) == 1 else None


def extract_nwchem_projected_imaginary_mode(
    text: str,
    atom_count: int,
    imaginary_frequency_cutoff_cm1: float = (
        DEFAULT_IMAGINARY_FREQUENCY_CUTOFF_CM1
    ),
) -> dict[str, Any] | None:
    """Strictly extract one significant mode from NWChem's final 3N table.

    The API deliberately returns ``None`` unless the output contains exactly
    one frequency below the configured cutoff and a complete final projected
    Cartesian eigenvector table whose frequencies exactly match the parsed
    ``P.Frequency`` sequence.  It is therefore safe as an eligibility input,
    but its result alone never establishes a minimum or transition state.
    """

    if isinstance(atom_count, bool) or not isinstance(atom_count, int) or atom_count <= 0:
        return None
    raw_frequencies: list[float] = []
    final_section = _last_projected_frequency_section(text)
    for line in re.findall(
        r"^\s*P\.Frequency\s+(.+)$", final_section, flags=re.MULTILINE
    ):
        tokens = re.findall(
            r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?",
            line,
        )
        values = [_finite_float(token) for token in tokens]
        if not values or any(value is None for value in values):
            return None
        raw_frequencies.extend(float(value) for value in values if value is not None)
    if len(raw_frequencies) != 3 * atom_count:
        return None
    return _projected_imaginary_mode_from_frequencies(
        text,
        raw_frequencies,
        atom_count,
        imaginary_frequency_cutoff_cm1,
    )


def extract_nwchem_projected_imaginary_modes(
    text: str,
    atom_count: int,
    imaginary_frequency_cutoff_cm1: float = (
        DEFAULT_IMAGINARY_FREQUENCY_CUTOFF_CM1
    ),
) -> list[dict[str, Any]] | None:
    """Extract all significant imaginary modes from the final complete table."""

    if isinstance(atom_count, bool) or not isinstance(atom_count, int) or atom_count <= 0:
        return None
    parsed = parse_nwchem_output(
        text,
        atom_count=atom_count,
        imaginary_frequency_cutoff_cm1=imaginary_frequency_cutoff_cm1,
    )
    value = parsed.get("projected_imaginary_modes")
    return list(value) if isinstance(value, list) else None


def parse_nwchem_output(
    text: str,
    temperature_K: float = 298.15,
    atom_count: int | None = None,
    imaginary_frequency_cutoff_cm1: float = (
        DEFAULT_IMAGINARY_FREQUENCY_CUTOFF_CM1
    ),
) -> dict[str, Any]:
    """Parse the compact, stable fields needed by hfauto from NWChem stdout."""
    imaginary_cutoff = resolve_imaginary_frequency_cutoff(
        imaginary_frequency_cutoff_cm1
    )
    energy_pattern = r"Total\s+DFT\s+energy\s*=\s*(-?\d+\.\d+(?:[Ee][+-]?\d+)?)"
    last_reported_energy = _last_float(energy_pattern, text)
    optimization_markers = list(re.finditer(r"Optimization converged", text))
    optimized_energy = (
        _last_float(energy_pattern, text[: optimization_markers[-1].end()])
        if optimization_markers
        else None
    )
    # Numerical frequencies print energies for displaced geometries.  The
    # optimized equilibrium energy is the chemically meaningful endpoint value.
    energy = optimized_energy if optimized_energy is not None else last_reported_energy
    version_match = re.search(r"NWChem\)\s+([0-9][\w.\-]+)", text)
    raw_frequencies: list[float] = []
    final_frequency_section = _last_projected_frequency_section(text)
    for line in re.findall(
        r"^\s*P\.Frequency\s+(.+)$",
        final_frequency_section,
        flags=re.MULTILINE,
    ):
        raw_frequencies.extend(float(value) for value in re.findall(r"-?\d+(?:\.\d+)?", line))
    # Projected translations/rotations are printed as signed zero.  Frequencies
    # below 10 cm-1 remain excluded from the backward-compatible spectral list,
    # but thermochemistry must retain real positive soft modes (for example a
    # +6.83 cm-1 intermolecular mode).  The 1 cm-1 threshold removes only the
    # projected near-zero translations/rotations; imaginary modes remain a
    # separate minimum/TS validation concern and are not fed to qRRHO.
    frequency_analysis_present = bool(raw_frequencies)
    vibrational = [value for value in raw_frequencies if abs(value) >= 10.0]
    thermochemistry_frequencies = [
        value
        for value in raw_frequencies
        if value > THERMOCHEMISTRY_FREQUENCY_CUTOFF_CM1
    ]
    imaginary = significant_imaginary_frequencies(
        raw_frequencies, imaginary_cutoff
    )
    raw_imaginary = [value for value in raw_frequencies if value < 0.0]
    projected_imaginary_modes = _projected_imaginary_modes_from_frequencies(
        text,
        raw_frequencies,
        atom_count,
        imaginary_frequency_cutoff_cm1=imaginary_cutoff,
    )
    projected_imaginary_mode = (
        projected_imaginary_modes[0]
        if projected_imaginary_modes is not None
        and len(projected_imaginary_modes) == 1
        else None
    )
    imaginary_mode_displacements = (
        projected_imaginary_mode["cartesian_displacements"]
        if projected_imaginary_mode is not None
        else None
    )
    mulliken_population = _last_population_analysis(text, "Mulliken")
    if not mulliken_population:
        mulliken_population = _last_property_mulliken_analysis(text)
    lowdin_population = _last_population_analysis(text, "Lowdin")
    zpe = _last_float(r"Zero-Point correction to Energy\s*=.*?\(\s*(-?\d+\.\d+)\s+au\)", text)
    enthalpy_correction = _last_float(r"Thermal correction to Enthalpy\s*=.*?\(\s*(-?\d+\.\d+)\s+au\)", text)
    entropy_cal = _last_float(r"^\s*Total Entropy\s*=\s*(-?\d+\.\d+)\s+cal/mol-K", text)
    reported_temperature = _last_float(
        r"^\s*Temperature\s*=\s*(\d+(?:\.\d+)?)\s*K", text
    )
    temperature_consistent = bool(
        reported_temperature is not None
        and abs(float(reported_temperature) - float(temperature_K)) <= 0.05
    )
    enthalpy = energy + enthalpy_correction if energy is not None and enthalpy_correction is not None else None
    gibbs_correction = None
    gibbs = None
    if (
        enthalpy_correction is not None
        and entropy_cal is not None
        and reported_temperature is not None
    ):
        gibbs_correction = (
            enthalpy_correction
            - float(reported_temperature) * entropy_cal / CAL_PER_HARTREE
        )
        gibbs = energy + gibbs_correction if energy is not None else None
    is_298K = bool(
        reported_temperature is not None
        and abs(float(reported_temperature) - 298.15) <= 0.05
    )
    normal = "Total times" in text and "There is an error in the input file" not in text
    dispersion_correction = _last_float(
        r"Dispersion correction\s*=\s*(-?\d+\.\d+(?:[Ee][+-]?\d+)?)", text
    )
    return {
        "program_version": version_match.group(1) if version_match else None,
        "electronic_energy_hartree": energy,
        "optimized_electronic_energy_hartree": optimized_energy,
        "last_reported_dft_energy_hartree": last_reported_energy,
        "zpe_hartree": zpe,
        "thermochemistry_temperature_K": reported_temperature,
        "requested_temperature_K": float(temperature_K),
        "temperature_consistent": temperature_consistent,
        "enthalpy_hartree": enthalpy,
        "gibbs_hartree": gibbs,
        "enthalpy_298K_hartree": enthalpy if is_298K else None,
        "gibbs_298K_hartree": gibbs if is_298K else None,
        "thermal_correction_enthalpy_hartree": enthalpy_correction,
        "thermal_correction_gibbs_hartree": gibbs_correction,
        "entropy_298K_cal_mol_K": entropy_cal,
        "raw_frequencies_cm1": raw_frequencies,
        "raw_frequency_count": len(raw_frequencies),
        "frequencies_cm1": vibrational,
        "thermochemistry_frequencies_cm1": thermochemistry_frequencies,
        "thermochemistry_frequency_cutoff_cm1": (
            THERMOCHEMISTRY_FREQUENCY_CUTOFF_CM1
        ),
        "frequency_analysis_present": frequency_analysis_present,
        "imaginary_frequency_cutoff_cm1": imaginary_cutoff,
        "n_imag": len(imaginary) if frequency_analysis_present else None,
        "n_imag_raw": len(raw_imaginary) if frequency_analysis_present else None,
        "lowest_freq_cm1": min(vibrational) if vibrational else None,
        "imag_freq_cm1": min(imaginary) if imaginary else None,
        "imaginary_mode_displacements": imaginary_mode_displacements,
        "projected_imaginary_mode": projected_imaginary_mode,
        "projected_imaginary_modes": projected_imaginary_modes,
        "hf_stretch_cm1": max(vibrational) if vibrational else None,
        "mulliken_population_analysis": mulliken_population,
        "lowdin_population_analysis": lowdin_population,
        "population_analysis_present": bool(mulliken_population),
        "scf_converged": energy is not None,
        "geometry_converged": "Optimization converged" in text,
        "normal_termination": normal,
        "dft_d3_applied": "DFT-D3 Model" in text and dispersion_correction is not None,
        "dispersion_correction_hartree": dispersion_correction,
    }


_FREQUENCY_ANALYSIS_FIELDS = (
    "zpe_hartree",
    "thermochemistry_temperature_K",
    "temperature_consistent",
    "enthalpy_hartree",
    "gibbs_hartree",
    "enthalpy_298K_hartree",
    "gibbs_298K_hartree",
    "thermal_correction_enthalpy_hartree",
    "thermal_correction_gibbs_hartree",
    "entropy_298K_cal_mol_K",
    "raw_frequencies_cm1",
    "raw_frequency_count",
    "frequencies_cm1",
    "thermochemistry_frequencies_cm1",
    "frequency_analysis_present",
    "n_imag",
    "n_imag_raw",
    "lowest_freq_cm1",
    "imag_freq_cm1",
    "imaginary_mode_displacements",
    "projected_imaginary_mode",
    "projected_imaginary_modes",
    "hf_stretch_cm1",
)


def _frequency_analysis_snapshot(parsed: dict[str, Any]) -> dict[str, Any]:
    """Copy spectral evidence without presenting it as final-state evidence."""

    return {key: parsed.get(key) for key in _FREQUENCY_ANALYSIS_FIELDS}


def _clear_frequency_analysis(parsed: dict[str, Any]) -> None:
    """Remove final-state claims when no final frequency calculation exists."""

    parsed.update(
        {
            "zpe_hartree": None,
            "thermochemistry_temperature_K": None,
            "temperature_consistent": False,
            "enthalpy_hartree": None,
            "gibbs_hartree": None,
            "enthalpy_298K_hartree": None,
            "gibbs_298K_hartree": None,
            "thermal_correction_enthalpy_hartree": None,
            "thermal_correction_gibbs_hartree": None,
            "entropy_298K_cal_mol_K": None,
            "raw_frequencies_cm1": [],
            "raw_frequency_count": 0,
            "frequencies_cm1": [],
            "thermochemistry_frequencies_cm1": [],
            "frequency_analysis_present": False,
            "n_imag": None,
            "n_imag_raw": None,
            "lowest_freq_cm1": None,
            "imag_freq_cm1": None,
            "imaginary_mode_displacements": None,
            "projected_imaginary_mode": None,
            "projected_imaginary_modes": None,
            "hf_stretch_cm1": None,
        }
    )


def _scope_saddle_frequency_analysis(
    parsed: dict[str, Any],
    text: str,
) -> None:
    """Expose only a frequency analysis performed after saddle convergence.

    A mode-following input may evaluate an initial Hessian before the saddle
    optimization.  If the optimization fails, that Hessian describes the seed,
    not a transition state.  Preserve it as provenance but keep all final-TS
    spectral and thermochemical fields empty.
    """

    convergence_offset = text.rfind("Optimization converged")
    final_section = text[convergence_offset:] if convergence_offset >= 0 else ""
    if "P.Frequency" in final_section:
        return
    if parsed.get("frequency_analysis_present"):
        parsed["preoptimization_hessian_analysis"] = _frequency_analysis_snapshot(
            parsed
        )
    _clear_frequency_analysis(parsed)


def backfill_nwchem_thermochemistry_frequencies(
    calculation: Artifact,
    expected_file_hashes: dict[str, str] | None,
) -> tuple[list[float] | None, list[str]]:
    """Reparse a hash-bound schema-1 output for the new additive mode field."""

    current_hashes, reasons = source_file_hashes(calculation)
    if not isinstance(expected_file_hashes, dict) or current_hashes != expected_file_hashes:
        reasons.append("source_frequency_file_hashes_changed_or_missing")
    output_path = Path(str(calculation.paths.get("output") or ""))
    stderr_path = Path(str(calculation.paths.get("stderr") or ""))
    if reasons or not output_path.is_file() or not stderr_path.is_file():
        if not output_path.is_file():
            reasons.append("source_frequency_output_missing")
        if not stderr_path.is_file():
            reasons.append("source_frequency_stderr_missing")
        return None, list(dict.fromkeys(reasons))
    try:
        text = output_path.read_text(encoding="utf-8", errors="ignore")
        stderr = stderr_path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return None, ["source_frequency_output_unreadable"]
    if stderr.strip():
        text += "\nSTDERR\n" + stderr
    parsed = parse_nwchem_output(
        text,
        temperature_K=float(
            calculation.data.get("requested_temperature_K", 298.15)
        ),
        imaginary_frequency_cutoff_cm1=calculation.data.get(
            "imaginary_frequency_cutoff_cm1",
            DEFAULT_IMAGINARY_FREQUENCY_CUTOFF_CM1,
        ),
    )
    decision_fields = (
        "program_version",
        "electronic_energy_hartree",
        "optimized_electronic_energy_hartree",
        "last_reported_dft_energy_hartree",
        "zpe_hartree",
        "enthalpy_hartree",
        "gibbs_hartree",
        "enthalpy_298K_hartree",
        "gibbs_298K_hartree",
        "thermal_correction_enthalpy_hartree",
        "thermal_correction_gibbs_hartree",
        "entropy_298K_cal_mol_K",
        "thermochemistry_temperature_K",
        "requested_temperature_K",
        "temperature_consistent",
        "raw_frequencies_cm1",
        "raw_frequency_count",
        "frequencies_cm1",
        "frequency_analysis_present",
        "n_imag",
        "n_imag_raw",
        "lowest_freq_cm1",
        "imag_freq_cm1",
        "hf_stretch_cm1",
        "scf_converged",
        "geometry_converged",
        "normal_termination",
        "dft_d3_applied",
        "dispersion_correction_hartree",
    )
    for key in decision_fields:
        if calculation.data.get(key) != parsed.get(key):
            reasons.append(f"source_frequency_raw_{key}_mismatch")
    expected_count = calculation.data.get("expected_raw_frequency_count")
    if (
        calculation.data.get("frequency_count_complete") is not True
        or not isinstance(expected_count, int)
        or isinstance(expected_count, bool)
        or parsed.get("raw_frequency_count") != expected_count
    ):
        reasons.append("source_frequency_count_not_exact")
    if reasons:
        return None, list(dict.fromkeys(reasons))
    return list(parsed["thermochemistry_frequencies_cm1"]), []


class NWChemInputRenderer:
    """Render conventional, reviewable NWChem DFT input decks."""

    @staticmethod
    def _geometry(name: str, xyz: XYZ) -> list[str]:
        lines = [f"geometry {name} units angstrom nocenter noautoz noautosym"] if name else ["geometry units angstrom nocenter noautoz noautosym"]
        lines.extend(f"  {symbol:2s} {x: .10f} {y: .10f} {z: .10f}" for symbol, (x, y, z) in zip(xyz.symbols, xyz.coords))
        lines.append("end")
        return lines

    @staticmethod
    def method_block(method: dict[str, Any], multiplicity: int) -> list[str]:
        functional = str(method.get("functional", method.get("xc", "pbe0")))
        lines = ["dft", f"  xc {functional}", f"  mult {multiplicity}"]
        if multiplicity > 1:
            lines.append("  odft")
        lines += [
            f"  grid {method.get('grid', 'fine')}",
            f"  convergence energy {float(method.get('scf_energy_tolerance', 1.0e-7)):.1e}",
            f"  iterations {int(method.get('scf_maxiter', 150))}",
        ]
        if method.get("disp_vdw") is not None:
            lines.append(f"  disp vdw {int(method['disp_vdw'])}")
        lines.append("end")
        return lines

    @classmethod
    def render(cls, species: Artifact, method: dict[str, Any], task: str) -> str:
        xyz = read_xyz(species_xyz_path(species))
        electronic_state = resolve_electronic_state(species.data, method, xyz.symbols)
        charge = int(electronic_state["charge"])
        multiplicity = int(electronic_state["multiplicity"])
        basis = str(method.get("basis", "def2-svp"))
        lines = [
            "start hfauto_job",
            f'title "hfauto {task} {canonical_species_id(species)}"',
            f"memory total {int(method.get('memory_mb', 2000))} mb",
            "",
            *cls._geometry("", xyz),
            f"charge {charge}",
            "",
            "basis spherical",
            f"  * library {basis}",
            "end",
            "",
            *cls.method_block(method, multiplicity),
        ]
        for block in method.get("extra_blocks", []) or []:
            lines += ["", str(block).rstrip()]
        if task in {"optimize", "opt_freq", "saddle_freq"}:
            frequency_block = [
                "",
                "freq",
                f"  temperature 1 {float(method.get('temperature_K', 298.15)):.2f}",
                "end",
            ]
            initial_saddle_hessian = bool(
                task == "saddle_freq"
                and method.get("saddle_initial_hessian", False)
            )
            if initial_saddle_hessian:
                lines += (
                    ["", "task dft hessian"]
                    if method.get("saddle_initial_hessian_only", False)
                    else [*frequency_block, "task dft frequencies"]
                )
            driver = ["", "driver"]
            thresholds = method.get("driver_thresholds", {}) or {}
            if thresholds:
                for key in ("gmax", "grms", "xmax", "xrms"):
                    if key in thresholds:
                        driver.append(f"  {key} {float(thresholds[key]):.8g}")
            else:
                driver.append(
                    f"  {method.get('optimization_convergence', 'tight')}"
                )
            if method.get("driver_eprec") is not None:
                driver.append(f"  eprec {float(method['driver_eprec']):.8g}")
            if method.get("driver_trust") is not None:
                driver.append(f"  trust {float(method['driver_trust']):.8g}")
            if task == "saddle_freq" and method.get("driver_saddle_step") is not None:
                driver.append(
                    f"  sadstp {float(method['driver_saddle_step']):.8g}"
                )
            if initial_saddle_hessian:
                driver.append("  inhess 2")
                mode_number = method.get("saddle_mode_number")
                if mode_number is not None:
                    if isinstance(mode_number, bool) or int(mode_number) < 1:
                        raise ValueError("saddle_mode_number must be a positive integer")
                    driver.append(f"  moddir {int(mode_number)}")
                driver.append(
                    "  firstneg"
                    if method.get("saddle_follow_first_negative", True)
                    else "  nofirstneg"
                )
            driver += [
                f"  maxiter {int(method.get('geometry_maxiter', 150))}",
                "  xyz final",
                "end",
                "",
                f"task dft {'saddle' if task == 'saddle_freq' else 'optimize'}",
            ]
            lines += driver
            if method.get("population_analysis"):
                lines += ["", "property", "  mulliken", "end", "task dft property"]
            if task in {"opt_freq", "saddle_freq"}:
                lines += [*frequency_block, "task dft frequencies"]
        elif task == "single_point":
            lines += ["", "task dft energy"]
            if method.get("population_analysis"):
                lines += ["", "property", "  mulliken", "end", "task dft property"]
        elif task == "frequency":
            lines += [
                "",
                "freq",
                f"  temperature 1 {float(method.get('temperature_K', 298.15)):.2f}",
                "end",
                "task dft frequencies",
            ]
        else:
            raise ValueError(f"Unsupported NWChem task: {task}")
        return "\n".join(lines) + "\n"


def nwchem_environment(executable: str, configured: dict[str, Any] | None = None) -> dict[str, str]:
    env = {str(k): str(v) for k, v in (configured or {}).items()}
    prefix = Path(executable).resolve().parent.parent
    if (prefix / "share" / "nwchem" / "libraries").exists():
        env.setdefault("NWCHEM_BASIS_LIBRARY", str(prefix / "share" / "nwchem" / "libraries") + "/")
        env.setdefault("NWCHEM_NWPW_LIBRARY", str(prefix / "share" / "nwchem" / "libraryps") + "/")
        env["PATH"] = str(prefix / "bin") + os.pathsep + os.environ.get("PATH", "")
        env["LD_LIBRARY_PATH"] = str(prefix / "lib") + os.pathsep + os.environ.get("LD_LIBRARY_PATH", "")
    return env


def run_nwchem_input(input_path: Path, method: dict[str, Any], stdout_name: str = "nwchem.out") -> tuple[CommandResult, str]:
    executable = resolve_executable("nwchem", method.get("executable"))
    if executable is None:
        raise FileNotFoundError("NWChem executable not found")
    ncores = max(1, int(method.get("ncores", 1)))
    command = [executable, input_path.name]
    if ncores > 1:
        mpi = resolve_executable("mpirun", method.get("mpi_executable"))
        if mpi is None:
            sibling = Path(executable).resolve().parent / "mpirun"
            mpi = str(sibling) if sibling.exists() else None
        if mpi is None:
            raise FileNotFoundError("mpirun is required when ncores > 1")
        command = [mpi, "-np", str(ncores), executable, input_path.name]
    result = run_command(
        command,
        cwd=input_path.parent,
        timeout_s=int(method.get("timeout_s", 86_400)),
        env=nwchem_environment(executable, method.get("env")),
        stdout_name=stdout_name,
        stderr_name=stdout_name.replace(".out", ".err"),
    )
    text = Path(result.stdout_path).read_text(encoding="utf-8", errors="ignore")
    stderr = Path(result.stderr_path).read_text(encoding="utf-8", errors="ignore")
    if stderr.strip():
        text += "\nSTDERR\n" + stderr
    return result, text


class NWChemEngine:
    name = "nwchem"

    def __init__(self, **kwargs: Any):
        self.config = kwargs

    def _method(self, method: dict[str, Any]) -> dict[str, Any]:
        return {**self.config, **(method or {})}

    def render_input(self, species: Artifact, method: dict[str, Any], workdir: str | Path, task: str) -> Path:
        wd = Path(workdir)
        wd.mkdir(parents=True, exist_ok=True)
        path = wd / "nwchem.nw"
        path.write_text(NWChemInputRenderer.render(species, self._method(method), task), encoding="utf-8")
        return path

    @staticmethod
    def _final_xyz(
        workdir: Path,
        species: Artifact,
        input_path: Path,
        task: str,
    ) -> tuple[Path | None, str | None]:
        """Validate and promote the highest numeric NWChem geometry frame."""

        if task in {"single_point", "frequency"}:
            return species_xyz_path(species), None
        try:
            expected_symbols = nwchem_input_geometry_symbols(input_path)
        except (OSError, UnicodeError, ValueError) as exc:
            return None, f"final_xyz_input_geometry_unreadable:{type(exc).__name__}"
        if not expected_symbols:
            return None, "final_xyz_input_geometry_missing"
        candidates: list[tuple[int, str, Path]] = []
        for path in workdir.glob("final-*.xyz"):
            match = re.fullmatch(r"final-(\d+)\.xyz", path.name)
            if match:
                candidates.append((int(match.group(1)), path.name, path))
        if not candidates:
            return None, "final_xyz_numeric_frame_missing"
        _, _, source = max(candidates, key=lambda item: (item[0], item[1]))
        try:
            geometry = read_xyz(source)
        except (OSError, UnicodeError, ValueError) as exc:
            return None, (
                f"final_xyz_highest_numeric_frame_unreadable:{source.name}:"
                f"{type(exc).__name__}"
            )
        if list(geometry.symbols) != expected_symbols:
            return None, (
                "final_xyz_highest_numeric_frame_atom_order_or_count_mismatch:"
                f"{source.name}"
            )
        if any(
            not math.isfinite(float(value))
            for row in geometry.coords
            for value in row
        ):
            return None, (
                f"final_xyz_highest_numeric_frame_nonfinite_coordinates:{source.name}"
            )
        final = workdir / "final.xyz"
        try:
            shutil.copyfile(source, final)
        except OSError as exc:
            return None, f"final_xyz_promotion_failed:{type(exc).__name__}"
        return final, None

    def _fallback(self, species: Artifact, method: dict[str, Any], workdir: Path, task: str, reason: str, input_path: Path) -> Artifact:
        dummy = DummyQMEngine().single_point(species, method, str(workdir / "dummy_fallback")) if task == "single_point" else DummyQMEngine().optimize_frequency(species, method, str(workdir / "dummy_fallback"))
        dummy.artifact_id = "calc_" + fingerprint_dict({"species": species.artifact_id, "method": method, "task": task, "engine": self.name, "fallback": "dummy"})
        dummy.paths["input"] = str(input_path)
        dummy.method = {"engine": self.name, "method_id": method.get("method_id", "nwchem"), "task": task, "backend_fallback": "dummy"}
        dummy.qc.update({"fallback_dummy": True, "real_qm_executed": False, "fallback_reason": reason})
        dummy.data.update({"engine": self.name, "task": task, "real_qm_executed": False, "calculation_level": "dummy_fallback_for_nwchem"})
        return _stamp_backend_schema(dummy)

    def _run(self, species: Artifact, method_in: dict[str, Any], workdir: str, task: str) -> Artifact:
        method = self._method(method_in)
        wd = Path(workdir)
        calc_identity: dict[str, Any] = {
            "species": species.artifact_id,
            "method": method,
            "task": task,
            "engine": self.name,
        }
        try:
            imaginary_cutoff = imaginary_frequency_cutoff_from_method(method)
            source_xyz = species_xyz_path(species)
            source_structure = read_xyz(source_xyz)
            electronic_state = resolve_electronic_state(
                species.data, method, source_structure.symbols
            )
            calc_identity["input_xyz_sha256"] = sha256_file(source_xyz)
            calc_identity["electronic_state"] = electronic_state
            calc_id = "calc_" + fingerprint_dict(calc_identity)
            input_path = self.render_input(species, method, wd, task)
        except (KeyError, OSError, TypeError, ValueError) as exc:
            calc_id = "calc_" + fingerprint_dict(
                {**calc_identity, "input_xyz_sha256": "unavailable"}
            )
            fail = Artifact.failure(
                calc_id,
                "calculation",
                str(exc),
                category="invalid_qm_input",
                parents=[species.artifact_id],
                recoverable=False,
            )
            fail.method = {
                "engine": self.name,
                "method_id": method.get("method_id", "nwchem"),
                "task": task,
            }
            return _stamp_backend_schema(fail)
        executable = resolve_executable("nwchem", method.get("executable"))
        unavailable = not allow_nwchem_subprocess(method) or executable is None or method.get("dry_run", False)
        if unavailable:
            reason = "dry_run=True" if method.get("dry_run", False) else ("NWChem subprocess execution is disabled" if not allow_nwchem_subprocess(method) else "NWChem executable not found")
            if method.get("fallback_to_dummy", False):
                return self._fallback(species, method, wd, task, reason, input_path)
            fail = Artifact.failure(
                calc_id,
                "calculation",
                reason,
                category="nwchem_not_run",
                parents=[species.artifact_id],
                recommended_fallback="configure NWChem or explicitly enable fallback_to_dummy",
                data={
                    "species_id": canonical_species_id(species),
                    "task": task,
                    "engine": self.name,
                    "resolved_charge": electronic_state["charge"],
                    "resolved_multiplicity": electronic_state["multiplicity"],
                    "electron_count": electronic_state["electron_count"],
                },
            )
            fail.paths = {"input": str(input_path)}
            fail.method = {
                "engine": self.name,
                "method_id": method.get("method_id", "nwchem"),
                "task": task,
            }
            return _stamp_backend_schema(fail)
        try:
            result, output = run_nwchem_input(input_path, method)
        except (OSError, TypeError, ValueError) as exc:
            fail = Artifact.failure(
                calc_id,
                "calculation",
                str(exc),
                category="nwchem_failed",
                parents=[species.artifact_id],
            )
            fail.paths = {"input": str(input_path)}
            fail.data.update(
                {
                    "resolved_charge": electronic_state["charge"],
                    "resolved_multiplicity": electronic_state["multiplicity"],
                    "electron_count": electronic_state["electron_count"],
                }
            )
            fail.method = {
                "engine": self.name,
                "method_id": method.get("method_id", "nwchem"),
                "task": task,
            }
            return _stamp_backend_schema(fail)
        atom_count = len(source_structure.symbols)
        parsed = parse_nwchem_output(
            output,
            float(method.get("temperature_K", 298.15)),
            atom_count=atom_count,
            imaginary_frequency_cutoff_cm1=imaginary_cutoff,
        )
        if task == "saddle_freq":
            _scope_saddle_frequency_analysis(parsed, output)
        # A single-point calculation has no geometry convergence or vibrational
        # analysis.  Reporting zero imaginary modes here used to make sparse SP
        # checks look like validated minima to downstream gates.
        if task in {"single_point", "optimize"}:
            parsed.update(
                {
                    "frequencies_cm1": [],
                    "thermochemistry_frequencies_cm1": [],
                    "raw_frequencies_cm1": [],
                    "raw_frequency_count": 0,
                    "frequency_analysis_present": False,
                    "n_imag": None,
                    "n_imag_raw": None,
                    "lowest_freq_cm1": None,
                    "imag_freq_cm1": None,
                    "imaginary_mode_displacements": None,
                    "projected_imaginary_mode": None,
                    "projected_imaginary_modes": None,
                    "hf_stretch_cm1": None,
                }
            )
            if task == "single_point":
                parsed["geometry_converged"] = None
        expected_raw_frequency_count = 3 * atom_count
        frequency_count_complete = bool(
            parsed.get("frequency_analysis_present")
            and int(parsed.get("raw_frequency_count", 0))
            == expected_raw_frequency_count
        )
        parsed["expected_raw_frequency_count"] = expected_raw_frequency_count
        parsed["frequency_count_complete"] = frequency_count_complete
        if task in {"frequency", "opt_freq", "saddle_freq"} and not frequency_count_complete:
            parsed["n_imag"] = None
            parsed["n_imag_raw"] = None
            parsed["imaginary_mode_displacements"] = None
            parsed["projected_imaginary_mode"] = None
            parsed["projected_imaginary_modes"] = None
        final_xyz, final_xyz_error = self._final_xyz(wd, species, input_path, task)
        final_structure = read_xyz(final_xyz) if final_xyz is not None else None
        atom_order_ok = bool(
            final_structure is not None
            and list(final_structure.symbols) == list(source_structure.symbols)
        )
        success = result.ok and parsed["electronic_energy_hartree"] is not None and parsed["normal_termination"]
        required_program_version = method.get("required_program_version")
        if required_program_version is not None:
            success = success and str(parsed.get("program_version")) == str(
                required_program_version
            )
        if method.get("disp_vdw") is not None:
            success = success and parsed["dft_d3_applied"]
        if method.get("population_analysis"):
            success = success and parsed["population_analysis_present"]
        if task in {"optimize", "opt_freq", "saddle_freq"}:
            success = (
                success
                and final_xyz is not None
                and parsed["geometry_converged"]
                and atom_order_ok
            )
        if task in {"frequency", "opt_freq", "saddle_freq"}:
            success = (
                success
                and frequency_count_complete
                and parsed.get("temperature_consistent") is True
            )
        if not success:
            fail = Artifact.failure(
                calc_id,
                "calculation",
                f"NWChem returncode={result.returncode}, energy={parsed['electronic_energy_hartree']}, "
                f"normal={parsed['normal_termination']}, version={parsed.get('program_version')}, "
                f"required_version={required_program_version}, atom_order_ok={atom_order_ok}, "
                f"final_xyz_error={final_xyz_error}",
                category="nwchem_failed",
                parents=[species.artifact_id],
                recommended_fallback="inspect NWChem output and retry with safer SCF or geometry settings",
                data={
                    "species_id": canonical_species_id(species),
                    "task": task,
                    "engine": self.name,
                    "resolved_charge": electronic_state["charge"],
                    "resolved_multiplicity": electronic_state["multiplicity"],
                    "electron_count": electronic_state["electron_count"],
                    "final_xyz_error": final_xyz_error,
                    **parsed,
                },
            )
            fail.paths = {
                "input": str(input_path),
                "output": result.stdout_path,
                "stderr": result.stderr_path,
                "final_xyz": str(final_xyz) if final_xyz else "",
            }
            fail.method = {
                "engine": self.name,
                "method_id": method.get("method_id", "nwchem"),
                "task": task,
                "imaginary_frequency_cutoff_cm1": imaginary_cutoff,
            }
            fail.qc = {"real_qm_executed": True, "nwchem_output_parsed": bool(parsed), "fallback_dummy": False, **{key: parsed.get(key) for key in ["scf_converged", "geometry_converged", "normal_termination", "n_imag"]}}
            fail.provenance = _execution_provenance(
                result,
                input_xyz_sha256=str(calc_identity["input_xyz_sha256"]),
                electronic_state=electronic_state,
            )
            return _stamp_backend_schema(fail)
        species_data = dict(species.data)
        if final_xyz is not None:
            species_data["xyz_path"] = str(final_xyz)
        geometry_qc = geometry_qc_from_xyz(species_data, final_xyz) if final_xyz else {"geometry_sane": None, "geometry_qc_status": "missing_final_xyz"}
        geometry_qc["atom_order_ok"] = atom_order_ok
        geometry_qc["input_symbols"] = list(source_structure.symbols)
        geometry_qc["final_symbols"] = (
            list(final_structure.symbols) if final_structure is not None else []
        )
        descriptors = hf_descriptors_from_species(species_data, parsed)
        data = {
            "calc_id": calc_id,
            "species_id": canonical_species_id(species),
            "source_species_artifact_id": species.artifact_id,
            "state": species.data.get("state"),
            "task": task,
            "engine": self.name,
            "method_id": method.get("method_id", "nwchem"),
            "calculation_level": "nwchem_real",
            "real_qm_executed": True,
            "input_xyz_sha256": calc_identity["input_xyz_sha256"],
            "resolved_charge": electronic_state["charge"],
            "resolved_multiplicity": electronic_state["multiplicity"],
            "electron_count": electronic_state["electron_count"],
            **parsed,
            **{key: value for key, value in geometry_qc.items() if value is not None},
            **{key: value for key, value in descriptors.items() if value is not None},
        }
        qc = {
            "scf_converged": bool(parsed["scf_converged"]),
            "geometry_converged": (
                bool(parsed["geometry_converged"])
                if task not in {"single_point", "frequency"}
                else None
            ),
            "normal_termination": bool(parsed["normal_termination"]),
            "n_imag": parsed["n_imag"],
            "is_minimum": (
                parsed["n_imag"] == 0 and frequency_count_complete
                if task == "opt_freq"
                else None
            ),
            "fallback_dummy": False,
            "engine_is_dummy": False,
            "real_qm_executed": True,
            "nwchem_output_parsed": True,
            "dispersion_requested": method.get("disp_vdw") is not None,
            "dispersion_applied": (
                bool(parsed["dft_d3_applied"]) if method.get("disp_vdw") is not None else None
            ),
            "geometry_sane": geometry_qc.get("geometry_sane"),
            "geometry_qc": geometry_qc,
        }
        artifact = Artifact(
            artifact_id=calc_id,
            artifact_type="calculation",
            parents=[species.artifact_id],
            paths={"input": str(input_path), "output": result.stdout_path, "stderr": result.stderr_path, "final_xyz": str(final_xyz) if final_xyz else ""},
            method={
                "engine": self.name,
                "method_id": method.get("method_id", "nwchem"),
                "task": task,
                "functional": method.get("functional", "pbe0"),
                "basis": method.get("basis", "def2-svp"),
                "disp_vdw": method.get("disp_vdw"),
                "grid": method.get("grid", "fine"),
                "scf_energy_tolerance": method.get("scf_energy_tolerance"),
                "optimization_convergence": method.get(
                    "optimization_convergence"
                ),
                "charge": electronic_state["charge"],
                "multiplicity": electronic_state["multiplicity"],
                "imaginary_frequency_cutoff_cm1": imaginary_cutoff,
            },
            data=data,
            qc=qc,
            provenance=_execution_provenance(
                result,
                input_xyz_sha256=str(calc_identity["input_xyz_sha256"]),
                electronic_state=electronic_state,
            ),
        )
        return _stamp_backend_schema(artifact)

    def optimize_frequency(self, species: Artifact, method: dict[str, Any], workdir: str) -> Artifact:
        return self._run(species, method, workdir, "opt_freq")

    def optimize(self, species: Artifact, method: dict[str, Any], workdir: str) -> Artifact:
        """Optimize geometry without implying a completed frequency analysis."""
        return self._run(species, method, workdir, "optimize")

    def saddle_frequency(self, species: Artifact, method: dict[str, Any], workdir: str) -> Artifact:
        return self._run(species, method, workdir, "saddle_freq")

    def frequency(self, species: Artifact, method: dict[str, Any], workdir: str) -> Artifact:
        """Evaluate a full Hessian at a fixed geometry without optimization."""
        return self._run(species, method, workdir, "frequency")

    def single_point(self, species: Artifact, method: dict[str, Any], workdir: str) -> Artifact:
        return self._run(species, method, workdir, "single_point")

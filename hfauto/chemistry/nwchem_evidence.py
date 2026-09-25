"""Fail-closed NWChem state/method evidence for endpoint calculations.

This module deliberately contains only small, deterministic readers and
comparators.  It does not run NWChem and does not classify chemistry; callers
remain responsible for deciding whether a geometrical endpoint is meaningful.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from hfauto.chemistry.electronic_state import resolve_electronic_state
from hfauto.core.hashing import sha256_file
from hfauto.core.schemas.artifact import Artifact

REQUIRED_ENDPOINT_NWCHEM_METHOD: dict[str, Any] = {
    "engine": "nwchem",
    "program_version": "7.2.3",
    "functional": "pbe0",
    "basis": "def2-svpd",
    "disp_vdw": 3,
}
SOURCE_FILE_KEYS = ("input", "output", "stderr", "final_xyz")
_FATAL_STDERR_MARKERS = (
    "segmentation fault",
    "core dumped",
    "signal 11",
    "mpi_abort",
    "armci abort",
    "ga_error",
    "error termination",
    "nwchem fatal",
    "aborting",
    "fatal error",
    "fortran runtime error",
    "terminate called after throwing",
    "traceback (most recent call last)",
)


def exact_int(value: Any) -> int | None:
    """Accept only an actual integer, never bool, float, or numeric text."""

    if isinstance(value, bool) or not isinstance(value, int):
        return None
    return value


def _single_match(pattern: str, text: str) -> str | None:
    matches = re.findall(pattern, text, flags=re.MULTILINE | re.IGNORECASE)
    return str(matches[0]) if len(matches) == 1 else None


def _single_int(pattern: str, text: str) -> int | None:
    value = _single_match(pattern, text)
    return int(value) if value is not None else None


def _single_float(pattern: str, text: str) -> float | None:
    value = _single_match(pattern, text)
    if value is None:
        return None
    try:
        return float(value.replace("d", "e").replace("D", "E"))
    except ValueError:
        return None


def nwchem_input_geometry_symbols(path: str | Path) -> list[str]:
    """Read the first explicit geometry block without NWChem internals."""

    symbols: list[str] = []
    in_geometry = False
    for line in Path(path).read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = line.strip()
        if not in_geometry:
            if re.match(r"^geometry(?:\s|$)", stripped, flags=re.IGNORECASE):
                in_geometry = True
            continue
        if stripped.lower() == "end":
            break
        fields = stripped.split()
        if len(fields) >= 4 and re.fullmatch(r"[A-Za-z]{1,3}", fields[0]):
            try:
                [float(value) for value in fields[1:4]]
            except ValueError:
                continue
            symbols.append(fields[0].capitalize())
    return symbols


def nwchem_scientific_input_signature(text: str) -> tuple[str, ...]:
    """Normalize an input deck while excluding title and geometry coordinates."""

    signature: list[str] = []
    in_geometry = False
    for raw_line in text.splitlines():
        line = raw_line.strip()
        lowered = line.lower()
        if not line or lowered.startswith("title "):
            continue
        if not in_geometry and re.match(r"^geometry(?:\s|$)", line, re.IGNORECASE):
            signature.extend((lowered, "<geometry-coordinates>"))
            in_geometry = True
            continue
        if in_geometry:
            if lowered == "end":
                signature.append("end")
                in_geometry = False
            continue
        signature.append(lowered)
    if in_geometry:
        signature.append("<unterminated-geometry>")
    return tuple(signature)


def parse_nwchem_input_evidence(path: str | Path) -> dict[str, Any]:
    """Extract the decision-relevant state and method from one input deck.

    A directive is returned as ``None`` when it is absent *or duplicated*.
    This makes ambiguous hand-edited decks fail closed.
    """

    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    functional = _single_match(r"^\s*xc\s+(\S+)\s*$", text)
    basis = _single_match(r"^\s*\*\s+library\s+(\S+)\s*$", text)
    grid = _single_match(r"^\s*grid\s+(\S+)\s*$", text)
    scf_tolerance = _single_float(
        r"^\s*convergence\s+energy\s+(\S+)\s*$", text
    )
    driver_matches = re.findall(
        r"^\s*driver\s*$([\s\S]*?)^\s*end\s*$",
        text,
        flags=re.MULTILINE | re.IGNORECASE,
    )
    driver_text = driver_matches[0] if len(driver_matches) == 1 else ""
    convergence_matches = re.findall(
        r"^\s*(tight|default|loose)\s*$",
        driver_text,
        flags=re.MULTILINE | re.IGNORECASE,
    )
    maxiter_matches = re.findall(
        r"^\s*maxiter\s+(\d+)\s*$",
        driver_text,
        flags=re.MULTILINE | re.IGNORECASE,
    )
    return {
        "geometry_symbols": nwchem_input_geometry_symbols(path),
        "charge": _single_int(r"^\s*charge\s+(-?\d+)\s*$", text),
        "multiplicity": _single_int(r"^\s*mult\s+(\d+)\s*$", text),
        "functional": functional.lower() if functional is not None else None,
        "basis": basis.lower() if basis is not None else None,
        "grid": grid.lower() if grid is not None else None,
        "scf_energy_tolerance": scf_tolerance,
        "disp_vdw": _single_int(r"^\s*disp\s+vdw\s+(\d+)\s*$", text),
        "optimization_convergence": (
            convergence_matches[0].lower()
            if len(convergence_matches) == 1
            else None
        ),
        "geometry_maxiter": (
            int(maxiter_matches[0]) if len(maxiter_matches) == 1 else None
        ),
    }


def paths_equal(first: str | Path | None, second: str | Path | None) -> bool:
    if not first or not second:
        return False
    try:
        return Path(first).resolve(strict=False) == Path(second).resolve(strict=False)
    except OSError:
        return False


def source_file_hashes(calculation: Artifact) -> tuple[dict[str, str], list[str]]:
    """Hash the complete source-file set required for an accepted minimum."""

    hashes: dict[str, str] = {}
    reasons: list[str] = []
    for key in SOURCE_FILE_KEYS:
        value = str((calculation.paths or {}).get(key) or "")
        path = Path(value)
        if not value or not path.is_file():
            reasons.append(f"source_{key}_file_missing")
            continue
        try:
            hashes[key] = sha256_file(path)
        except OSError:
            reasons.append(f"source_{key}_file_unreadable")
    return hashes, reasons


def _state_values(mapping: dict[str, Any] | None) -> dict[str, int | None]:
    values = mapping if isinstance(mapping, dict) else {}
    return {
        "charge": exact_int(values.get("charge")),
        "multiplicity": exact_int(values.get("multiplicity")),
        "electron_count": exact_int(values.get("electron_count")),
    }


def assess_endpoint_nwchem_evidence(
    calculation: Artifact,
    expected_electronic_state: dict[str, Any] | None,
) -> dict[str, Any]:
    """Cross-check raw input, artifact state/method, command, and file evidence."""

    reasons: list[str] = []
    expected = _state_values(expected_electronic_state)
    if any(expected[key] is None for key in expected):
        reasons.append("expected_electronic_state_incomplete")

    input_path = str((calculation.paths or {}).get("input") or "")
    input_evidence: dict[str, Any] = {}
    if not input_path or not Path(input_path).is_file():
        reasons.append("nwchem_input_missing")
    else:
        try:
            input_evidence = parse_nwchem_input_evidence(input_path)
        except (OSError, UnicodeError, ValueError):
            reasons.append("nwchem_input_unreadable")

    symbols = list(input_evidence.get("geometry_symbols") or [])
    if not symbols:
        reasons.append("nwchem_input_geometry_missing")
    input_charge = exact_int(input_evidence.get("charge"))
    input_multiplicity = exact_int(input_evidence.get("multiplicity"))
    if input_charge is None:
        reasons.append("nwchem_input_charge_missing_or_ambiguous")
    if input_multiplicity is None:
        reasons.append("nwchem_input_multiplicity_missing_or_ambiguous")

    computed_state: dict[str, Any] = {}
    if symbols and input_charge is not None and input_multiplicity is not None:
        try:
            computed_state = resolve_electronic_state(
                {"charge": input_charge, "multiplicity": input_multiplicity},
                {},
                symbols,
            )
        except ValueError:
            reasons.append("nwchem_input_electron_parity_invalid")
    if computed_state.get("electron_count") is None:
        reasons.append("nwchem_input_electron_count_unresolved")

    for key in ("charge", "multiplicity", "electron_count"):
        if expected.get(key) is not None and computed_state.get(key) != expected[key]:
            reasons.append(f"nwchem_input_{key}_does_not_match_expected")

    data_state = {
        "charge": exact_int((calculation.data or {}).get("resolved_charge")),
        "multiplicity": exact_int(
            (calculation.data or {}).get("resolved_multiplicity")
        ),
        "electron_count": exact_int((calculation.data or {}).get("electron_count")),
    }
    method_state = {
        "charge": exact_int((calculation.method or {}).get("charge")),
        "multiplicity": exact_int((calculation.method or {}).get("multiplicity")),
    }
    provenance_state = _state_values(
        (calculation.provenance or {}).get("electronic_state")
    )
    for label, observed in (
        ("data", data_state),
        ("method", method_state),
        ("provenance", provenance_state),
    ):
        required_keys = (
            ("charge", "multiplicity")
            if label == "method"
            else ("charge", "multiplicity", "electron_count")
        )
        for key in required_keys:
            if observed.get(key) is None:
                reasons.append(f"artifact_{label}_{key}_missing_or_not_exact_integer")
            elif expected.get(key) is not None and observed[key] != expected[key]:
                reasons.append(f"artifact_{label}_{key}_does_not_match_expected")

    required = REQUIRED_ENDPOINT_NWCHEM_METHOD
    method = calculation.method or {}
    data = calculation.data or {}
    observed_method = {
        "engine": str(method.get("engine") or "").lower(),
        "program_version": str(data.get("program_version") or ""),
        "functional": str(method.get("functional") or "").lower(),
        "basis": str(method.get("basis") or "").lower(),
        "disp_vdw": exact_int(method.get("disp_vdw")),
    }
    for key, required_value in required.items():
        if observed_method.get(key) != required_value:
            reasons.append(f"artifact_method_{key}_is_not_required_value")
    for key in ("functional", "basis", "disp_vdw"):
        if input_evidence.get(key) != required[key]:
            reasons.append(f"nwchem_input_{key}_is_not_required_value")
    if str(data.get("engine") or "").lower() != required["engine"]:
        reasons.append("artifact_data_engine_is_not_nwchem")
    if str(data.get("calculation_level") or "") != "nwchem_real":
        reasons.append("artifact_calculation_level_is_not_nwchem_real")
    if data.get("dft_d3_applied") is not True:
        reasons.append("artifact_data_dft_d3_applied_is_not_true")
    if str((calculation.provenance or {}).get("created_by") or "") != "NWChemEngine":
        reasons.append("artifact_provenance_created_by_is_not_nwchem_engine")

    command = (calculation.provenance or {}).get("command")
    if not isinstance(command, dict):
        reasons.append("source_command_evidence_missing")
        command = {}
    if exact_int(command.get("returncode")) != 0:
        reasons.append("source_command_returncode_is_not_zero")
    if command.get("timed_out") is not False:
        reasons.append("source_command_timed_out_is_not_false")
    for command_key, path_key in (
        ("stdout_path", "output"),
        ("stderr_path", "stderr"),
    ):
        if not paths_equal(command.get(command_key), calculation.paths.get(path_key)):
            reasons.append(f"source_command_{command_key}_path_mismatch")

    hashes, hash_reasons = source_file_hashes(calculation)
    reasons.extend(hash_reasons)
    stderr_path = str((calculation.paths or {}).get("stderr") or "")
    stderr_nonempty = False
    fatal_stderr_markers: list[str] = []
    if stderr_path and Path(stderr_path).is_file():
        try:
            stderr_text = Path(stderr_path).read_text(
                encoding="utf-8", errors="ignore"
            )
            stderr_nonempty = bool(stderr_text.strip())
            lowered_stderr = stderr_text.lower()
            fatal_stderr_markers = [
                marker for marker in _FATAL_STDERR_MARKERS if marker in lowered_stderr
            ]
            if fatal_stderr_markers:
                reasons.append("source_stderr_contains_fatal_marker")
        except OSError:
            reasons.append("source_stderr_file_unreadable")

    return {
        "accepted": not reasons,
        "reasons": list(dict.fromkeys(reasons)),
        "expected_electronic_state": expected,
        "input_electronic_state": computed_state,
        "input_geometry_symbols": symbols,
        "artifact_data_electronic_state": data_state,
        "artifact_method_electronic_state": method_state,
        "artifact_provenance_electronic_state": provenance_state,
        "required_method": dict(required),
        "observed_method": observed_method,
        "input_method": {
            key: input_evidence.get(key)
            for key in ("functional", "basis", "disp_vdw")
        },
        "source_file_sha256": hashes,
        "stderr_nonempty": stderr_nonempty,
        "stderr_fatal_markers": fatal_stderr_markers,
        "command": {
            key: command.get(key)
            for key in ("returncode", "timed_out", "stdout_path", "stderr_path")
        },
    }

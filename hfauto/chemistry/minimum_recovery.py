"""Fail-closed eligibility for recovering a rejected endpoint stationary point."""

from __future__ import annotations

from typing import Any

import numpy as np

from hfauto.chemistry.nwchem_evidence import (
    exact_int,
    parse_nwchem_input_evidence,
    paths_equal,
    source_file_hashes,
)
from hfauto.chemistry.xyz import read_xyz
from hfauto.core.hashing import fingerprint_dict
from hfauto.core.schemas.artifact import Artifact


def assess_minimum_mode_following_source(
    calculation: Artifact,
) -> dict[str, Any]:
    """Accept one real, hash-bound, one-mode rejected minimum attempt.

    The source is a stationary-point calculation, not a minimum.  This gate
    authorizes only construction of two displacement seeds.  It never promotes
    the source or either seed to an endpoint basin.
    """

    reasons: list[str] = []
    if calculation.status.status != "success":
        reasons.append("source_calculation_not_successful")
    if calculation.data.get("task") != "opt_freq":
        reasons.append("source_task_is_not_opt_freq")
    if calculation.qc.get("real_qm_executed") is not True:
        reasons.append("source_real_qm_not_executed")
    if calculation.qc.get("fallback_dummy") is not False:
        reasons.append("source_is_dummy_or_unknown")
    for key in ("scf_converged", "geometry_converged", "normal_termination"):
        if calculation.qc.get(key) is not True:
            reasons.append(f"source_{key}_is_not_true")
    if calculation.data.get("frequency_count_complete") is not True:
        reasons.append("source_frequency_count_not_complete")
    raw_count = exact_int(calculation.data.get("raw_frequency_count"))
    expected_count = exact_int(
        calculation.data.get("expected_raw_frequency_count")
    )
    if raw_count is None or expected_count is None or raw_count != expected_count:
        reasons.append("source_raw_frequency_count_mismatch")
    if exact_int(calculation.data.get("n_imag")) != 1:
        reasons.append("source_does_not_have_exactly_one_imaginary_mode")
    if exact_int(calculation.qc.get("n_imag")) != 1:
        reasons.append("source_qc_does_not_have_exactly_one_imaginary_mode")
    if calculation.qc.get("minimum_accepted") is not False:
        reasons.append("source_is_not_an_explicitly_rejected_minimum")

    final_xyz = str(calculation.paths.get("final_xyz") or "")
    mode_array: np.ndarray | None = None
    geometry_symbols: list[str] = []
    projected_mode = calculation.data.get("projected_imaginary_mode")
    if not isinstance(projected_mode, dict):
        projected_mode = {}
        reasons.append("source_projected_imaginary_mode_missing")
    try:
        geometry = read_xyz(final_xyz)
        geometry_symbols = list(geometry.symbols)
        if expected_count != 3 * len(geometry.symbols):
            reasons.append("source_expected_frequency_count_is_not_three_n")
        mode_array = np.asarray(
            calculation.data.get("imaginary_mode_displacements"), dtype=float
        )
        if mode_array.shape != (len(geometry.symbols), 3):
            reasons.append("source_imaginary_mode_shape_mismatch")
        elif not np.isfinite(mode_array).all() or np.linalg.norm(mode_array) <= 1.0e-14:
            reasons.append("source_imaginary_mode_is_nonfinite_or_zero")
        projected_array = np.asarray(
            projected_mode.get("cartesian_displacements"), dtype=float
        )
        if projected_array.shape != mode_array.shape or not np.allclose(
            projected_array,
            mode_array,
            rtol=0.0,
            atol=1.0e-12,
        ):
            reasons.append("source_projected_mode_displacements_mismatch")
    except (OSError, TypeError, ValueError):
        reasons.append("source_geometry_or_imaginary_mode_unreadable")

    frequency = calculation.data.get("imag_freq_cm1")
    try:
        frequency_value = float(frequency)
        projected_frequency = float(projected_mode.get("frequency_cm1"))
    except (TypeError, ValueError):
        reasons.append("source_imaginary_mode_frequency_unreadable")
    else:
        if not (
            np.isfinite(frequency_value)
            and frequency_value < 0.0
            and np.isclose(
                frequency_value,
                projected_frequency,
                rtol=0.0,
                atol=1.0e-8,
            )
        ):
            reasons.append("source_projected_mode_frequency_mismatch")
    if "cartesian" not in str(
        projected_mode.get("coordinate_convention") or ""
    ).lower():
        reasons.append("source_mode_is_not_explicitly_cartesian")
    if not projected_mode.get("component_units"):
        reasons.append("source_mode_component_units_missing")
    if not projected_mode.get("mass_weighting_handling"):
        reasons.append("source_mode_mass_weighting_handling_missing")

    method = calculation.method or {}
    data = calculation.data or {}
    input_evidence: dict[str, Any] = {}
    if str(method.get("engine") or "").lower() != "nwchem":
        reasons.append("source_mode_provider_is_not_supported")
    else:
        try:
            input_evidence = parse_nwchem_input_evidence(
                calculation.paths.get("input", "")
            )
        except (OSError, TypeError, UnicodeError, ValueError):
            reasons.append("source_nwchem_input_unreadable")
        expected_input: dict[str, Any] = {
            "charge": exact_int(data.get("resolved_charge")),
            "multiplicity": exact_int(data.get("resolved_multiplicity")),
            "functional": str(method.get("functional") or "").lower()
            or None,
            "basis": str(method.get("basis") or "").lower() or None,
            "disp_vdw": exact_int(method.get("disp_vdw")),
            "grid": str(method.get("grid") or "").lower() or None,
            "optimization_convergence": str(
                method.get("optimization_convergence") or ""
            ).lower()
            or None,
            "geometry_maxiter": exact_int(method.get("geometry_maxiter")),
        }
        for key, expected in expected_input.items():
            missing_required = expected is None and key != "disp_vdw"
            if missing_required or input_evidence.get(key) != expected:
                reasons.append(f"source_nwchem_input_{key}_mismatch")
        if expected_input["optimization_convergence"] != "tight":
            reasons.append("source_geometry_optimization_is_not_tight")
        try:
            expected_tolerance = float(method.get("scf_energy_tolerance"))
            input_tolerance = float(
                input_evidence.get("scf_energy_tolerance")
            )
        except (TypeError, ValueError):
            reasons.append("source_nwchem_input_scf_energy_tolerance_mismatch")
        else:
            if not np.isclose(
                input_tolerance,
                expected_tolerance,
                rtol=1.0e-12,
                atol=0.0,
            ):
                reasons.append(
                    "source_nwchem_input_scf_energy_tolerance_mismatch"
                )
        if input_evidence.get("geometry_symbols") != geometry_symbols:
            reasons.append("source_nwchem_input_final_atom_order_mismatch")
        if str(data.get("calculation_level") or "") != "nwchem_real":
            reasons.append("source_calculation_level_is_not_nwchem_real")
        if str(data.get("program_version") or "").strip() == "":
            reasons.append("source_program_version_missing")
        if (exact_int(method.get("disp_vdw")) or 0) > 0 and data.get(
            "dft_d3_applied"
        ) is not True:
            reasons.append("source_requested_dispersion_not_applied")
        if str((calculation.provenance or {}).get("created_by") or "") != (
            "NWChemEngine"
        ):
            reasons.append("source_provider_provenance_mismatch")

    provenance_state = (calculation.provenance or {}).get("electronic_state")
    if not isinstance(provenance_state, dict):
        provenance_state = {}
    for key, data_key in (
        ("charge", "resolved_charge"),
        ("multiplicity", "resolved_multiplicity"),
        ("electron_count", "electron_count"),
    ):
        expected = exact_int(data.get(data_key))
        if expected is None or exact_int(provenance_state.get(key)) != expected:
            reasons.append(f"source_provenance_{key}_mismatch")

    hashes, hash_reasons = source_file_hashes(calculation)
    reasons.extend(hash_reasons)
    command = (calculation.provenance or {}).get("command")
    if not isinstance(command, dict):
        command = {}
        reasons.append("source_command_evidence_missing")
    if exact_int(command.get("returncode")) != 0:
        reasons.append("source_command_returncode_is_not_zero")
    if command.get("timed_out") is not False:
        reasons.append("source_command_timed_out_is_not_false")
    for command_key, path_key in (
        ("stdout_path", "output"),
        ("stderr_path", "stderr"),
    ):
        if not paths_equal(
            command.get(command_key), calculation.paths.get(path_key)
        ):
            reasons.append(f"source_command_{command_key}_path_mismatch")

    return {
        "accepted": not reasons,
        "reasons": list(dict.fromkeys(reasons)),
        "source_calculation_artifact_id": calculation.artifact_id,
        "source_species_id": calculation.data.get("species_id"),
        "imaginary_frequency_cm1": calculation.data.get("imag_freq_cm1"),
        "projected_imaginary_mode": projected_mode,
        "nwchem_input_evidence": input_evidence,
        "source_file_sha256": hashes,
        "imaginary_mode_fingerprint": (
            fingerprint_dict({"cartesian_mode": mode_array.tolist()})
            if mode_array is not None and np.isfinite(mode_array).all()
            else None
        ),
    }

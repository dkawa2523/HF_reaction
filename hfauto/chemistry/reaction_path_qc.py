"""Reaction-path QC helpers for proton-transfer TS/IRC automation.

These helpers are deliberately small and dependency-light.  They do not replace
manual inspection of important candidates; they provide machine-readable gates
that prevent proxy/fallback geometries from being silently treated as validated
activation barriers.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from hfauto.chemistry.basin_identity import match_geometry_to_basin
from hfauto.chemistry.geometry_qc import geometry_qc_from_xyz
from hfauto.chemistry.proton_transfer import (
    hf_endpoint_metrics_from_xyz,
    reaction_coordinate_atoms,
    reaction_spectator_hf_pairs,
)
from hfauto.chemistry.reactions import (
    bond_change_coordinate_value,
    bond_change_coordinate_vector,
    validate_bond_changes_between_geometries,
    validate_reaction_coordinate_between_geometries,
)
from hfauto.chemistry.xyz import XYZ, read_xyz, write_xyz
from hfauto.core.frequency_qc import (
    DEFAULT_IMAGINARY_FREQUENCY_CUTOFF_CM1,
    is_significant_imaginary_frequency,
    resolve_imaginary_frequency_cutoff,
)

# Most abundant isotope masses (amu, NIST AME2016) used to mass-weight modes.
_ISOTOPE_MASS_AMU: dict[str, float] = {
    "H": 1.00782503, "He": 4.00260325, "Li": 7.01600344, "Be": 9.01218307,
    "B": 11.00930536, "C": 12.0, "N": 14.00307400, "O": 15.99491462,
    "F": 18.99840316, "Ne": 19.99244018, "Na": 22.98976928, "Mg": 23.98504170,
    "Al": 26.98153853, "Si": 27.97692653, "P": 30.97376200, "S": 31.97207117,
    "Cl": 34.96885268, "Ar": 39.96238312, "K": 38.96370649, "Ca": 39.96259086,
    "Sc": 44.95590828, "Ti": 47.94794198, "V": 50.94395704, "Cr": 51.94050623,
    "Mn": 54.93804391, "Fe": 55.93493633, "Co": 58.93319429, "Ni": 57.93534241,
    "Cu": 62.92959772, "Zn": 63.92914201, "Ga": 68.92557350, "Ge": 73.92117776,
    "As": 74.92159457, "Se": 79.91652180, "Br": 78.91833760, "Kr": 83.91149773,
    "I": 126.90447190,
}


def _isotope_masses(symbols: list[str]) -> np.ndarray:
    """Return an (n, 1) column of isotope masses; unknown elements raise."""
    unknown = sorted(set(symbols) - _ISOTOPE_MASS_AMU.keys())
    if unknown:
        raise ValueError(f"no isotope mass tabulated for elements {unknown}")
    return np.array([_ISOTOPE_MASS_AMU[symbol] for symbol in symbols]).reshape((-1, 1))


def reaction_coordinate_value(
    species_data: dict[str, Any], xyz_path: str | Path | None
) -> float | None:
    """Return a declared path coordinate, with legacy HF transfer fallback."""
    generic = bond_change_coordinate_value(species_data, xyz_path)
    if generic is not None:
        return generic
    if not xyz_path:
        return None
    qc = geometry_qc_from_xyz(species_data, xyz_path)
    r_bh = qc.get("B_H_distance_A")
    r_hf = qc.get("H_F_distance_A")
    if r_bh is None or r_hf is None:
        return None
    return float(r_hf) - float(r_bh)


def proton_transfer_coordinate_vector(
    species_data: dict[str, Any],
    xyz_path: str | Path | None,
) -> np.ndarray | None:
    """Return the Cartesian gradient of q = r(H-F) - r(B-H).

    The vector is the reference used for a literal normal-mode projection.  It
    has one three-vector per atom and is zero away from the configured B/H/F
    atoms.  ``None`` means the coordinate or geometry could not be resolved.
    """

    if not xyz_path:
        return None
    try:
        xyz = read_xyz(xyz_path)
        atoms = reaction_coordinate_atoms(species_data)
        base = int(atoms["base_atom"])
        proton = int(atoms["transfer_h"])
        fluorine = int(atoms["leaving_f"])
        b_to_h = np.asarray(xyz.coords[proton] - xyz.coords[base], dtype=float)
        h_to_f = np.asarray(xyz.coords[fluorine] - xyz.coords[proton], dtype=float)
    except (FileNotFoundError, IndexError, KeyError, TypeError, ValueError):
        return None
    r_bh = float(np.linalg.norm(b_to_h))
    r_hf = float(np.linalg.norm(h_to_f))
    if r_bh <= 1.0e-12 or r_hf <= 1.0e-12:
        return None
    u_bh = b_to_h / r_bh
    u_hf = h_to_f / r_hf
    gradient = np.zeros_like(np.asarray(xyz.coords, dtype=float))
    gradient[base] = u_bh
    gradient[proton] = -(u_bh + u_hf)
    gradient[fluorine] = u_hf
    return gradient


def reaction_coordinate_vector(
    species_data: dict[str, Any],
    xyz_path: str | Path | None,
) -> np.ndarray | None:
    """Return the declared bond-change vector, falling back to legacy HF PT."""

    generic = bond_change_coordinate_vector(species_data, xyz_path)
    if generic is not None:
        return generic
    return proton_transfer_coordinate_vector(species_data, xyz_path)


def make_displaced_reaction_coordinate_seed(
    reaction_data: dict[str, Any],
    endpoint_xyz: str | Path,
    target_xyz: str | Path,
    out_xyz: str | Path,
    *,
    displacement_A: float = 0.15,
) -> dict[str, Any]:
    """Move an endpoint off its stationary point toward the opposite basin.

    A dimer started exactly at a minimum can satisfy force convergence while
    retaining positive curvature.  The seed therefore follows the local
    gradient of the declared reaction coordinate, with its sign chosen from
    the endpoint-to-target progress.  This is chemistry-declared and does not
    depend on a particular TS backend.
    """

    requested_amplitude = float(displacement_A)
    if not 0.0 < requested_amplitude <= 1.0:
        raise ValueError("dimer seed displacement must be in (0, 1] angstrom")
    endpoint = read_xyz(endpoint_xyz)
    target = read_xyz(target_xyz)
    if endpoint.symbols != target.symbols:
        raise ValueError("endpoint and target atom symbols are not identical")
    direction = reaction_coordinate_vector(reaction_data, endpoint_xyz)
    if (
        direction is None
        or direction.shape != endpoint.coords.shape
        or not np.all(np.isfinite(direction))
    ):
        raise ValueError("declared reaction-coordinate vector is unavailable")

    endpoint_q = reaction_coordinate_value(reaction_data, endpoint_xyz)
    target_q = reaction_coordinate_value(reaction_data, target_xyz)
    if endpoint_q is not None and target_q is not None:
        progress = float(target_q - endpoint_q)
    else:
        progress = float(
            np.sum((target.coords - endpoint.coords) * direction)
        )
    if abs(progress) <= 1.0e-12:
        raise ValueError("reaction-coordinate direction does not distinguish basins")
    sign = 1.0 if progress > 0.0 else -1.0
    seed_path = Path(out_xyz)
    accepted: tuple[float, float | None, dict[str, Any]] | None = None
    for reduction in range(12):
        amplitude = requested_amplitude / (2**reduction)
        seed = XYZ(
            list(endpoint.symbols),
            endpoint.coords + sign * amplitude * direction,
            comment=(
                "state=ts_seed generated_by=hfauto "
                f"reaction_coordinate_displacement_A={amplitude:.8f}"
            ),
        )
        write_xyz(seed, seed_path)
        seed_q = reaction_coordinate_value(reaction_data, seed_path)
        geometry = geometry_qc_from_xyz(reaction_data, seed_path)
        if endpoint_q is not None and target_q is not None and seed_q is not None:
            fraction = (seed_q - endpoint_q) / (target_q - endpoint_q)
            coordinate_progress_ok = 0.0 < fraction <= 0.75
        else:
            coordinate_progress_ok = True
        if coordinate_progress_ok and geometry.get("geometry_sane") is True:
            accepted = (amplitude, seed_q, geometry)
            break
    if accepted is None:
        raise ValueError("could not build a sane dimer seed toward the target basin")
    amplitude, seed_q, geometry = accepted
    return {
        "path": str(seed_path),
        "displacement_A": amplitude,
        "requested_displacement_A": requested_amplitude,
        "direction_sign": sign,
        "endpoint_coordinate": endpoint_q,
        "seed_coordinate": seed_q,
        "target_coordinate": target_q,
        "geometry_qc": geometry,
    }


def make_midpoint_ts_xyz(
    reactant_xyz: str | Path, product_xyz: str | Path, out_xyz: str | Path
) -> Path:
    """Write a midpoint TS guess from atom-matched reactant/product endpoints."""
    rc = read_xyz(reactant_xyz)
    ip = read_xyz(product_xyz)
    if rc.symbols != ip.symbols:
        raise ValueError("reactant/product atom symbols are not identical")
    coords = 0.5 * (rc.coords + ip.coords)
    return write_xyz(
        XYZ(
            list(rc.symbols),
            coords,
            comment="state=transition_state generated_by=hfauto midpoint_ts",
        ),
        out_xyz,
    )


def estimate_reaction_mode_overlap(
    species_data: dict[str, Any],
    xyz_path: str | Path | None,
    n_imag: int | None = None,
    imag_freq_cm1: float | None = None,
    output_text: str | None = None,
    mode_displacements: Any | None = None,
    mode_component_units: str | None = None,
    allow_legacy_geometry_heuristic: bool = False,
    imaginary_frequency_cutoff_cm1: float = (DEFAULT_IMAGINARY_FREQUENCY_CUTOFF_CM1),
) -> tuple[float, str]:
    """Project a literal imaginary-mode displacement onto the reaction coordinate.

    Scientific validation is fail-closed when a Cartesian normal-mode vector is
    unavailable.  The old geometry-only estimate is retained solely behind an
    explicit non-production compatibility switch and never invents the former
    fixed 0.72 score when the coordinate itself is unavailable.
    """

    cutoff = resolve_imaginary_frequency_cutoff(imaginary_frequency_cutoff_cm1)
    if n_imag != 1 or not is_significant_imaginary_frequency(imag_freq_cm1, cutoff):
        return 0.0, "frequency_gate_failed"

    if mode_displacements is not None:
        return _best_mode_reference_overlap(
            species_data,
            xyz_path,
            mode_displacements,
            mode_component_units,
        )

    if not allow_legacy_geometry_heuristic:
        return 0.0, "normal_mode_displacement_unavailable"

    # Compatibility markers and geometry-only scores are explicitly legacy and
    # must not be enabled by production stages.
    if output_text:
        import re

        m = re.search(
            r"HFAUTO_MODE_OVERLAP\s*[:=]\s*(0(?:\.\d+)?|1(?:\.0+)?)",
            output_text,
            flags=re.IGNORECASE,
        )
        if m:
            return float(m.group(1)), "legacy_explicit_hfauto_marker"
    q = reaction_coordinate_value(species_data, xyz_path)
    if q is None:
        return 0.0, "legacy_geometry_coordinate_unavailable"
    # A proton-transfer TS should be near the crossing region.  q exactly zero is
    # ideal; q around +/-0.6 A is still plausible for a loose acid-base TS.
    score = 0.90 - min(abs(float(q)), 1.0) * 0.40
    return max(0.0, min(0.95, score)), "legacy_geometry_coordinate_heuristic"


def _mode_reference_vectors(
    species_data: dict[str, Any],
    xyz_path: str | Path | None,
) -> list[tuple[str, np.ndarray]]:
    """Return explicit coordinate and path-tangent references without guessing."""

    references: list[tuple[str, np.ndarray]] = []
    declared = reaction_coordinate_vector(species_data, xyz_path)
    if declared is not None:
        references.append(("declared_reaction_coordinate", declared))
    configured = species_data.get("reaction_mode_references")
    if not isinstance(configured, list):
        single = species_data.get("reaction_mode_reference")
        configured = [single] if isinstance(single, dict) else []
    for reference in configured:
        if not isinstance(reference, dict):
            continue
        try:
            vector = np.asarray(reference.get("cartesian_displacements"), dtype=float)
        except (TypeError, ValueError):
            continue
        if vector.ndim != 2 or vector.shape[1:] != (3,) or not np.isfinite(vector).all():
            continue
        label = str(reference.get("kind") or "path_tangent")
        references.append((label, vector))
    return references


def _best_mode_reference_overlap(
    species_data: dict[str, Any],
    xyz_path: str | Path | None,
    mode_displacements: Any,
    mode_component_units: str | None,
) -> tuple[float, str]:
    references = _mode_reference_vectors(species_data, xyz_path)
    if not references:
        return 0.0, "reaction_mode_reference_unavailable"
    try:
        displacement = np.asarray(mode_displacements, dtype=float)
    except (TypeError, ValueError):
        return 0.0, "normal_mode_displacement_invalid"
    expected_shape = references[0][1].shape
    if displacement.ndim == 1 and displacement.size == int(np.prod(expected_shape)):
        displacement = displacement.reshape(expected_shape)
    if displacement.shape != expected_shape or not np.isfinite(displacement).all():
        return 0.0, "normal_mode_displacement_invalid"

    masses: np.ndarray | None = None
    projection_prefix = "cartesian"
    if mode_component_units == "amu^-1/2":
        if xyz_path is None:
            return 0.0, "atomic_mass_weighting_unavailable"
        masses = _isotope_masses(read_xyz(xyz_path).symbols)
        if masses.shape[0] != expected_shape[0]:
            return 0.0, "atomic_mass_weighting_unavailable"
        # NWChem prints a Cartesian displacement for each mode.  Move both
        # that displacement and the Cartesian reaction reference into the
        # same mass-weighted space; applying inverse weights to only the
        # reference would change the physical metric.
        displacement = displacement * np.sqrt(masses)
        projection_prefix = "mass_weighted"

    best_score = 0.0
    best_label = "reaction_mode_reference_unavailable"
    for label, raw_reference in references:
        if raw_reference.shape != displacement.shape:
            continue
        reference = raw_reference * np.sqrt(masses) if masses is not None else raw_reference
        denominator = float(np.linalg.norm(reference) * np.linalg.norm(displacement))
        if denominator <= 1.0e-15:
            continue
        score = abs(float(np.vdot(reference.ravel(), displacement.ravel()))) / denominator
        if score > best_score:
            best_score = score
            best_label = (
                f"{projection_prefix}_normal_mode_projection"
                if label == "declared_reaction_coordinate"
                else f"{projection_prefix}_{label}_projection"
            )
    return max(0.0, min(1.0, best_score)), best_label


def select_reaction_coordinate_imaginary_mode(
    species_data: dict[str, Any],
    xyz_path: str | Path | None,
    modes: list[dict[str, Any]] | None,
) -> dict[str, Any] | None:
    """Select one evidenced negative mode by its reaction-reference overlap."""

    candidates: list[dict[str, Any]] = []
    for mode in modes or []:
        if not isinstance(mode, dict):
            continue
        score, method = _best_mode_reference_overlap(
            species_data,
            xyz_path,
            mode.get("cartesian_displacements"),
            str(mode.get("component_units") or "") or None,
        )
        try:
            mode_number = int(mode["mode_number"])
            frequency = float(mode["frequency_cm1"])
        except (KeyError, TypeError, ValueError):
            continue
        candidates.append(
            {
                "mode_number": mode_number,
                "frequency_cm1": frequency,
                "mode_overlap_score": score,
                "mode_overlap_method": method,
            }
        )
    if not candidates:
        return None
    ordered = sorted(
        candidates,
        key=lambda item: item["mode_overlap_score"],
        reverse=True,
    )
    selected = dict(ordered[0])
    runner_up = float(ordered[1]["mode_overlap_score"]) if len(ordered) > 1 else 0.0
    selected.update(
        {
            "runner_up_overlap_score": runner_up,
            "overlap_margin": float(selected["mode_overlap_score"]) - runner_up,
        }
    )
    return selected


def _identity_invariant_endpoint_match(
    observed_xyz: str | Path,
    expected_xyz: str | Path,
    species_data: dict[str, Any],
    *,
    expected_class: str,
    q_tolerance_A: float,
) -> dict[str, Any]:
    """Compare one IRC endpoint to a validated endpoint in proton-state space."""

    try:
        atoms = reaction_coordinate_atoms(species_data)
        spectator_pairs = reaction_spectator_hf_pairs(species_data, expected_xyz)
        expected = hf_endpoint_metrics_from_xyz(read_xyz(expected_xyz), atoms, spectator_pairs)
        observed = hf_endpoint_metrics_from_xyz(read_xyz(observed_xyz), atoms, spectator_pairs)
    except (IndexError, KeyError, OSError, TypeError, ValueError) as exc:
        return {
            "accepted": False,
            "reasons": [f"identity_invariant_endpoint_geometry_unreadable:{exc}"],
            "expected_endpoint_class": expected_class,
        }

    expected_q = float(expected["q_HF_minus_BH_identity_invariant_A"])
    observed_q = float(observed["q_HF_minus_BH_identity_invariant_A"])
    q_delta = abs(observed_q - expected_q)
    reasons: list[str] = []
    if expected.get("endpoint_class") != expected_class:
        reasons.append(f"validated_endpoint_is_{expected.get('endpoint_class')}")
    if expected.get("endpoint_geometry_sane") is not True:
        reasons.append("validated_endpoint_geometry_is_not_sane")
    if observed.get("endpoint_class") != expected_class:
        reasons.append(f"observed_endpoint_is_{observed.get('endpoint_class')}")
    if observed.get("endpoint_geometry_sane") is not True:
        reasons.append("observed_endpoint_geometry_is_not_sane")
    if q_delta > float(q_tolerance_A):
        reasons.append("identity_invariant_q_mismatch")
    return {
        "accepted": not reasons,
        "reasons": reasons,
        "expected_endpoint_class": expected_class,
        "observed_endpoint_class": observed.get("endpoint_class"),
        "validated_endpoint_class": expected.get("endpoint_class"),
        "observed_q_HF_minus_BH_A": observed_q,
        "validated_q_HF_minus_BH_A": expected_q,
        "identity_invariant_q_delta_A": q_delta,
        "identity_invariant_q_tolerance_A": float(q_tolerance_A),
        "observed_endpoint_geometry": observed,
        "validated_endpoint_geometry": expected,
    }


def endpoint_match_qc(
    observed_xyz: str | Path | None,
    expected_xyz: str | Path | None,
    species_data: dict[str, Any],
    rmsd_threshold_A: float = 0.75,
    *,
    expected_class: str,
    q_tolerance_A: float = 0.30,
    require_identity_invariant_geometry: bool = True,
    permutation_rmsd_threshold_A: float = 0.20,
    distance_spectrum_threshold_A: float = 0.08,
) -> dict[str, Any]:
    if not observed_xyz or not expected_xyz:
        return {
            "endpoint_match": False,
            "endpoint_match_status": "missing_endpoint_xyz",
            "endpoint_rmsd_A": None,
        }
    basin_match = match_geometry_to_basin(
        observed_xyz,
        expected_xyz,
        ordered_rmsd_threshold_A=rmsd_threshold_A,
        permutation_rmsd_threshold_A=permutation_rmsd_threshold_A,
        distance_spectrum_threshold_A=distance_spectrum_threshold_A,
    )
    rmsd = basin_match.get("ordered_rmsd_A")
    geometry_match = _identity_invariant_endpoint_match(
        observed_xyz,
        expected_xyz,
        species_data,
        expected_class=expected_class,
        q_tolerance_A=q_tolerance_A,
    )
    rmsd_ok = bool(basin_match["accepted"])
    geometry_ok = bool(geometry_match.get("accepted"))
    accepted = bool(rmsd_ok and (geometry_ok or not require_identity_invariant_geometry))
    reasons: list[str] = []
    if not rmsd_ok:
        reasons.append(str(basin_match["reason"]))
    if require_identity_invariant_geometry:
        reasons.extend(geometry_match.get("reasons", []))
    return {
        "endpoint_match": accepted,
        "endpoint_match_status": "ok" if accepted else reasons[0],
        "endpoint_match_reasons": reasons,
        "endpoint_rmsd_A": rmsd,
        "endpoint_q_observed_A": geometry_match.get("observed_q_HF_minus_BH_A"),
        "endpoint_q_expected_A": geometry_match.get("validated_q_HF_minus_BH_A"),
        "rmsd_threshold_A": rmsd_threshold_A,
        "basin_geometry_match": basin_match,
        "identity_invariant_geometry_required": require_identity_invariant_geometry,
        "identity_invariant_endpoint_match": geometry_match,
    }


def endpoint_pair_match_qc(
    forward_xyz: str | Path | None,
    backward_xyz: str | Path | None,
    reactant_xyz: str | Path | None,
    product_xyz: str | Path | None,
    species_data: dict[str, Any],
    rmsd_threshold_A: float = 0.75,
    q_tolerance_A: float = 0.30,
    require_identity_invariant_geometry: bool = True,
    permutation_rmsd_threshold_A: float = 0.20,
    distance_spectrum_threshold_A: float = 0.08,
) -> dict[str, Any]:
    """Match two IRC endpoint files to reactant/product, accepting either direction."""
    coordinate_terms = (species_data.get("reaction_coordinate") or {}).get("terms") or []
    if species_data.get("bond_changes") or coordinate_terms:

        def basin_match(
            observed: str | Path | None,
            expected: str | Path | None,
        ) -> dict[str, Any]:
            match = match_geometry_to_basin(
                observed,
                expected,
                ordered_rmsd_threshold_A=rmsd_threshold_A,
                permutation_rmsd_threshold_A=permutation_rmsd_threshold_A,
                distance_spectrum_threshold_A=distance_spectrum_threshold_A,
            )
            return {
                "endpoint_match": match["accepted"],
                "endpoint_match_status": match["reason"],
                "endpoint_rmsd_A": match.get("ordered_rmsd_A"),
                "rmsd_threshold_A": rmsd_threshold_A,
                "basin_geometry_match": match,
                "identity_invariant_geometry_required": False,
            }

        forward_product = basin_match(forward_xyz, product_xyz)
        backward_reactant = basin_match(backward_xyz, reactant_xyz)
        forward_reactant = basin_match(forward_xyz, reactant_xyz)
        backward_product = basin_match(backward_xyz, product_xyz)
        direct_ok = bool(forward_product["endpoint_match"] and backward_reactant["endpoint_match"])
        swapped_ok = bool(forward_reactant["endpoint_match"] and backward_product["endpoint_match"])
        if direct_ok or not swapped_ok:
            orientation = "forward_product_backward_reactant"
            chosen_forward, chosen_backward = forward_product, backward_reactant
        else:
            orientation = "forward_reactant_backward_product"
            chosen_forward, chosen_backward = forward_reactant, backward_product
        observed_reactant = backward_xyz if direct_ok or not swapped_ok else forward_xyz
        observed_product = forward_xyz if direct_ok or not swapped_ok else backward_xyz
        if not observed_reactant or not observed_product:
            path_qc = {"accepted": False, "reasons": ["missing_endpoint_xyz"]}
        elif species_data.get("bond_changes"):
            path_qc = validate_bond_changes_between_geometries(
                species_data, observed_reactant, observed_product
            )
        else:
            path_qc = validate_reaction_coordinate_between_geometries(
                species_data, observed_reactant, observed_product
            )
        accepted = bool((direct_ok or swapped_ok) and path_qc["accepted"])
        rmsd_values = [
            value
            for value in (
                chosen_forward.get("endpoint_rmsd_A"),
                chosen_backward.get("endpoint_rmsd_A"),
            )
            if value is not None
        ]
        return {
            "irc_validated": accepted,
            "forward_ok": bool(chosen_forward["endpoint_match"]),
            "backward_ok": bool(chosen_backward["endpoint_match"]),
            "reactant_endpoint_match": accepted,
            "product_endpoint_match": accepted,
            "endpoint_orientation": orientation if accepted else "unmatched",
            "identity_invariant_geometry_required": False,
            "endpoint_validation_scope": (
                "basin_assignment_and_declared_bond_change_graph"
                if species_data.get("bond_changes")
                else "basin_assignment_and_declared_reaction_coordinate"
            ),
            "declared_path_validated": path_qc["accepted"],
            "declared_path_qc": path_qc,
            "bond_change_path_validated": (
                path_qc["accepted"] if species_data.get("bond_changes") else None
            ),
            "bond_change_path_qc": (
                path_qc if species_data.get("bond_changes") else None
            ),
            "endpoint_rmsd_A": max(rmsd_values) if len(rmsd_values) == 2 else None,
            "forward_endpoint_qc": chosen_forward,
            "backward_endpoint_qc": chosen_backward,
        }
    common = {
        "q_tolerance_A": q_tolerance_A,
        "require_identity_invariant_geometry": require_identity_invariant_geometry,
        "permutation_rmsd_threshold_A": permutation_rmsd_threshold_A,
        "distance_spectrum_threshold_A": distance_spectrum_threshold_A,
    }
    f_to_p = endpoint_match_qc(
        forward_xyz,
        product_xyz,
        species_data,
        rmsd_threshold_A,
        expected_class="ion_pair",
        **common,
    )
    b_to_r = endpoint_match_qc(
        backward_xyz,
        reactant_xyz,
        species_data,
        rmsd_threshold_A,
        expected_class="neutral_complex",
        **common,
    )
    f_to_r = endpoint_match_qc(
        forward_xyz,
        reactant_xyz,
        species_data,
        rmsd_threshold_A,
        expected_class="neutral_complex",
        **common,
    )
    b_to_p = endpoint_match_qc(
        backward_xyz,
        product_xyz,
        species_data,
        rmsd_threshold_A,
        expected_class="ion_pair",
        **common,
    )

    direct_ok = bool(f_to_p["endpoint_match"] and b_to_r["endpoint_match"])
    swapped_ok = bool(f_to_r["endpoint_match"] and b_to_p["endpoint_match"])
    if direct_ok or not swapped_ok:
        orientation = "forward_product_backward_reactant"
        chosen_f = f_to_p
        chosen_b = b_to_r
    else:
        orientation = "forward_reactant_backward_product"
        chosen_f = f_to_r
        chosen_b = b_to_p
    ok = bool(direct_ok or swapped_ok)
    return {
        "irc_validated": ok,
        "forward_ok": bool(chosen_f["endpoint_match"]),
        "backward_ok": bool(chosen_b["endpoint_match"]),
        "reactant_endpoint_match": ok,
        "product_endpoint_match": ok,
        "endpoint_orientation": orientation if ok else "unmatched",
        "identity_invariant_geometry_required": require_identity_invariant_geometry,
        "identity_invariant_q_tolerance_A": float(q_tolerance_A),
        "endpoint_rmsd_A": max(
            x
            for x in [chosen_f.get("endpoint_rmsd_A"), chosen_b.get("endpoint_rmsd_A")]
            if x is not None
        )
        if chosen_f.get("endpoint_rmsd_A") is not None
        and chosen_b.get("endpoint_rmsd_A") is not None
        else None,
        "forward_endpoint_qc": chosen_f,
        "backward_endpoint_qc": chosen_b,
    }

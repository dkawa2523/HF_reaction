"""Minimum-basin identity and endpoint assignment.

A stationary geometry is not, by itself, a reaction endpoint identity.  This
module combines frequency evidence, electronic state, geometry, energy, and
the declared reaction coordinate into one method-neutral assessment.  Quantum
chemistry backends only provide evidence; they do not decide basin identity.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from hfauto.chemistry.connectivity import covalent_adjacency
from hfauto.chemistry.geometry import kabsch_rmsd, xyz_file_rmsd
from hfauto.chemistry.reactions import validate_reaction_endpoint_pair
from hfauto.chemistry.xyz import XYZ, read_xyz
from hfauto.core.artifacts import species_xyz_path
from hfauto.core.hashing import fingerprint_dict
from hfauto.core.schemas.artifact import Artifact


def element_pair_distance_signature(xyz: XYZ) -> dict[tuple[str, str], list[float]]:
    """Return a rigid- and same-element-permutation-invariant distance spectrum."""

    signature: defaultdict[tuple[str, str], list[float]] = defaultdict(list)
    for index, symbol in enumerate(xyz.symbols):
        for other_index in range(index):
            pair = tuple(sorted((symbol, xyz.symbols[other_index])))
            delta = xyz.coords[index] - xyz.coords[other_index]
            signature[pair].append(float(math.sqrt(float(delta @ delta))))
    return {pair: sorted(distances) for pair, distances in signature.items()}


def permutation_invariant_distance_rmsd(
    first_path: str | Path,
    second_path: str | Path,
) -> float | None:
    """Compare geometries without depending on same-element atom numbering."""

    try:
        first = read_xyz(first_path)
        second = read_xyz(second_path)
    except (OSError, UnicodeError, ValueError):
        return None
    if Counter(first.symbols) != Counter(second.symbols):
        return None
    first_signature = element_pair_distance_signature(first)
    second_signature = element_pair_distance_signature(second)
    if first_signature.keys() != second_signature.keys():
        return None
    squared_differences: list[float] = []
    for pair, first_distances in first_signature.items():
        second_distances = second_signature[pair]
        if len(first_distances) != len(second_distances):
            return None
        squared_differences.extend(
            (first_value - second_value) ** 2
            for first_value, second_value in zip(first_distances, second_distances)
        )
    return math.sqrt(sum(squared_differences) / len(squared_differences)) if squared_differences else 0.0


def _atom_graph_label(
    xyz: XYZ,
    adjacency: np.ndarray,
    index: int,
) -> tuple[str, int, tuple[tuple[str, int], ...]]:
    neighbours = Counter(
        xyz.symbols[other]
        for other in range(len(xyz.symbols))
        if adjacency[index, other]
    )
    return (
        xyz.symbols[index],
        int(adjacency[index].sum()),
        tuple(sorted(neighbours.items())),
    )


def permutation_invariant_graph_rmsd(
    first_path: str | Path,
    second_path: str | Path,
    *,
    max_isomorphisms: int = 4096,
) -> float | None:
    """Return the best proper-rotation RMSD over covalent graph isomorphisms.

    Unlike a distance spectrum, this distinguishes mirror images.  Atom
    permutations are allowed only when element, degree, neighbour elements,
    and all mapped bond/non-bond relationships agree.
    """

    try:
        first = read_xyz(first_path)
        second = read_xyz(second_path)
    except (OSError, UnicodeError, ValueError):
        return None
    if Counter(first.symbols) != Counter(second.symbols):
        return None
    first_graph = covalent_adjacency(first)
    second_graph = covalent_adjacency(second)
    first_labels = [
        _atom_graph_label(first, first_graph, index)
        for index in range(len(first.symbols))
    ]
    second_labels = [
        _atom_graph_label(second, second_graph, index)
        for index in range(len(second.symbols))
    ]
    if Counter(first_labels) != Counter(second_labels):
        return None
    candidates = {
        index: [
            other
            for other, label in enumerate(second_labels)
            if label == first_labels[index]
        ]
        for index in range(len(first.symbols))
    }
    order = sorted(
        range(len(first.symbols)),
        key=lambda index: (len(candidates[index]), -first_labels[index][1]),
    )
    mapping: dict[int, int] = {}
    used: set[int] = set()
    best: float | None = None
    completed = 0

    def search(position: int) -> None:
        nonlocal best, completed
        if completed >= int(max_isomorphisms):
            return
        if position == len(order):
            completed += 1
            reordered = np.empty_like(first.coords)
            for first_index, second_index in mapping.items():
                reordered[second_index] = first.coords[first_index]
            rmsd = kabsch_rmsd(reordered, second.coords)
            if rmsd is not None and (best is None or rmsd < best):
                best = rmsd
            return
        first_index = order[position]
        for second_index in candidates[first_index]:
            if second_index in used:
                continue
            if any(
                first_graph[first_index, mapped_first]
                != second_graph[second_index, mapped_second]
                for mapped_first, mapped_second in mapping.items()
            ):
                continue
            mapping[first_index] = second_index
            used.add(second_index)
            search(position + 1)
            used.remove(second_index)
            del mapping[first_index]

    search(0)
    return best


def _state(species: Artifact) -> tuple[int, int]:
    return (
        int(species.data.get("charge", species.data.get("resolved_charge", 0)) or 0),
        int(
            species.data.get(
                "multiplicity", species.data.get("resolved_multiplicity", 1)
            )
            or 1
        ),
    )


def _minimum_evidence(
    species: Artifact,
    calculation: Artifact | None,
) -> dict[str, Any]:
    n_imag = None
    frequency_complete = None
    energy = None
    if calculation is not None:
        n_imag = calculation.data.get("n_imag", calculation.qc.get("n_imag"))
        frequency_complete = calculation.data.get(
            "frequency_count_complete", calculation.qc.get("frequency_count_complete")
        )
        energy = calculation.data.get(
            "electronic_energy_hartree",
            calculation.data.get("optimized_electronic_energy_hartree"),
        )
    if n_imag is None:
        n_imag = species.data.get("n_imag", species.qc.get("n_imag"))
    minimum_accepted = bool(
        species.qc.get("minimum_accepted") is True
        or species.qc.get("is_minimum") is True
        or (
            calculation is not None
            and calculation.status.status == "success"
            and calculation.qc.get("is_minimum") is True
        )
    )
    frequency_validated = bool(
        minimum_accepted
        and n_imag == 0
        and frequency_complete is True
    )
    return {
        "species_artifact_id": species.artifact_id,
        "calculation_artifact_id": calculation.artifact_id if calculation else None,
        "minimum_accepted": minimum_accepted,
        "frequency_validated": frequency_validated,
        "n_imag": n_imag,
        "frequency_count_complete": frequency_complete,
        "electronic_energy_hartree": energy,
        "charge": _state(species)[0],
        "multiplicity": _state(species)[1],
        "real_qm_executed": bool(
            species.qc.get("real_qm_executed") is True
            or (
                calculation is not None
                and calculation.qc.get("real_qm_executed") is True
            )
        ),
        "fallback_dummy": bool(
            species.qc.get("fallback_dummy") is True
            or (
                calculation is not None
                and calculation.qc.get("fallback_dummy") is True
            )
        ),
        "method_evidence_validated": bool(
            species.qc.get("state_method_evidence_validated") is True
            or (
                calculation is not None
                and calculation.qc.get("method_evidence_validated") is True
            )
        ),
    }


def _basin_fingerprint(species: Artifact, xyz_path: str | Path) -> str:
    xyz = read_xyz(xyz_path)
    signature = element_pair_distance_signature(xyz)
    return "basin_" + fingerprint_dict(
        {
            "elements": dict(sorted(Counter(xyz.symbols).items())),
            "charge": _state(species)[0],
            "multiplicity": _state(species)[1],
            "distance_spectrum_0p01A": {
                "-".join(pair): [round(value, 2) for value in distances]
                for pair, distances in sorted(signature.items())
            },
        }
    )


def assess_minimum_pair_basin_identity(
    first: Artifact,
    second: Artifact,
    first_calculation: Artifact | None = None,
    second_calculation: Artifact | None = None,
    *,
    same_basin_distance_rmsd_A: float = 0.02,
    distinct_basin_distance_rmsd_A: float = 0.05,
    same_basin_energy_tolerance_hartree: float = 1.0e-5,
) -> dict[str, Any]:
    """Compare two frequency-validated minima without assigning reaction roles.

    A close pair is merged only when invariant geometry and energy both agree.
    A separated pair is a pair of distinct minima, not yet proof of a reaction
    or of a transition-state connection.
    """

    try:
        first_xyz = species_xyz_path(first)
        second_xyz = species_xyz_path(second)
        first_basin_id = _basin_fingerprint(first, first_xyz)
        second_basin_id = _basin_fingerprint(second, second_xyz)
    except (OSError, TypeError, ValueError) as exc:
        return {
            "status": "invalid_endpoints",
            "accepted": False,
            "reasons": [f"endpoint_geometry_unreadable:{exc}"],
        }

    first_evidence = _minimum_evidence(first, first_calculation)
    second_evidence = _minimum_evidence(second, second_calculation)
    state_match = _state(first) == _state(second)
    ordered_rmsd = xyz_file_rmsd(first_xyz, second_xyz)
    graph_rmsd = permutation_invariant_graph_rmsd(first_xyz, second_xyz)
    distance_spectrum_rmsd = permutation_invariant_distance_rmsd(
        first_xyz, second_xyz
    )
    first_energy = first_evidence["electronic_energy_hartree"]
    second_energy = second_evidence["electronic_energy_hartree"]
    energy_difference = (
        abs(float(second_energy) - float(first_energy))
        if first_energy is not None and second_energy is not None
        else None
    )
    both_minima = bool(
        first_evidence["frequency_validated"]
        and second_evidence["frequency_validated"]
    )
    geometry_close = bool(
        graph_rmsd is not None
        and graph_rmsd <= float(same_basin_distance_rmsd_A)
    )
    energy_close = bool(
        energy_difference is not None
        and energy_difference <= float(same_basin_energy_tolerance_hartree)
    )
    geometry_distinct = bool(
        distance_spectrum_rmsd is not None
        and distance_spectrum_rmsd
        >= float(distinct_basin_distance_rmsd_A)
    )

    reasons: list[str] = []
    if not state_match:
        reasons.append("electronic_state_mismatch")
    if not both_minima:
        reasons.append("frequency_validated_minimum_pair_missing")

    if state_match and both_minima and geometry_close and energy_close:
        status = "same_basin"
        accepted = True
        canonical_basin_id = min(first_basin_id, second_basin_id)
        first_basin_id = canonical_basin_id
        second_basin_id = canonical_basin_id
    elif state_match and both_minima and geometry_distinct:
        status = "distinct_basin"
        accepted = True
    else:
        status = "invalid_endpoints" if not state_match or not both_minima else "unresolved"
        accepted = False
        if not geometry_close and not geometry_distinct:
            reasons.append("basin_geometry_separation_ambiguous")
        if geometry_close and energy_difference is None:
            reasons.append("same_basin_energy_evidence_missing")
        elif geometry_close and not energy_close:
            reasons.append("close_geometries_have_distinct_energies")

    return {
        "status": status,
        "accepted": accepted,
        "first_basin_id": first_basin_id,
        "second_basin_id": second_basin_id,
        "same_basin": status == "same_basin",
        "distinct_basin": status == "distinct_basin",
        "electronic_state_match": state_match,
        "ordered_rmsd_A": ordered_rmsd,
        "permutation_invariant_graph_rmsd_A": graph_rmsd,
        "distance_spectrum_rmsd_A": distance_spectrum_rmsd,
        "electronic_energy_difference_hartree": energy_difference,
        "thresholds": {
            "same_basin_distance_rmsd_A": float(same_basin_distance_rmsd_A),
            "distinct_basin_distance_rmsd_A": float(
                distinct_basin_distance_rmsd_A
            ),
            "same_basin_energy_tolerance_hartree": float(
                same_basin_energy_tolerance_hartree
            ),
        },
        "first_evidence": first_evidence,
        "second_evidence": second_evidence,
        "validation_scope": (
            "same_pes_frequency_minima_geometry_and_energy"
        ),
        "reasons": list(dict.fromkeys(reasons)),
    }


def assess_basin_pair(
    reaction: Artifact,
    reactant: Artifact,
    product: Artifact,
    reactant_calculation: Artifact | None = None,
    product_calculation: Artifact | None = None,
    *,
    same_basin_distance_rmsd_A: float = 0.02,
    distinct_basin_distance_rmsd_A: float = 0.05,
    same_basin_energy_tolerance_hartree: float = 1.0e-5,
) -> dict[str, Any]:
    """Assign a minimum pair to declared reactant/product chemistry."""

    comparison = assess_minimum_pair_basin_identity(
        reactant,
        product,
        reactant_calculation,
        product_calculation,
        same_basin_distance_rmsd_A=same_basin_distance_rmsd_A,
        distinct_basin_distance_rmsd_A=distinct_basin_distance_rmsd_A,
        same_basin_energy_tolerance_hartree=(
            same_basin_energy_tolerance_hartree
        ),
    )
    if "first_basin_id" not in comparison:
        return comparison

    endpoint_qc = validate_reaction_endpoint_pair(reaction, reactant, product)
    reasons = list(comparison["reasons"])
    status = str(comparison["status"])
    accepted = bool(comparison["accepted"])
    if endpoint_qc.get("accepted") is not True:
        reasons.extend(endpoint_qc.get("reasons", []))
        if status == "distinct_basin":
            status = "unresolved"
            accepted = False

    return {
        **comparison,
        "status": status,
        "accepted": accepted,
        "reaction_id": reaction.data.get("reaction_id", reaction.artifact_id),
        "reactant_basin_id": comparison["first_basin_id"],
        "product_basin_id": comparison["second_basin_id"],
        "same_basin": status == "same_basin",
        "distinct_basin": status == "distinct_basin",
        "reactant_evidence": comparison["first_evidence"],
        "product_evidence": comparison["second_evidence"],
        "endpoint_chemistry_qc": endpoint_qc,
        "validation_scope": (
            "same_pes_frequency_minima_geometry_energy_and_declared_chemistry"
        ),
        "reasons": list(dict.fromkeys(reasons)),
    }


def match_geometry_to_basin(
    observed_xyz: str | Path | None,
    reference_xyz: str | Path | None,
    *,
    ordered_rmsd_threshold_A: float = 0.75,
    permutation_rmsd_threshold_A: float = 0.20,
    distance_spectrum_threshold_A: float = 0.08,
) -> dict[str, Any]:
    """Assign a trajectory endpoint to one known basin geometry.

    Ordered RMSD handles atom-mapped paths.  Covalent-graph isomorphisms provide
    a symmetry-tolerant, chirality-preserving fallback.  The distance spectrum
    is diagnostic only because it cannot distinguish mirror images.
    """

    if not observed_xyz or not reference_xyz:
        return {
            "accepted": False,
            "reason": "missing_endpoint_xyz",
            "ordered_rmsd_A": None,
            "distance_spectrum_rmsd_A": None,
        }
    ordered = xyz_file_rmsd(observed_xyz, reference_xyz)
    graph_rmsd = permutation_invariant_graph_rmsd(observed_xyz, reference_xyz)
    distance_spectrum = permutation_invariant_distance_rmsd(
        observed_xyz, reference_xyz
    )
    ordered_ok = ordered is not None and ordered <= float(ordered_rmsd_threshold_A)
    permutation_ok = (
        graph_rmsd is not None
        and graph_rmsd <= float(permutation_rmsd_threshold_A)
    )
    accepted = bool(ordered_ok or permutation_ok)
    return {
        "accepted": accepted,
        "reason": "ok" if accepted else "geometry_outside_reference_basin",
        "ordered_rmsd_A": ordered,
        "ordered_rmsd_threshold_A": float(ordered_rmsd_threshold_A),
        "permutation_invariant_graph_rmsd_A": graph_rmsd,
        "permutation_rmsd_threshold_A": float(permutation_rmsd_threshold_A),
        "distance_spectrum_rmsd_A": distance_spectrum,
        "distance_spectrum_threshold_A": float(distance_spectrum_threshold_A),
        "match_method": (
            "ordered_rmsd"
            if ordered_ok
            else "covalent_graph_isomorphism_rmsd"
            if permutation_ok
            else None
        ),
    }

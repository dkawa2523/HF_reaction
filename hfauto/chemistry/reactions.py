"""Chemistry-neutral reaction endpoint and mode-coordinate helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from hfauto.chemistry.geometry import align_coordinates, xyz_file_rmsd
from hfauto.chemistry.xyz import read_xyz
from hfauto.core.artifacts import species_xyz_path
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.chemistry import BondChangeRecord, ReactionCoordinateTermRecord

PROTON_TRANSFER_TYPES = {"proton_transfer", "hf_proton_transfer"}


def is_proton_transfer_reaction(reaction: Artifact | dict[str, Any]) -> bool:
    data = reaction.data if isinstance(reaction, Artifact) else reaction
    family = str(data.get("mechanism_family") or data.get("reaction_type") or "")
    if family:
        normalized = family.lower()
        if normalized not in PROTON_TRANSFER_TYPES:
            return False
        if normalized == "hf_proton_transfer":
            return True
        atoms = (data.get("reaction_coordinate", {}) or {}).get("atoms", {}) or {}
        if data.get("bond_changes") or data.get("hypothesis"):
            return {"base_atom", "transfer_h", "leaving_f"} <= set(atoms)
        return True
    # Artifacts produced before reaction families were introduced were all
    # proton-transfer reactions. Keep that compatibility at this one boundary
    # instead of spreading legacy checks through TS and IRC stages.
    return not bool(data.get("bond_changes") or data.get("hypothesis"))


def is_proton_transfer_state(species: Artifact | dict[str, Any]) -> bool:
    data = species.data if isinstance(species, Artifact) else species
    atoms = (data.get("reaction_coordinate", {}) or {}).get("atoms", {}) or {}
    if {"base_atom", "transfer_h", "leaving_f"} <= set(atoms):
        return True
    if data.get("bond_changes") or data.get("chemical_state"):
        return False
    # ``reactant_complex`` and ``ion_pair`` were historically HF-specific.
    # New generic states carry ``chemical_state`` and therefore do not enter
    # this compatibility branch.
    return str(data.get("state", "")) in {"reactant_complex", "ion_pair"}


def _bond_changes(data: dict[str, Any]) -> list[BondChangeRecord]:
    return [
        item if isinstance(item, BondChangeRecord) else BondChangeRecord.model_validate(item)
        for item in (data.get("bond_changes", []) or [])
    ]


def reaction_coordinate_terms(
    data: dict[str, Any],
) -> list[ReactionCoordinateTermRecord]:
    """Return explicit terms, or the coordinate implied by bond changes.

    Keeping this normalization public gives path initialization, endpoint
    validation, and mode validation one shared interpretation of a reaction
    coordinate.
    """
    coordinate = data.get("reaction_coordinate") or {}
    if hasattr(coordinate, "model_dump"):
        coordinate = coordinate.model_dump()
    explicit = list(coordinate.get("terms", []) or [])
    if explicit:
        return [ReactionCoordinateTermRecord.model_validate(term) for term in explicit]
    terms: list[ReactionCoordinateTermRecord] = []
    for change in _bond_changes(data):
        sign = -1.0 if change.kind in {"form", "increase_order"} else 1.0
        terms.append(
            ReactionCoordinateTermRecord(
                kind="distance",
                atoms=change.atoms,
                coefficient=sign * change.weight,
                label=change.label,
            )
        )
    return terms


def coordinate_term_value(
    term: ReactionCoordinateTermRecord, coordinates: np.ndarray
) -> float:
    """Evaluate one distance, angle, or periodic dihedral coordinate."""
    points = coordinates[list(term.atoms)]
    if term.kind == "distance":
        return float(np.linalg.norm(points[0] - points[1]))
    first = points[0] - points[1]
    second = points[2] - points[1]
    if term.kind == "angle":
        denominator = float(np.linalg.norm(first) * np.linalg.norm(second))
        if denominator <= 1.0e-15:
            raise ValueError("angle coordinate contains a zero-length bond")
        cosine = float(np.clip(np.dot(first, second) / denominator, -1.0, 1.0))
        return float(np.arccos(cosine))
    first_normal = np.cross(points[1] - points[0], points[2] - points[1])
    second_normal = np.cross(points[2] - points[1], points[3] - points[2])
    axis = points[2] - points[1]
    normalizer = float(
        np.linalg.norm(first_normal) * np.linalg.norm(second_normal) * np.linalg.norm(axis)
    )
    if normalizer <= 1.0e-15:
        raise ValueError("dihedral coordinate contains collinear or coincident atoms")
    x = float(np.dot(first_normal, second_normal))
    y = float(np.dot(np.cross(first_normal, second_normal), axis / np.linalg.norm(axis)))
    return float(np.arctan2(y, x))


def _term_gradient(
    term: ReactionCoordinateTermRecord,
    coordinates: np.ndarray,
    step_A: float = 1.0e-5,
) -> np.ndarray:
    """Numerical internal-coordinate gradient, including periodic dihedrals."""

    gradient = np.zeros_like(coordinates, dtype=float)
    for atom_index in term.atoms:
        for axis in range(3):
            plus = coordinates.copy()
            minus = coordinates.copy()
            plus[atom_index, axis] += step_A
            minus[atom_index, axis] -= step_A
            delta = coordinate_term_value(term, plus) - coordinate_term_value(
                term, minus
            )
            if term.kind == "dihedral":
                delta = float(np.arctan2(np.sin(delta), np.cos(delta)))
            gradient[atom_index, axis] = delta / (2.0 * step_A)
    return gradient


def validate_reaction_endpoint_pair(
    reaction: Artifact | dict[str, Any],
    reactant: Artifact,
    product: Artifact,
) -> dict[str, Any]:
    """Validate generic NEB endpoints from atom order and declared bond changes.

    No element-specific distance cutoff is guessed. Each change is validated by
    its direction and configured minimum distance change, which keeps the rule
    usable for molecular, ionic, and surface-bound systems.
    """

    data = reaction.data if isinstance(reaction, Artifact) else reaction
    reasons: list[str] = []
    try:
        reactant_path = species_xyz_path(reactant)
        product_path = species_xyz_path(product)
        reactant_xyz = read_xyz(reactant_path)
        product_xyz = read_xyz(product_path)
    except (OSError, TypeError, ValueError) as exc:
        return {"accepted": False, "reasons": [f"endpoint_geometry_unreadable:{exc}"]}

    atom_order_ok = reactant_xyz.symbols == product_xyz.symbols
    if not atom_order_ok:
        reasons.append("endpoint_atom_order_mismatch")
    if int(reactant.data.get("charge", 0)) != int(product.data.get("charge", 0)):
        reasons.append("endpoint_charge_mismatch")
    if int(reactant.data.get("multiplicity", 1)) != int(product.data.get("multiplicity", 1)):
        reasons.append("endpoint_multiplicity_mismatch")
    reactant_order_key = reactant.data.get("atom_order_key")
    product_order_key = product.data.get("atom_order_key")
    if reactant_order_key and product_order_key and reactant_order_key != product_order_key:
        reasons.append("endpoint_atom_mapping_mismatch")

    bond_change_qc = validate_bond_changes_between_geometries(
        data, reactant_path, product_path
    )
    coordinate_qc = validate_reaction_coordinate_between_geometries(
        data, reactant_path, product_path
    )
    derived_coordinate = endpoint_derived_coordinate_basis(
        data, reactant_path, product_path
    )
    if data.get("bond_changes"):
        reasons.extend(bond_change_qc["reasons"])
        validation_scope = "atom_order_charge_spin_and_declared_bond_changes"
    else:
        reasons.extend(coordinate_qc["reasons"])
        validation_scope = "atom_order_charge_spin_and_declared_reaction_coordinate"

    rmsd = xyz_file_rmsd(reactant_path, product_path) if atom_order_ok else None
    if rmsd is not None and rmsd < float(data.get("minimum_endpoint_rmsd_A", 0.02)):
        reasons.append("endpoints_are_geometrically_indistinguishable")
    return {
        "accepted": not reasons,
        "reasons": list(dict.fromkeys(reasons)),
        "validation_scope": validation_scope,
        "atom_order_ok": atom_order_ok,
        "endpoint_rmsd_A": rmsd,
        "bond_change_evidence": bond_change_qc.get("bond_change_evidence", []),
        "reaction_coordinate_evidence": coordinate_qc,
        "endpoint_derived_coordinate": derived_coordinate,
    }


def endpoint_derived_coordinate_basis(
    reaction_data: dict[str, Any],
    reactant_xyz_path: str | Path,
    product_xyz_path: str | Path,
    *,
    minimum_distance_change_A: float = 0.05,
    maximum_terms: int = 12,
    maximum_atoms: int = 256,
) -> dict[str, Any]:
    """Audit a declared coordinate against the full endpoint displacement.

    The ranked distance changes are hypotheses for path initialization, not
    inferred bonds.  Automatic bond relabeling from Cartesian distances would
    be unsafe for ionic, coordination, and weakly bound molecular complexes.
    """

    try:
        reactant = read_xyz(reactant_xyz_path)
        product = read_xyz(product_xyz_path)
        terms = reaction_coordinate_terms(reaction_data)
    except (OSError, TypeError, ValueError) as exc:
        return {
            "accepted": False,
            "reasons": [f"geometry_unreadable:{exc}"],
        }
    if reactant.symbols != product.symbols:
        return {"accepted": False, "reasons": ["atom_order_mismatch"]}
    atom_count = len(reactant.symbols)
    if atom_count > int(maximum_atoms):
        return {
            "accepted": False,
            "reasons": ["endpoint_coordinate_audit_atom_budget_exceeded"],
            "atom_count": atom_count,
            "maximum_atoms": int(maximum_atoms),
        }
    changes: list[dict[str, Any]] = []
    for first in range(atom_count):
        for second in range(first):
            reactant_distance = float(
                np.linalg.norm(
                    reactant.coords[first] - reactant.coords[second]
                )
            )
            product_distance = float(
                np.linalg.norm(product.coords[first] - product.coords[second])
            )
            delta = product_distance - reactant_distance
            if abs(delta) < float(minimum_distance_change_A):
                continue
            changes.append(
                {
                    "kind": "distance",
                    "atoms": [second, first],
                    "symbols": [
                        reactant.symbols[second],
                        reactant.symbols[first],
                    ],
                    "reactant_value_A": reactant_distance,
                    "product_value_A": product_distance,
                    "delta_A": delta,
                    "coefficient_for_positive_progress": (
                        1.0 if delta > 0.0 else -1.0
                    ),
                    "interpretation": (
                        "distance_increase"
                        if delta > 0.0
                        else "distance_decrease"
                    ),
                }
            )
    changes.sort(key=lambda item: (-abs(item["delta_A"]), item["atoms"]))
    ranked = changes[: max(0, int(maximum_terms))]

    try:
        aligned_product = align_coordinates(
            reactant.coords, product.coords
        )
        displacement = (
            np.asarray(aligned_product, dtype=float)
            - np.asarray(reactant.coords, dtype=float)
        ).reshape(-1)
        midpoint = 0.5 * (
            np.asarray(reactant.coords, dtype=float)
            + np.asarray(aligned_product, dtype=float)
        )
        gradients = [
            _term_gradient(term, midpoint).reshape(-1) for term in terms
        ]
    except (TypeError, ValueError):
        gradients = []
        displacement = np.asarray([], dtype=float)

    displacement_norm = float(np.linalg.norm(displacement))
    nonzero_gradients = [
        gradient
        for gradient in gradients
        if float(np.linalg.norm(gradient)) > 1.0e-12
    ]
    span_overlap: float | None = None
    direction_cosine: float | None = None
    if displacement_norm > 1.0e-12 and nonzero_gradients:
        normalized = np.column_stack(
            [gradient / np.linalg.norm(gradient) for gradient in nonzero_gradients]
        )
        basis, _ = np.linalg.qr(normalized)
        projection = basis @ (basis.T @ displacement)
        span_overlap = float(np.linalg.norm(projection) / displacement_norm)
        combined = sum(
            (
                float(term.coefficient) * gradient
                for term, gradient in zip(terms, gradients)
            ),
            start=np.zeros_like(displacement),
        )
        combined_norm = float(np.linalg.norm(combined))
        if combined_norm > 1.0e-12:
            direction_cosine = float(
                np.dot(combined, displacement)
                / (combined_norm * displacement_norm)
            )

    declared_distance_pairs = {
        tuple(sorted(term.atoms))
        for term in terms
        if term.kind == "distance"
    }
    ranked_pairs = [tuple(item["atoms"]) for item in ranked]
    return {
        "accepted": True,
        "reasons": [],
        "definition": (
            "ranked_atom_mapped_endpoint_distance_changes_after_independent_"
            "minimum_validation"
        ),
        "automatic_bond_claims_allowed": False,
        "minimum_distance_change_A": float(minimum_distance_change_A),
        "significant_distance_change_count": len(changes),
        "recommended_coordinate_terms": ranked,
        "declared_term_count": len(terms),
        "declared_coordinate_span_overlap": span_overlap,
        "declared_coordinate_direction_cosine": direction_cosine,
        "primary_distance_change_declared": bool(
            ranked_pairs and ranked_pairs[0] in declared_distance_pairs
        ),
    }


def validate_reaction_coordinate_between_geometries(
    reaction_data: dict[str, Any],
    reactant_xyz_path: str | Path,
    product_xyz_path: str | Path,
) -> dict[str, Any]:
    """Validate progress along an explicit, weighted internal coordinate.

    This is the bond-change-free path for conformational reactions.  The term
    coefficients define the positive reaction direction; dihedral differences
    are wrapped to ``[-pi, pi]`` before they are combined.
    """

    try:
        reactant = read_xyz(reactant_xyz_path)
        product = read_xyz(product_xyz_path)
        terms = reaction_coordinate_terms(reaction_data)
        coordinate = reaction_data.get("reaction_coordinate") or {}
        if hasattr(coordinate, "model_dump"):
            coordinate = coordinate.model_dump()
        min_change = float(coordinate.get("min_change", 0.05))
    except (OSError, TypeError, ValueError) as exc:
        return {"accepted": False, "reasons": [f"geometry_unreadable:{exc}"]}
    if reactant.symbols != product.symbols:
        return {"accepted": False, "reasons": ["atom_order_mismatch"]}
    if not terms:
        return {"accepted": False, "reasons": ["reaction_coordinate_terms_missing"]}

    evidence: list[dict[str, Any]] = []
    progress = 0.0
    for term in terms:
        if max(term.atoms) >= len(reactant.symbols):
            return {
                "accepted": False,
                "reasons": [f"reaction_coordinate_atom_out_of_range:{max(term.atoms)}"],
            }
        try:
            reactant_value = coordinate_term_value(term, reactant.coords)
            product_value = coordinate_term_value(term, product.coords)
        except ValueError as exc:
            return {
                "accepted": False,
                "reasons": [f"reaction_coordinate_unreadable:{exc}"],
            }
        delta = product_value - reactant_value
        if term.kind == "dihedral":
            delta = float(np.arctan2(np.sin(delta), np.cos(delta)))
        contribution = float(term.coefficient) * delta
        progress += contribution
        evidence.append(
            {
                **term.model_dump(),
                "reactant_value": reactant_value,
                "product_value": product_value,
                "wrapped_delta": delta,
                "progress_contribution": contribution,
            }
        )
    accepted = progress >= min_change
    return {
        "accepted": accepted,
        "reasons": [] if accepted else ["reaction_coordinate_progress_not_observed"],
        "progress": progress,
        "min_change": min_change,
        "term_evidence": evidence,
    }


def validate_bond_changes_between_geometries(
    reaction_data: dict[str, Any],
    reactant_xyz_path: str | Path,
    product_xyz_path: str | Path,
) -> dict[str, Any]:
    """Check declared bond-change direction on any atom-mapped geometry pair."""

    reasons: list[str] = []
    try:
        reactant = read_xyz(reactant_xyz_path)
        product = read_xyz(product_xyz_path)
        changes = _bond_changes(reaction_data)
    except (OSError, TypeError, ValueError) as exc:
        return {"accepted": False, "reasons": [f"geometry_unreadable:{exc}"]}
    if reactant.symbols != product.symbols:
        return {"accepted": False, "reasons": ["atom_order_mismatch"]}
    if not changes:
        return {"accepted": False, "reasons": ["bond_changes_missing"]}
    evidence: list[dict[str, Any]] = []
    for change in changes:
        first, second = change.atoms
        if max(first, second) >= len(reactant.symbols):
            reasons.append(f"bond_change_atom_out_of_range:{first}-{second}")
            continue
        reactant_distance = float(
            np.linalg.norm(reactant.coords[first] - reactant.coords[second])
        )
        product_distance = float(
            np.linalg.norm(product.coords[first] - product.coords[second])
        )
        delta = product_distance - reactant_distance
        if change.kind in {"form", "increase_order"}:
            direction_ok = delta <= -change.min_distance_change_A
        elif change.kind in {"break", "decrease_order"}:
            direction_ok = delta >= change.min_distance_change_A
        else:
            direction_ok = abs(delta) >= change.min_distance_change_A
        if not direction_ok:
            reasons.append(f"bond_change_not_observed:{change.label or f'{first}-{second}'}")
        evidence.append(
            {
                **change.model_dump(),
                "reactant_distance_A": reactant_distance,
                "product_distance_A": product_distance,
                "distance_change_A": delta,
                "direction_ok": direction_ok,
            }
        )
    return {
        "accepted": not reasons,
        "reasons": list(dict.fromkeys(reasons)),
        "bond_change_evidence": evidence,
    }


def bond_change_coordinate_value(
    species_data: dict[str, Any],
    xyz_path: str | Path | None,
) -> float | None:
    """Return one distance-based progress coordinate for declared bond changes.

    Forming bonds use negative distance and breaking bonds positive distance, so
    the coordinate increases from reactant to product for either change type.
    The mean keeps the scale stable when a hypothesis contains several changes.
    """

    if not xyz_path:
        return None
    try:
        xyz = read_xyz(xyz_path)
        terms = reaction_coordinate_terms(species_data)
    except (OSError, TypeError, ValueError):
        return None
    values: list[float] = []
    for term in terms:
        if max(term.atoms) >= len(xyz.coords):
            return None
        try:
            values.append(
                term.coefficient * coordinate_term_value(term, xyz.coords)
            )
        except ValueError:
            return None
    return float(sum(values)) if values else None


def bond_change_coordinate_vector(
    species_data: dict[str, Any],
    xyz_path: str | Path | None,
) -> np.ndarray | None:
    """Return a normalized path direction from declared bond changes.

    Forming bonds contribute in the distance-shortening direction; breaking
    bonds contribute in the distance-lengthening direction.
    """

    if not xyz_path:
        return None
    try:
        xyz = read_xyz(xyz_path)
        terms = reaction_coordinate_terms(species_data)
    except (OSError, TypeError, ValueError):
        return None
    if not terms:
        return None
    vector = np.zeros_like(np.asarray(xyz.coords, dtype=float))
    for term in terms:
        if max(term.atoms) >= len(xyz.coords):
            return None
        try:
            vector += term.coefficient * _term_gradient(term, xyz.coords)
        except ValueError:
            return None
    norm = float(np.linalg.norm(vector))
    return vector / norm if norm > 1.0e-12 else None

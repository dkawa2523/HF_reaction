"""Chemistry-aware initial paths for double-ended reaction searches.

This module only constructs geometry paths.  It neither selects an electronic
structure method nor runs NEB/string.  Unsupported coordinates return
``None`` so method backends can choose another documented strategy without
pretending that a chemically informed path was built.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import pairwise
from typing import Any

import numpy as np

from hfauto.chemistry.connectivity import (
    covalent_adjacency,
    covalent_radius_A,
    fragments_after_bond_cut,
)
from hfauto.chemistry.geometry import (
    align_coordinates,
    align_coordinates_by_indices,
    kabsch_rmsd,
)
from hfauto.chemistry.reactions import (
    coordinate_term_value,
    reaction_coordinate_terms,
)
from hfauto.chemistry.xyz import XYZ
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.chemistry import (
    BondChangeRecord,
    ReactionCoordinateTermRecord,
)


@dataclass(frozen=True)
class InitializedReactionPath:
    """A generated trajectory and the evidence needed to audit it."""

    images: list[XYZ]
    metadata: dict[str, Any]


def _wrapped_angle(value: float) -> float:
    return float(np.arctan2(np.sin(value), np.cos(value)))


def _rotate_fragment(
    coordinates: np.ndarray,
    atom_indices: list[int],
    axis_start: np.ndarray,
    axis_end: np.ndarray,
    angle: float,
) -> np.ndarray:
    """Rotate one covalent fragment around a bond with Rodrigues' formula."""

    axis = np.asarray(axis_end, dtype=float) - np.asarray(axis_start, dtype=float)
    norm = float(np.linalg.norm(axis))
    if norm <= 1.0e-12:
        raise ValueError("torsional path has a zero-length central bond")
    unit = axis / norm
    result = np.asarray(coordinates, dtype=float).copy()
    shifted = result[atom_indices] - axis_start
    cosine, sine = float(np.cos(angle)), float(np.sin(angle))
    rotated = (
        shifted * cosine
        + np.cross(unit, shifted) * sine
        + np.outer(shifted @ unit, unit) * (1.0 - cosine)
    )
    result[atom_indices] = rotated + axis_start
    return result


def _rotation_choices(
    term: ReactionCoordinateTermRecord,
    start: XYZ,
    aligned_end: np.ndarray,
    target_delta: float,
) -> list[tuple[list[int], float, np.ndarray, float]]:
    first, center_left, center_right, last = term.atoms
    fragments = fragments_after_bond_cut(start, center_left, center_right)
    if fragments is None:
        raise ValueError(
            "dihedral central bond is absent or cyclic; rigid-fragment "
            "initialization is not valid"
        )
    left, right = fragments
    if first not in left or last not in right:
        raise ValueError(
            "dihedral atom order does not span the two central-bond fragments"
        )

    target_value = coordinate_term_value(term, aligned_end)
    candidates: list[tuple[float, float, int, list[int], float, np.ndarray]] = []
    for fragment in (left, right):
        for trial_angle in (target_delta, -target_delta):
            rotated = _rotate_fragment(
                start.coords,
                fragment,
                start.coords[center_left],
                start.coords[center_right],
                trial_angle,
            )
            coordinate_error = abs(
                _wrapped_angle(coordinate_term_value(term, rotated) - target_value)
            )
            residual_rmsd = float(
                np.sqrt(np.mean(np.sum((aligned_end - rotated) ** 2, axis=1)))
            )
            candidates.append(
                (
                    coordinate_error,
                    residual_rmsd,
                    len(fragment),
                    fragment,
                    trial_angle,
                    rotated,
                )
            )
    valid = [candidate for candidate in candidates if candidate[0] <= 1.0e-5]
    if not valid:
        raise ValueError("rigid rotation cannot reproduce the endpoint dihedral")
    return [
        (fragment, angle, rotated_end, residual_rmsd)
        for (
            _coordinate_error,
            residual_rmsd,
            _size,
            fragment,
            angle,
            rotated_end,
        ) in sorted(valid, key=lambda item: item[:3])
    ]


def _path_geometry_qc(images: list[XYZ], adjacency: np.ndarray) -> dict[str, float]:
    nonbonded_distances: list[float] = []
    maximum_bond_stretch_ratio = 1.0
    minimum_bond_compression_ratio = 1.0
    segment_rmsds: list[float] = []
    endpoint_bond_ranges: dict[tuple[int, int], tuple[float, float]] = {}
    for first in range(len(images[0].symbols)):
        for second in range(first):
            if adjacency[first, second]:
                endpoint_lengths = [
                    float(
                        np.linalg.norm(
                            image.coords[first] - image.coords[second]
                        )
                    )
                    for image in (images[0], images[-1])
                ]
                endpoint_bond_ranges[(first, second)] = (
                    min(endpoint_lengths),
                    max(endpoint_lengths),
                )
    for image in images:
        for first in range(len(image.symbols)):
            for second in range(first):
                distance = float(
                    np.linalg.norm(image.coords[first] - image.coords[second])
                )
                if adjacency[first, second]:
                    lower, upper = endpoint_bond_ranges[(first, second)]
                    compression = distance / lower
                    stretch = distance / upper
                    minimum_bond_compression_ratio = min(
                        minimum_bond_compression_ratio, compression
                    )
                    maximum_bond_stretch_ratio = max(
                        maximum_bond_stretch_ratio, stretch
                    )
                else:
                    nonbonded_distances.append(
                        distance
                    )
    for first, second in pairwise(images):
        rmsd = kabsch_rmsd(first.coords, second.coords)
        if rmsd is not None:
            segment_rmsds.append(rmsd)
    minimum_distance = min(nonbonded_distances, default=float("inf"))
    if minimum_distance < 0.55:
        raise ValueError(
            "generated torsional path contains a nonbonded atom collision "
            f"({minimum_distance:.3f} angstrom)"
        )
    if minimum_bond_compression_ratio < 0.70:
        raise ValueError(
            "generated path collapses a covalent bond "
            f"(minimum endpoint-relative ratio {minimum_bond_compression_ratio:.3f})"
        )
    if maximum_bond_stretch_ratio > 1.50:
        raise ValueError(
            "generated path over-stretches a covalent bond "
            f"(maximum endpoint-relative ratio {maximum_bond_stretch_ratio:.3f})"
        )
    return {
        "minimum_nonbonded_distance_A": minimum_distance,
        "maximum_aligned_segment_rmsd_A": max(segment_rmsds, default=0.0),
        "minimum_covalent_bond_compression_ratio": minimum_bond_compression_ratio,
        "maximum_covalent_bond_stretch_ratio": maximum_bond_stretch_ratio,
    }


def _bonded_transfer_pattern(
    data: dict[str, Any], start: XYZ
) -> tuple[int, int, int] | None:
    """Return ``(moving_atom, donor, acceptor)`` for one unambiguous transfer."""

    changes = [
        BondChangeRecord.model_validate(change)
        for change in (data.get("bond_changes") or [])
    ]
    forming = [
        change for change in changes if change.kind in {"form", "increase_order"}
    ]
    breaking = [
        change for change in changes if change.kind in {"break", "decrease_order"}
    ]
    patterns: list[tuple[int, int, int]] = []
    adjacency = covalent_adjacency(start)
    for formed in forming:
        for broken in breaking:
            shared = set(formed.atoms) & set(broken.atoms)
            if len(shared) != 1:
                continue
            moving = shared.pop()
            acceptor = next(atom for atom in formed.atoms if atom != moving)
            donor = next(atom for atom in broken.atoms if atom != moving)
            if (
                adjacency[moving, donor]
                and not adjacency[moving, acceptor]
                and adjacency[donor, acceptor]
            ):
                patterns.append((moving, donor, acceptor))
    unique = list(dict.fromkeys(patterns))
    return unique[0] if len(unique) == 1 else None


def _deterministic_perpendicular(axis: np.ndarray) -> np.ndarray:
    basis = np.eye(3)[int(np.argmin(np.abs(axis)))]
    perpendicular = np.cross(axis, basis)
    return perpendicular / np.linalg.norm(perpendicular)


def _transfer_arc_images(
    start: XYZ,
    aligned_end: np.ndarray,
    *,
    moving: int,
    donor: int,
    acceptor: int,
    image_count: int,
    direction: float,
) -> list[XYZ]:
    axis_vector = start.coords[acceptor] - start.coords[donor]
    axis_length = float(np.linalg.norm(axis_vector))
    if axis_length <= 1.0e-12:
        raise ValueError("atom-transfer centers are coincident")
    axis = axis_vector / axis_length
    perpendicular = _deterministic_perpendicular(axis) * float(direction)
    start_center = 0.5 * (start.coords[donor] + start.coords[acceptor])
    end_center = 0.5 * (aligned_end[donor] + aligned_end[acceptor])
    start_radius = float(np.linalg.norm(start.coords[moving] - start_center))
    end_radius = float(np.linalg.norm(aligned_end[moving] - end_center))
    target_contact = 1.25 * max(
        covalent_radius_A(start.symbols[moving])
        + covalent_radius_A(start.symbols[donor]),
        covalent_radius_A(start.symbols[moving])
        + covalent_radius_A(start.symbols[acceptor]),
    )
    middle_radius = max(
        0.80,
        float(
            np.sqrt(
                max(target_contact * target_contact - (0.5 * axis_length) ** 2, 0.0)
            )
        ),
    )
    images: list[XYZ] = []
    for index, fraction in enumerate(np.linspace(0.0, 1.0, image_count)):
        value = float(fraction)
        smooth = value * value * (3.0 - 2.0 * value)
        coordinates = start.coords + smooth * (aligned_end - start.coords)
        center = 0.5 * (coordinates[donor] + coordinates[acceptor])
        linear_radius = (1.0 - value) * start_radius + value * end_radius
        radius = linear_radius - np.sin(np.pi * value) * (
            linear_radius - middle_radius
        )
        angle = np.pi * value
        arc_direction = -np.cos(angle) * axis + np.sin(angle) * perpendicular
        ideal_start = start_center - start_radius * axis
        ideal_end = end_center + end_radius * axis
        endpoint_residual = (
            (1.0 - value) * (start.coords[moving] - ideal_start)
            + value * (aligned_end[moving] - ideal_end)
        )
        coordinates[moving] = center + radius * arc_direction + endpoint_residual
        if index == 0:
            coordinates = start.coords.copy()
        elif index == image_count - 1:
            coordinates = aligned_end.copy()
        images.append(
            XYZ(
                symbols=list(start.symbols),
                coords=coordinates,
                comment=(
                    "initial_path=bonded_center_transfer_arc "
                    f"image={index + 1}/{image_count} fraction={value:.8f}"
                ),
            )
        )
    return images


def _initialize_bonded_atom_transfer(
    data: dict[str, Any], start: XYZ, end: XYZ, image_count: int
) -> InitializedReactionPath | None:
    pattern = _bonded_transfer_pattern(data, start)
    if pattern is None:
        return None
    moving, donor, acceptor = pattern
    scaffold = [index for index in range(len(start.symbols)) if index != moving]
    aligned_end = align_coordinates_by_indices(
        start.coords, end.coords, scaffold
    )
    adjacency = covalent_adjacency(start)
    candidates: list[tuple[float, list[XYZ], dict[str, float]]] = []
    errors: list[ValueError] = []
    for direction in (-1.0, 1.0):
        images = _transfer_arc_images(
            start,
            aligned_end,
            moving=moving,
            donor=donor,
            acceptor=acceptor,
            image_count=image_count,
            direction=direction,
        )
        try:
            qc = _path_geometry_qc(images, adjacency)
        except ValueError as exc:
            errors.append(exc)
            continue
        candidates.append((qc["minimum_nonbonded_distance_A"], images, qc))
    if not candidates:
        raise errors[0] if errors else ValueError("no collision-free transfer arc")
    _clearance, images, geometry_qc = max(candidates, key=lambda item: item[0])
    return InitializedReactionPath(
        images=images,
        metadata={
            "strategy": "bonded_center_transfer_arc",
            "moving_atom_index": moving,
            "donor_atom_index": donor,
            "acceptor_atom_index": acceptor,
            "image_count": image_count,
            **geometry_qc,
        },
    )


def initialize_reaction_path(
    reaction: Artifact | dict[str, Any],
    start: XYZ,
    end: XYZ,
    *,
    image_count: int,
) -> InitializedReactionPath | None:
    """Build a smooth path for one unambiguous acyclic dihedral change.

    The moving covalent fragment remains rigid during the torsional part.  A
    smooth endpoint-residual interpolation then accommodates the small bond
    and angle relaxation that distinguishes independently optimized minima.
    This gives NEB/string a chemically meaningful initial topology while the
    electronic-structure optimizer remains responsible for the final MEP.
    """

    count = int(image_count)
    if count < 3:
        raise ValueError("an initialized reaction path requires at least three images")
    if list(start.symbols) != list(end.symbols):
        raise ValueError("reaction path endpoints have different atom order")
    data = reaction.data if isinstance(reaction, Artifact) else reaction
    terms = reaction_coordinate_terms(dict(data or {}))
    if len(terms) != 1 or terms[0].kind != "dihedral":
        return _initialize_bonded_atom_transfer(
            dict(data or {}), start, end, count
        )

    term = terms[0]
    if max(term.atoms) >= len(start.symbols):
        raise ValueError("reaction-coordinate atom index is outside the geometry")
    aligned_end = align_coordinates(start.coords, end.coords)
    start_value = coordinate_term_value(term, start.coords)
    end_value = coordinate_term_value(term, aligned_end)
    target_delta = _wrapped_angle(end_value - start_value)
    if abs(target_delta) < 1.0e-6:
        return None

    candidates: list[
        tuple[float, list[XYZ], list[int], float, dict[str, float]]
    ] = []
    errors: list[ValueError] = []
    for moving_atoms, rotation_angle, rotated_end, residual_rmsd in _rotation_choices(
        term, start, aligned_end, target_delta
    ):
        residual = aligned_end - rotated_end
        images: list[XYZ] = []
        for index, fraction in enumerate(np.linspace(0.0, 1.0, count)):
            rotated = _rotate_fragment(
                start.coords,
                moving_atoms,
                start.coords[term.atoms[1]],
                start.coords[term.atoms[2]],
                rotation_angle * float(fraction),
            )
            smooth_fraction = float(
                fraction * fraction * (3.0 - 2.0 * fraction)
            )
            coordinates = rotated + smooth_fraction * residual
            if index == 0:
                coordinates = start.coords.copy()
            elif index == count - 1:
                coordinates = aligned_end.copy()
            images.append(
                XYZ(
                    symbols=list(start.symbols),
                    coords=coordinates,
                    comment=(
                        "initial_path=acyclic_dihedral_rotation "
                        f"image={index + 1}/{count} fraction={fraction:.8f}"
                    ),
                )
            )
        try:
            geometry_qc = _path_geometry_qc(
                images, covalent_adjacency(start)
            )
        except ValueError as exc:
            errors.append(exc)
            continue
        candidates.append(
            (
                residual_rmsd,
                images,
                moving_atoms,
                rotation_angle,
                geometry_qc,
            )
        )
    if not candidates:
        raise errors[0] if errors else ValueError(
            "no connectivity-preserving dihedral path was generated"
        )
    (
        final_residual_rmsd,
        images,
        moving_atoms,
        rotation_angle,
        geometry_qc,
    ) = min(candidates, key=lambda item: item[0])
    return InitializedReactionPath(
        images=images,
        metadata={
            "strategy": "acyclic_dihedral_rotation",
            "coordinate_kind": "dihedral",
            "coordinate_atoms": list(term.atoms),
            "central_bond_atoms": list(term.atoms[1:3]),
            "rotated_atom_indices": moving_atoms,
            "source_value_rad": start_value,
            "target_value_rad": end_value,
            "wrapped_target_delta_rad": target_delta,
            "applied_rotation_rad": rotation_angle,
            "endpoint_residual_rmsd_A": final_residual_rmsd,
            "image_count": count,
            **geometry_qc,
        },
    )

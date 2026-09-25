"""Deterministic rigid-fragment seeds for recovering molecular endpoint basins."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np

from hfauto.chemistry.basin_identity import permutation_invariant_distance_rmsd
from hfauto.chemistry.connectivity import covalent_fragments
from hfauto.chemistry.xyz import XYZ


def _unit(vector: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if norm <= 1.0e-12:
        raise ValueError("cannot normalize a near-zero vector")
    return np.asarray(vector, dtype=float) / norm


def _rotation(axis: np.ndarray, angle_degrees: float) -> np.ndarray:
    direction = _unit(axis)
    angle = math.radians(float(angle_degrees))
    cross = np.array(
        [
            [0.0, -direction[2], direction[1]],
            [direction[2], 0.0, -direction[0]],
            [-direction[1], direction[0], 0.0],
        ]
    )
    return (
        np.eye(3) * math.cos(angle)
        + (1.0 - math.cos(angle)) * np.outer(direction, direction)
        + math.sin(angle) * cross
    )


def _perpendicular_axis(direction: np.ndarray, azimuth_degrees: float) -> np.ndarray:
    radial = _unit(direction)
    axes = np.eye(3)
    seed = min(axes, key=lambda value: abs(float(np.dot(value, radial))))
    first = _unit(seed - float(np.dot(seed, radial)) * radial)
    second = _unit(np.cross(radial, first))
    azimuth = math.radians(float(azimuth_degrees))
    return _unit(math.cos(azimuth) * first + math.sin(azimuth) * second)


def _minimum_interfragment_distance(
    coordinates: np.ndarray, fragments: list[list[int]]
) -> float | None:
    distances: list[float] = []
    for first_index, first in enumerate(fragments):
        for second in fragments[first_index + 1 :]:
            delta = (
                coordinates[np.asarray(first), None, :]
                - coordinates[None, np.asarray(second), :]
            )
            distances.append(float(np.min(np.linalg.norm(delta, axis=2))))
    return min(distances) if distances else None


def reaction_atom_indices(reaction_data: dict[str, Any]) -> list[int]:
    """Return stable atom indices mentioned by bond or coordinate terms."""

    indices: list[int] = []
    for change in reaction_data.get("bond_changes", []) or []:
        indices.extend(int(value) for value in change.get("atoms", []))
    coordinate = reaction_data.get("reaction_coordinate") or {}
    for term in coordinate.get("terms", []) or []:
        indices.extend(int(value) for value in term.get("atoms", []))
    atoms = coordinate.get("atoms") or {}
    if isinstance(atoms, dict):
        indices.extend(int(value) for value in atoms.values())
    return list(dict.fromkeys(indices))


def generate_rigid_fragment_seeds(
    source: XYZ,
    reaction_data: dict[str, Any],
    *,
    angle_degrees: tuple[float, ...] = (25.0, -25.0, 50.0, -50.0),
    radial_shifts_A: tuple[float, ...] = (0.0,),
    minimum_interfragment_distance_A: float = 0.70,
    maximum_seeds: int = 4,
) -> list[tuple[XYZ, dict[str, Any]]]:
    """Rotate disconnected fragments around a reactive core without distortion.

    Connectivity is used only to identify rigid fragments.  The operation does
    not assign chemical states, bonds, minima, or mechanisms.
    """

    if maximum_seeds < 1:
        raise ValueError("maximum_seeds must be positive")
    fragments = covalent_fragments(source)
    if len(fragments) < 2:
        return []
    mentioned = set(reaction_atom_indices(reaction_data))
    core = max(
        fragments,
        key=lambda fragment: (
            len(mentioned.intersection(fragment)),
            len(fragment),
            -min(fragment),
        ),
    )
    # Prefer the largest reactive fragment when several fragments contain a
    # declared atom. This keeps the molecular scaffold fixed for complexes.
    largest = fragments[0]
    if mentioned.intersection(largest):
        core = largest
    mobile = [fragment for fragment in fragments if fragment != core]
    core_reactive = sorted(mentioned.intersection(core))
    anchor = (
        np.asarray(source.coords[core_reactive[0]], dtype=float)
        if core_reactive
        else np.asarray(source.coords[core], dtype=float).mean(axis=0)
    )
    templates = [
        (float(angle), float(radial_shift))
        for radial_shift in radial_shifts_A
        for angle in angle_degrees
    ]
    seeds: list[tuple[XYZ, dict[str, Any]]] = []
    for template_index, (angle, radial_shift) in enumerate(templates):
        coordinates = np.asarray(source.coords, dtype=float).copy()
        transforms: list[dict[str, Any]] = []
        for mobile_index, fragment in enumerate(mobile):
            fragment_indices = np.asarray(fragment, dtype=int)
            centroid = coordinates[fragment_indices].mean(axis=0)
            radial = centroid - anchor
            if float(np.linalg.norm(radial)) <= 1.0e-10:
                continue
            azimuth = (template_index * 137.507764 + mobile_index * 83.0) % 360.0
            axis = _perpendicular_axis(radial, azimuth)
            signed_angle = angle * (1.0 if mobile_index % 2 == 0 else -1.0)
            rotation = _rotation(axis, signed_angle)
            rotated = (coordinates[fragment_indices] - anchor) @ rotation.T
            rotated_radial = centroid - anchor
            rotated_center = rotated.mean(axis=0)
            if radial_shift:
                rotated += _unit(rotated_center) * radial_shift
            coordinates[fragment_indices] = rotated + anchor
            transforms.append(
                {
                    "fragment_atom_indices": fragment,
                    "angle_degrees": signed_angle,
                    "azimuth_degrees": azimuth,
                    "radial_shift_A": radial_shift,
                    "original_center_distance_A": float(
                        np.linalg.norm(rotated_radial)
                    ),
                }
            )
        minimum_distance = _minimum_interfragment_distance(
            coordinates, fragments
        )
        if (
            minimum_distance is not None
            and minimum_distance < float(minimum_interfragment_distance_A)
        ):
            continue
        seed = XYZ(
            symbols=list(source.symbols),
            coords=coordinates,
            comment=(
                "state=endpoint_minimum_candidate "
                f"fragment_orientation_seed={template_index}"
            ),
        )
        seeds.append(
            (
                seed,
                {
                    "template_index": template_index,
                    "fragment_count": len(fragments),
                    "core_atom_indices": core,
                    "mobile_fragment_count": len(mobile),
                    "minimum_interfragment_distance_A": minimum_distance,
                    "transforms": transforms,
                    "chemical_state_assigned": False,
                    "minimum_claimed": False,
                },
            )
        )
        if len(seeds) >= int(maximum_seeds):
            break
    return seeds


def select_diverse_seed_records(
    records: list[dict[str, Any]],
    *,
    maximum_selected: int = 2,
    energy_window_kcal_mol: float = 3.0,
    minimum_distance_rmsd_A: float = 0.03,
) -> list[dict[str, Any]]:
    """Select low-energy, geometrically distinct preoptimized seeds."""

    if maximum_selected < 1:
        raise ValueError("maximum_selected must be positive")
    eligible = [
        record
        for record in records
        if record.get("accepted") is True
        and record.get("energy_hartree") is not None
        and record.get("xyz_path")
        and Path(str(record["xyz_path"])).is_file()
    ]
    if not eligible:
        return []
    eligible.sort(
        key=lambda record: (
            float(record["energy_hartree"]), str(record["species_id"])
        )
    )
    minimum_energy = float(eligible[0]["energy_hartree"])
    for record in eligible:
        record["relative_energy_kcal_mol"] = (
            float(record["energy_hartree"]) - minimum_energy
        ) * 627.5094740631
    window = [
        record
        for record in eligible
        if float(record["relative_energy_kcal_mol"])
        <= float(energy_window_kcal_mol)
    ]
    selected = [window[0]]
    while len(selected) < int(maximum_selected):
        candidates: list[tuple[float, float, str, dict[str, Any]]] = []
        for record in window:
            if record in selected:
                continue
            distances = [
                permutation_invariant_distance_rmsd(
                    record["xyz_path"], existing["xyz_path"]
                )
                for existing in selected
            ]
            finite = [value for value in distances if value is not None]
            diversity = min(finite) if finite else 0.0
            if diversity < float(minimum_distance_rmsd_A):
                continue
            candidates.append(
                (
                    -diversity,
                    float(record["relative_energy_kcal_mol"]),
                    str(record["species_id"]),
                    record,
                )
            )
        if not candidates:
            break
        candidates.sort(key=lambda item: item[:3])
        selected.append(candidates[0][3])
    return selected

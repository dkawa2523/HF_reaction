"""Rigid-fragment assembly for chemistry-neutral encounter complexes."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from hfauto.chemistry.connectivity import covalent_adjacency
from hfauto.chemistry.xyz import XYZ, write_xyz


def _unit(vector: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if norm <= 1.0e-12:
        raise ValueError("placement direction must be non-zero")
    return np.asarray(vector, dtype=float) / norm


def _rotation_between(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    source_u = _unit(source)
    target_u = _unit(target)
    cross = np.cross(source_u, target_u)
    cosine = float(np.clip(np.dot(source_u, target_u), -1.0, 1.0))
    sine = float(np.linalg.norm(cross))
    if sine <= 1.0e-12:
        if cosine > 0.0:
            return np.eye(3)
        axis_seed = np.array([1.0, 0.0, 0.0])
        if abs(float(np.dot(source_u, axis_seed))) > 0.9:
            axis_seed = np.array([0.0, 1.0, 0.0])
        axis = _unit(np.cross(source_u, axis_seed))
        return 2.0 * np.outer(axis, axis) - np.eye(3)
    axis = cross / sine
    cross_matrix = np.array(
        [[0.0, -axis[2], axis[1]], [axis[2], 0.0, -axis[0]], [-axis[1], axis[0], 0.0]]
    )
    return np.eye(3) + sine * cross_matrix + (1.0 - cosine) * (cross_matrix @ cross_matrix)


def _axis_rotation(axis: np.ndarray, degrees: float) -> np.ndarray:
    direction = _unit(axis)
    angle = np.deg2rad(float(degrees))
    cross_matrix = np.array(
        [
            [0.0, -direction[2], direction[1]],
            [direction[2], 0.0, -direction[0]],
            [-direction[1], direction[0], 0.0],
        ]
    )
    return (
        np.eye(3) * np.cos(angle)
        + (1.0 - np.cos(angle)) * np.outer(direction, direction)
        + np.sin(angle) * cross_matrix
    )


def default_placement_direction(host: XYZ, host_anchor_index: int) -> np.ndarray:
    """Estimate an open coordination direction from the local geometry."""

    if not 0 <= int(host_anchor_index) < len(host.symbols):
        raise IndexError("host anchor index is outside the current complex")
    anchor = np.asarray(host.coords[int(host_anchor_index)], dtype=float)
    vectors = np.asarray(host.coords, dtype=float) - anchor
    distances = np.linalg.norm(vectors, axis=1)
    bonded = np.asarray(covalent_adjacency(host)[host_anchor_index], dtype=bool)
    if np.any(bonded):
        occupied = np.sum(vectors[bonded] / distances[bonded, None], axis=0)
        if float(np.linalg.norm(occupied)) > 0.20:
            return _unit(-occupied)
    center = np.asarray(host.coords, dtype=float).mean(axis=0)
    direction = anchor - center
    if float(np.linalg.norm(direction)) <= 1.0e-12:
        direction = np.array([1.0, 0.0, 0.0])
    return _unit(direction)


def _automatic_placement_directions(
    host: XYZ,
    host_anchor_index: int,
) -> list[np.ndarray]:
    """Return deterministic collision alternatives after the chemical guess."""

    preferred = default_placement_direction(host, host_anchor_index)
    anchor = np.asarray(host.coords[int(host_anchor_index)], dtype=float)
    center = np.asarray(host.coords, dtype=float).mean(axis=0)
    candidates = [preferred]
    if float(np.linalg.norm(anchor - center)) > 1.0e-12:
        candidates.append(_unit(anchor - center))
    candidates.extend(
        _unit(np.asarray(direction, dtype=float))
        for direction in (
            (1, 0, 0),
            (-1, 0, 0),
            (0, 1, 0),
            (0, -1, 0),
            (0, 0, 1),
            (0, 0, -1),
            (1, 1, 1),
            (1, 1, -1),
            (1, -1, 1),
            (1, -1, -1),
            (-1, 1, 1),
            (-1, 1, -1),
            (-1, -1, 1),
            (-1, -1, -1),
        )
    )
    unique: list[np.ndarray] = []
    for candidate in candidates:
        if not any(float(np.dot(candidate, item)) > 1.0 - 1.0e-8 for item in unique):
            unique.append(candidate)
    return unique


def append_fragment(
    host: XYZ,
    guest: XYZ,
    *,
    host_anchor_index: int,
    guest_anchor_index: int,
    distance_A: float,
    direction: list[float] | tuple[float, float, float] | np.ndarray | None = None,
    guest_orientation_atom_index: int | None = None,
    orientation: str = "toward_host",
    twist_degrees: float = 0.0,
    minimum_interfragment_distance_A: float | None = None,
) -> tuple[XYZ, dict[str, Any]]:
    """Append one rigid guest fragment at a deterministic anchored placement."""

    host_anchor_index = int(host_anchor_index)
    guest_anchor_index = int(guest_anchor_index)
    if not 0 <= guest_anchor_index < len(guest.symbols):
        raise IndexError("guest anchor index is outside the guest fragment")
    if float(distance_A) <= 0.0:
        raise ValueError("distance_A must be positive")
    centered = np.asarray(guest.coords, dtype=float) - guest.coords[guest_anchor_index]
    if guest_orientation_atom_index is not None:
        orientation_index = int(guest_orientation_atom_index)
        if not 0 <= orientation_index < len(guest.symbols) or orientation_index == guest_anchor_index:
            raise IndexError("guest orientation atom must differ from the guest anchor")
    else:
        orientation_index = None
    if orientation not in {"toward_host", "away_from_host"}:
        raise ValueError("orientation must be 'toward_host' or 'away_from_host'")

    directions = (
        _automatic_placement_directions(host, host_anchor_index)
        if direction is None
        else [_unit(np.asarray(direction, dtype=float))]
    )
    placements: list[tuple[float, np.ndarray, np.ndarray]] = []
    for placement_direction in directions:
        rotation = np.eye(3)
        if orientation_index is not None:
            target = (
                -placement_direction
                if orientation == "toward_host"
                else placement_direction
            )
            rotation = _rotation_between(centered[orientation_index], target)
        rotated = centered @ rotation.T
        if abs(float(twist_degrees)) > 1.0e-12:
            rotated = rotated @ _axis_rotation(
                placement_direction, twist_degrees
            ).T
        target_anchor = (
            host.coords[host_anchor_index]
            + placement_direction * float(distance_A)
        )
        placed = rotated + target_anchor
        interfragment = np.linalg.norm(
            np.asarray(host.coords)[:, None, :] - placed[None, :, :], axis=2
        )
        placements.append(
            (float(np.min(interfragment)), placement_direction, placed)
        )

    preferred_minimum, preferred_direction, preferred_placed = placements[0]
    required_clearance = (
        float(minimum_interfragment_distance_A)
        if minimum_interfragment_distance_A is not None
        else None
    )
    if required_clearance is None or preferred_minimum >= required_clearance:
        minimum_distance, placement_direction, placed = (
            preferred_minimum,
            preferred_direction,
            preferred_placed,
        )
        preferred_used = True
    else:
        minimum_distance, placement_direction, placed = max(
            placements, key=lambda item: item[0]
        )
        preferred_used = False
    combined = XYZ(
        symbols=[*host.symbols, *guest.symbols],
        coords=np.vstack([host.coords, placed]),
        comment=host.comment,
    )
    return combined, {
        "host_anchor_index": host_anchor_index,
        "guest_anchor_index": guest_anchor_index,
        "distance_A": float(distance_A),
        "direction": [float(value) for value in placement_direction],
        "minimum_interfragment_distance_A": minimum_distance,
        "preferred_direction_used": preferred_used,
        "direction_candidates_evaluated": len(placements),
    }


def element_counts(symbols: list[str]) -> dict[str, int]:
    return dict(sorted(Counter(symbols).items()))


def write_complex(xyz: XYZ, path: str | Path, *, system_id: str) -> Path:
    xyz.comment = f"state=assembled_complex system_id={system_id} generated_by=hfauto"
    return write_xyz(xyz, path)

"""Geometry-derived covalent connectivity used across chemistry workflows."""

from __future__ import annotations

import numpy as np

from hfauto.chemistry.xyz import XYZ

_COVALENT_RADII_A = {
    "H": 0.31,
    "B": 0.84,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "Si": 1.11,
    "P": 1.07,
    "S": 1.05,
    "Cl": 1.02,
    "Br": 1.20,
    "I": 1.39,
}


def covalent_radius_A(symbol: str) -> float:
    """Return the radius used by the geometry-derived connectivity model."""

    return float(_COVALENT_RADII_A.get(str(symbol), 0.77))


def covalent_adjacency(xyz: XYZ) -> np.ndarray:
    """Return a conservative geometry-derived covalent graph."""

    size = len(xyz.symbols)
    adjacency = np.zeros((size, size), dtype=bool)
    for index in range(size):
        for other_index in range(index):
            first_radius = covalent_radius_A(xyz.symbols[index])
            second_radius = covalent_radius_A(xyz.symbols[other_index])
            distance = float(
                np.linalg.norm(xyz.coords[index] - xyz.coords[other_index])
            )
            bonded = distance <= 1.25 * (first_radius + second_radius)
            adjacency[index, other_index] = bonded
            adjacency[other_index, index] = bonded
    return adjacency


def covalent_fragments(xyz: XYZ) -> list[list[int]]:
    """Return deterministic connected components of the covalent graph."""

    adjacency = covalent_adjacency(xyz)
    unseen = set(range(len(xyz.symbols)))
    fragments: list[list[int]] = []
    while unseen:
        start = min(unseen)
        stack = [start]
        component: list[int] = []
        unseen.remove(start)
        while stack:
            index = stack.pop()
            component.append(index)
            neighbours = [
                other for other in sorted(unseen) if adjacency[index, other]
            ]
            for other in neighbours:
                unseen.remove(other)
                stack.append(other)
        fragments.append(sorted(component))
    return sorted(fragments, key=lambda value: (-len(value), value))


def connectivity_changes(initial: XYZ, final: XYZ) -> dict[str, object]:
    """Describe atom-mapped covalent-graph changes between two geometries.

    This is a structural observation, not a bond-order calculation.  It is
    intended to keep an unconstrained relaxation that changes connectivity
    from being silently relabelled as the original chemical state.
    """

    if initial.symbols != final.symbols:
        return {
            "accepted": False,
            "composition_preserved": sorted(initial.symbols) == sorted(final.symbols),
            "atom_order_preserved": False,
            "topology_changed": None,
            "formed_bonds": [],
            "broken_bonds": [],
            "reason": "atom_order_or_composition_changed",
        }
    initial_graph = covalent_adjacency(initial)
    final_graph = covalent_adjacency(final)
    formed: list[list[int]] = []
    broken: list[list[int]] = []
    for first in range(len(initial.symbols)):
        for second in range(first):
            pair = [second, first]
            if not initial_graph[first, second] and final_graph[first, second]:
                formed.append(pair)
            elif initial_graph[first, second] and not final_graph[first, second]:
                broken.append(pair)
    return {
        "accepted": True,
        "composition_preserved": True,
        "atom_order_preserved": True,
        "topology_changed": bool(formed or broken),
        "formed_bonds": formed,
        "broken_bonds": broken,
        "initial_fragment_count": len(covalent_fragments(initial)),
        "final_fragment_count": len(covalent_fragments(final)),
        "reason": "connectivity_changed" if formed or broken else "same_connectivity",
    }


def fragments_after_bond_cut(
    xyz: XYZ, first: int, second: int
) -> tuple[list[int], list[int]] | None:
    """Return the two fragments made by cutting one acyclic covalent bond.

    ``None`` means the atom pair is not a detected bond or remains connected
    through a ring.  Path builders can therefore fail closed instead of
    rotating an arbitrary subset of a cyclic molecule.
    """

    size = len(xyz.symbols)
    first_index, second_index = int(first), int(second)
    if (
        first_index < 0
        or second_index < 0
        or first_index >= size
        or second_index >= size
        or first_index == second_index
    ):
        raise ValueError("bond-cut atom indices are invalid")
    adjacency = covalent_adjacency(xyz)
    if not adjacency[first_index, second_index]:
        return None
    adjacency[first_index, second_index] = False
    adjacency[second_index, first_index] = False

    def component(start: int) -> set[int]:
        visited = {start}
        stack = [start]
        while stack:
            atom = stack.pop()
            for neighbour in np.flatnonzero(adjacency[atom]):
                value = int(neighbour)
                if value not in visited:
                    visited.add(value)
                    stack.append(value)
        return visited

    first_fragment = component(first_index)
    if second_index in first_fragment:
        return None
    second_fragment = component(second_index)
    if first_fragment | second_fragment != set(range(size)):
        return None
    return sorted(first_fragment), sorted(second_fragment)

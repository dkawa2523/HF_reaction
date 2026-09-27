"""Covalent topology with hysteresis, fragments, labile H, acceptors, WL atom classes, state
labels (§5.5) and the value and gradient of a declared reaction coordinate.

A pair is bonded at r ≤ 1.15 Σr_cov and non-bonded at r ≥ 1.45 Σr_cov.  In between it keeps
its previous state; with no previous state it counts as bonded, so strong and symmetric
hydrogen bonds such as FHF⁻ stay one fragment.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from collections.abc import Collection, Sequence

import numpy as np

from hfauto.chemistry.elements import covalent_radius
from hfauto.chemistry.geometry import dihedral_deg
from hfauto.chemistry.xyz import hill_formula
from hfauto.core.records import CoordinateTerm

Bond = tuple[int, int]  # (i, j) with i < j
_FD_STEP_A = 1.0e-4  # central-difference step of a declared coordinate's gradient

BOND_MAX_RATIO = 1.15
NONBOND_MIN_RATIO = 1.45
_WL_ITERATIONS = 3

_LABILE_PARTNERS = frozenset({"N", "O", "F", "S", "Cl", "Br", "I", "P"})
# An atom with a lone pair left: at most this many covalent partners.
_ACCEPTOR_MAX_BONDS = {"N": 3, "P": 3, "O": 2, "S": 2, "F": 1, "Cl": 1, "Br": 1, "I": 1}


def _ratios(symbols: Sequence[str], coords: np.ndarray) -> np.ndarray:
    """r_ij / (r_cov,i + r_cov,j); the diagonal is +inf."""

    x = np.asarray(coords, dtype=float).reshape(-1, 3)
    if len(x) != len(symbols):
        raise ValueError("symbols and coordinates differ in atom count")
    radii = np.array([covalent_radius(s) for s in symbols])
    ratio = np.linalg.norm(x[:, None] - x[None], axis=-1) / (radii[:, None] + radii[None])
    np.fill_diagonal(ratio, np.inf)
    return ratio


def _pairs(mask: np.ndarray) -> frozenset[Bond]:
    i, j = np.nonzero(np.triu(mask, 1))
    return frozenset(zip(i.tolist(), j.tolist(), strict=True))


def bonds(
    symbols: Sequence[str], coords: np.ndarray, previous: Collection[Bond] | None = None
) -> frozenset[Bond]:
    ratio = _ratios(symbols, coords)
    loose = ratio < NONBOND_MIN_RATIO
    if previous is None:
        return _pairs(loose)
    kept = np.zeros_like(loose)
    for i, j in previous:
        kept[i, j] = kept[j, i] = True
    return _pairs((ratio <= BOND_MAX_RATIO) | (loose & kept))


def _components(n_atoms: int, bonded: Collection[Bond]) -> tuple[tuple[int, ...], ...]:
    parent = list(range(n_atoms))

    def root(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j in bonded:
        parent[root(i)] = root(j)
    groups: dict[int, list[int]] = {}
    for i in range(n_atoms):
        groups.setdefault(root(i), []).append(i)
    return tuple(sorted((tuple(g) for g in groups.values()), key=lambda g: (-len(g), g)))


def fragments(
    symbols: Sequence[str], coords: np.ndarray, previous: Collection[Bond] | None = None
) -> tuple[tuple[int, ...], ...]:
    """Connected components, largest first."""

    return _components(len(symbols), bonds(symbols, coords, previous))


def bond_changes(
    symbols: Sequence[str], a: np.ndarray, b: np.ndarray
) -> tuple[frozenset[Bond], frozenset[Bond]]:
    """(formed, broken) going from a to b, with a's bonds as the previous state of b."""

    before = bonds(symbols, a)
    after = bonds(symbols, b, previous=before)
    return after - before, before - after


def declared_coordinate(terms: Sequence[CoordinateTerm], x: np.ndarray) -> float:
    """Σ coefficient × (distance Å | angle ° | dihedral °)."""
    total = 0.0
    for term in terms:
        p = np.asarray(x, dtype=float).reshape(-1, 3)[list(term.atoms)]
        if term.kind == "distance":
            value = float(np.linalg.norm(p[1] - p[0]))
        elif term.kind == "angle":
            u, v = p[0] - p[1], p[2] - p[1]
            value = float(np.degrees(np.arccos(np.clip(
                u @ v / (np.linalg.norm(u) * np.linalg.norm(v)), -1.0, 1.0))))
        else:
            value = dihedral_deg(x, term.atoms)
        total += term.coefficient * value
    return total


def declared_coordinate_gradient(terms: Sequence[CoordinateTerm], x: np.ndarray) -> np.ndarray:
    """Central-difference gradient of ``declared_coordinate``; dihedral steps wrap at ±180°."""
    flat, grad = np.asarray(x, dtype=float).ravel(), np.zeros(np.size(x))
    for i in range(flat.size):
        step = np.zeros_like(flat)
        step[i] = _FD_STEP_A
        delta = declared_coordinate(terms, flat + step) - declared_coordinate(terms, flat - step)
        grad[i] = ((delta + 180.0) % 360.0 - 180.0) / (2 * _FD_STEP_A)
    return grad


def labile_hydrogens(symbols: Sequence[str], coords: np.ndarray) -> tuple[int, ...]:
    """H atoms bonded to N/O/F/S/Cl/Br/I/P."""

    labile = {
        h
        for pair in bonds(symbols, coords)
        for h, partner in (pair, pair[::-1])
        if symbols[h] == "H" and symbols[partner] in _LABILE_PARTNERS
    }
    return tuple(sorted(labile))


def acceptor_atoms(symbols: Sequence[str], coords: np.ndarray) -> tuple[int, ...]:
    """N/P/O/S/halogen atoms that still carry a lone pair (few enough covalent partners)."""

    degree = Counter(i for pair in bonds(symbols, coords) for i in pair)
    return tuple(
        i
        for i, s in enumerate(symbols)
        if s in _ACCEPTOR_MAX_BONDS and degree[i] <= _ACCEPTOR_MAX_BONDS[s]
    )


def _wl_rounds(symbols: Sequence[str], bonded: Collection[Bond]) -> list[list[str]]:
    """Weisfeiler–Lehman atom labels of the element-labelled graph, round 0 (elements) to 3."""

    neighbours: list[list[int]] = [[] for _ in symbols]
    for i, j in bonded:
        neighbours[i].append(j)
        neighbours[j].append(i)
    rounds = [list(symbols)]
    for _ in range(_WL_ITERATIONS):
        labels = rounds[-1]
        rounds.append([
            hashlib.sha256(
                f"{labels[i]}|{','.join(sorted(labels[j] for j in neighbours[i]))}".encode()
            ).hexdigest()[:16]
            for i in range(len(labels))
        ])
    return rounds


def wl_classes(symbols: Sequence[str], bonded: Collection[Bond]) -> tuple[int, ...]:
    """Atom equivalence class per atom (last WL round); the numbering is permutation invariant."""

    last = _wl_rounds(symbols, bonded)[-1]
    rank = {label: k for k, label in enumerate(sorted(set(last)))}
    return tuple(rank[label] for label in last)


def state_label(symbols: Sequence[str], coords: np.ndarray) -> str:
    """Sorted fragment formulas plus the first 8 hex digits of the hash of all WL labels."""

    bonded = bonds(symbols, coords)
    formulas = sorted(
        hill_formula([symbols[i] for i in group])
        for group in _components(len(symbols), bonded)
    )
    seen = Counter(label for labels in _wl_rounds(symbols, bonded) for label in labels)
    digest = hashlib.sha256(json.dumps(sorted(seen.items())).encode()).hexdigest()
    return f"{'+'.join(formulas)}_{digest[:8]}"

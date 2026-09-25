"""Covalent topology with hysteresis, fragments, labile H, acceptors, state labels (§5.5)
and the value and gradient of a declared reaction coordinate.

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

from hfauto.chemistry.xyz import hill_formula
from hfauto.core.records import CoordinateTerm

Bond = tuple[int, int]  # (i, j) with i < j
_FD_STEP_A = 1.0e-4  # central-difference step of a declared coordinate's gradient

BOND_MAX_RATIO = 1.15
NONBOND_MIN_RATIO = 1.45
PROTON_Q_MIN_A = 0.3
_WL_ITERATIONS = 3

# Cordero et al., Dalton Trans. 2008 (values shared with the legacy connectivity table).
_COVALENT_RADII_A = {
    "H": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71,
    "O": 0.66, "F": 0.57, "Ne": 0.58, "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "Ar": 1.06, "K": 2.03, "Ca": 1.76, "Sc": 1.70,
    "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39, "Fe": 1.32, "Co": 1.26, "Ni": 1.24,
    "Cu": 1.32, "Zn": 1.22, "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20,
    "Kr": 1.16, "I": 1.39,
}
# Bondi, J. Phys. Chem. 1964; Be, B, Al, Ca, Ge from Mantina et al., JPCA 2009.
_VDW_RADII_A = {
    "H": 1.20, "He": 1.40, "Li": 1.82, "Be": 1.53, "B": 1.92, "C": 1.70, "N": 1.55,
    "O": 1.52, "F": 1.47, "Ne": 1.54, "Na": 2.27, "Mg": 1.73, "Al": 1.84, "Si": 2.10,
    "P": 1.80, "S": 1.80, "Cl": 1.75, "Ar": 1.88, "K": 2.75, "Ca": 2.31, "Ni": 1.63,
    "Cu": 1.40, "Zn": 1.39, "Ga": 1.87, "Ge": 2.11, "As": 1.85, "Se": 1.90, "Br": 1.85,
    "Kr": 2.02, "I": 1.98,
}
_LABILE_PARTNERS = frozenset({"N", "O", "F", "S", "Cl", "Br", "I", "P"})
# An atom with a lone pair left: at most this many covalent partners.
_ACCEPTOR_MAX_BONDS = {"N": 3, "P": 3, "O": 2, "S": 2, "F": 1, "Cl": 1, "Br": 1, "I": 1}


def _lookup(table: dict[str, float], symbol: str, what: str) -> float:
    try:
        return table[symbol]
    except KeyError:
        raise ValueError(f"no {what} radius tabulated for element {symbol!r}") from None


def vdw_radius(symbol: str) -> float:
    return _lookup(_VDW_RADII_A, symbol, "van der Waals")


def _ratios(symbols: Sequence[str], coords: np.ndarray) -> np.ndarray:
    """r_ij / (r_cov,i + r_cov,j); the diagonal is +inf."""

    x = np.asarray(coords, dtype=float).reshape(-1, 3)
    if len(x) != len(symbols):
        raise ValueError("symbols and coordinates differ in atom count")
    radii = np.array([_lookup(_COVALENT_RADII_A, s, "covalent") for s in symbols])
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


def proton_coordinate(coords: np.ndarray, d: int, h: int, a: int) -> float:
    """q = r(D–H) − r(A–H) in Å."""

    x = np.asarray(coords, dtype=float).reshape(-1, 3)
    return float(np.linalg.norm(x[d] - x[h]) - np.linalg.norm(x[a] - x[h]))


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
            b0, b1, b2 = p[0] - p[1], p[2] - p[1], p[3] - p[2]
            b1 = b1 / np.linalg.norm(b1)
            v, w = b0 - (b0 @ b1) * b1, b2 - (b2 @ b1) * b1
            value = float(np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w)))
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


def transferred_hydrogens(symbols: Sequence[str], a: np.ndarray, b: np.ndarray) -> int:
    """Number of H whose covalent partners change.

    An H that swaps one donor D for one acceptor A counts only if q flips sign with
    |q| ≥ 0.3 Å at both ends.
    """

    formed, broken = bond_changes(symbols, a, b)
    count = 0
    for h in (i for i, s in enumerate(symbols) if s == "H"):
        gained = [i + j - h for i, j in formed if h in (i, j)]
        lost = [i + j - h for i, j in broken if h in (i, j)]
        if not gained and not lost:
            continue
        if len(gained) == 1 and len(lost) == 1:
            qa = proton_coordinate(a, lost[0], h, gained[0])
            qb = proton_coordinate(b, lost[0], h, gained[0])
            if min(abs(qa), abs(qb)) < PROTON_Q_MIN_A or qa * qb > 0:
                continue
        count += 1
    return count


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


def _wl_hash(symbols: Sequence[str], bonded: Collection[Bond]) -> str:
    """Weisfeiler–Lehman hash of the element-labelled graph (permutation invariant)."""

    neighbours: list[list[int]] = [[] for _ in symbols]
    for i, j in bonded:
        neighbours[i].append(j)
        neighbours[j].append(i)
    labels = list(symbols)
    seen = Counter(labels)
    for _ in range(_WL_ITERATIONS):
        labels = [
            hashlib.sha256(
                f"{labels[i]}|{','.join(sorted(labels[j] for j in neighbours[i]))}".encode()
            ).hexdigest()[:16]
            for i in range(len(labels))
        ]
        seen.update(labels)
    return hashlib.sha256(json.dumps(sorted(seen.items())).encode()).hexdigest()


def state_label(symbols: Sequence[str], coords: np.ndarray) -> str:
    """Sorted fragment formulas plus the first 8 hex digits of the WL hash."""

    bonded = bonds(symbols, coords)
    formulas = sorted(
        hill_formula([symbols[i] for i in group])
        for group in _components(len(symbols), bonded)
    )
    return f"{'+'.join(formulas)}_{_wl_hash(symbols, bonded)[:8]}"

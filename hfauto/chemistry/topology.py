"""Covalent bond graph, fragments, bond changes, WL atom classes and state labels (§5.5).

A pair is bonded when r < r_thr = r_cov,i + r_cov,j + 0.4 Å (additive tolerance, Meng & Lewis
1991; SCINE BondDetector, OpenBabel). The rule sees one structure only, so ``bond_changes`` is
symmetric: FHF⁻ (1.14 Å) and I3⁻ (2.92 Å) are bonded, a halogen bond I···N at 2.8 Å is not.
``resolved_bond_changes`` keeps only the changes that clear a band of ±RESOLVED_A around r_thr
on both sides, so a contact the threshold cuts within one basin is no change.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from collections.abc import Collection, Sequence

import numpy as np

from hfauto.chemistry.elements import covalent_radius
from hfauto.chemistry.geometry import neighbours
from hfauto.chemistry.xyz import hill_formula

Bond = tuple[int, int]  # (i, j) with i < j
BOND_TOLERANCE_A = 0.4
# Half-width of the band around r_thr that resolves a bond change: the GFN2 and PBE0 N···H of
# the S19 amine·HF basin differ by 0.15 Å, so a crossing within 0.1 Å of r_thr is method noise.
RESOLVED_A = 0.1
_WL_ITERATIONS = 3


def _excess(symbols: Sequence[str], coords: np.ndarray) -> np.ndarray:
    """r − r_thr (Å) per atom pair: negative where bonded."""

    x = np.asarray(coords, dtype=float).reshape(-1, 3)
    if len(x) != len(symbols):
        raise ValueError("symbols and coordinates differ in atom count")
    radii = np.array([covalent_radius(s) for s in symbols])
    r = np.linalg.norm(x[:, None] - x[None], axis=-1)
    return r - (np.add.outer(radii, radii) + BOND_TOLERANCE_A)


def _pairs(mask: np.ndarray) -> frozenset[Bond]:
    i, j = np.nonzero(np.triu(mask, 1))
    return frozenset(zip(i.tolist(), j.tolist(), strict=True))


def bonds(symbols: Sequence[str], coords: np.ndarray) -> frozenset[Bond]:
    return _pairs(_excess(symbols, coords) < 0)


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


def fragments(symbols: Sequence[str], coords: np.ndarray) -> tuple[tuple[int, ...], ...]:
    """Connected components, largest first."""

    return _components(len(symbols), bonds(symbols, coords))


def bond_changes(
    symbols: Sequence[str], a: np.ndarray, b: np.ndarray
) -> tuple[frozenset[Bond], frozenset[Bond]]:
    """(formed, broken) going from a to b: the two set differences of the bond graphs."""

    before, after = bonds(symbols, a), bonds(symbols, b)
    return after - before, before - after


def resolved_bond_changes(
    symbols: Sequence[str], x: np.ndarray, y: np.ndarray
) -> tuple[frozenset[Bond], frozenset[Bond]]:
    """(formed, broken) going from x to y, only the pairs whose r − r_thr is ≥ +RESOLVED_A in
    one structure and ≤ −RESOLVED_A in the other (both in one atom order)."""

    ex, ey = _excess(symbols, x), _excess(symbols, y)
    return (_pairs((ex >= RESOLVED_A) & (ey <= -RESOLVED_A)),
            _pairs((ex <= -RESOLVED_A) & (ey >= RESOLVED_A)))


def _wl_rounds(symbols: Sequence[str], bonded: Collection[Bond]) -> list[list[str]]:
    """Weisfeiler–Lehman atom labels of the element-labelled graph, round 0 (elements) to 3."""

    partners = neighbours(bonded)
    rounds = [list(symbols)]
    for _ in range(_WL_ITERATIONS):
        labels = rounds[-1]
        rounds.append([
            hashlib.sha256(
                f"{labels[i]}|{','.join(sorted(labels[j] for j in partners.get(i, ())))}".encode()
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

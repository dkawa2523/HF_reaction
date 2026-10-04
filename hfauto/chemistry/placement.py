"""Placement seeds of a complex and the small-rigid test of the conformers stage (design §8.2).

Each guest is placed as a rigid body in a random orientation and direction (one fixed RNG seed)
at van der Waals contact + 0.5 Å from the atoms placed so far, keeping ``n_seeds`` pairwise
distinct seeds; three or more fragments are stacked one at a time. Directions are taken in
molecule-fixed frames, so the seeds do not depend on the input orientation; ties between
symmetric atoms in a frame follow atom order, so for symmetric fragments they depend on atom
numbering. CREST ``--nci`` samples from seed00, and every seed is the composition's output when
no CREST candidate keeps its state label. Observed in VAL9 R2 / B1: CREST 3.0.2 failed on
CH3·O2 (rc -11) and OH·CH4 (rc 1), and the GFN2 screen optimized only the seeds that collapse
into another bonding state (CH3O2 from 1 of 6 seeds; CH3 + H2O from 2 of 6; the others fail
SCC or maxiter). Those collapses carry both systems' discovery: one contact seed fails the screen
for both.
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence

import numpy as np

from hfauto.chemistry.elements import vdw_radius
from hfauto.chemistry.identity import permutation_invariant_rmsd
from hfauto.chemistry.topology import bonds
from hfauto.chemistry.xyz import XYZ

CONTACT_PAD_A = 0.5  # closest contact = Σ r_vdw + 0.5 Å
MAX_TRIES = 50  # placements tried per cluster
DUPLICATE_RMSD_A = 0.1
RNG_SEED = 20260925
_TINY = 1e-3


def _unit(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v)


def _perpendicular(axis: np.ndarray, hint: np.ndarray | None) -> np.ndarray:
    """Unit vector ⟂ axis along the perpendicular part of ``hint``; a canonical one otherwise
    (only for fragments symmetric about the axis, where every choice is equivalent)."""
    if hint is not None:
        p = hint - (hint @ axis) * axis
        if np.linalg.norm(p) > _TINY:
            return _unit(p)
    trial = np.eye(3)[int(np.argmin(np.abs(axis)))]
    return _unit(trial - (trial @ axis) * axis)


def _first_max(values: np.ndarray) -> int:
    """First index within 1e-3 of the maximum: ties between symmetric atoms follow atom order,
    not rounding noise of the input orientation."""
    return int(np.argmax(values >= values.max() - _TINY))


def _body_frame(x: np.ndarray) -> np.ndarray:
    """Orthonormal rows (axis, e1, e2) at the centroid: the axis points to the (first) farthest
    atom, e1 to the atoms' side of it (or to the atom farthest from it when they balance)."""
    c = x - x.mean(axis=0)
    r = np.linalg.norm(c, axis=1)
    if r.max() < 0.1:
        return np.eye(3)
    axis = _unit(c[_first_max(r)])
    off = c - np.outer(c @ axis, axis)
    hint = off.sum(axis=0)
    if np.linalg.norm(hint) < 0.1:
        dist = np.linalg.norm(off, axis=1)
        hint = off[_first_max(dist)] if dist.max() > 0.1 else None
    e1 = _perpendicular(axis, hint)
    return np.array([axis, e1, np.cross(axis, e1)])


def _random_rotation(rng: np.random.Generator) -> np.ndarray:
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q = q * np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return q


def _placements(cluster: XYZ, guest: XYZ, rng: np.random.Generator) -> Iterator[XYZ]:
    """Random orientation and direction in the cluster frame; contact at Σ r_vdw + 0.5 Å (at
    least 2.9 Å for any pair of the element table, so no further clash test)."""
    frame = _body_frame(cluster.coords)
    center = cluster.coords.mean(axis=0)
    body = (guest.coords - guest.coords.mean(axis=0)) @ _body_frame(guest.coords).T
    radii = np.add.outer([vdw_radius(s) for s in cluster.symbols],
                         [vdw_radius(s) for s in guest.symbols]) + CONTACT_PAD_A
    for _ in range(MAX_TRIES):
        oriented = body @ _random_rotation(rng).T @ frame
        direction = _unit(rng.normal(size=3)) @ frame

        def gap(s: float, oriented: np.ndarray = oriented, direction: np.ndarray = direction
                ) -> float:
            y = oriented + center + s * direction
            return float(np.min(np.linalg.norm(cluster.coords[:, None] - y[None], axis=-1) - radii))

        lo, hi = 0.0, 1.0 + float(np.ptp(cluster.coords)) + float(np.ptp(guest.coords))
        while gap(hi) < 0:
            hi *= 2.0
        for _ in range(40):  # bisection keeps gap(hi) ≥ 0
            mid = 0.5 * (lo + hi)
            lo, hi = (lo, mid) if gap(mid) >= 0 else (mid, hi)
        moved = oriented + center + hi * direction
        yield XYZ([*cluster.symbols, *guest.symbols], np.vstack([cluster.coords, moved]))


def _novel(seed: XYZ, kept: Sequence[XYZ]) -> bool:
    return all(permutation_invariant_rmsd(seed.symbols, seed.coords, k.coords)[0]
               >= DUPLICATE_RMSD_A for k in kept)


def _grow(clusters: Sequence[XYZ], guest: XYZ, rng: np.random.Generator, n: int) -> list[XYZ]:
    """Attach ``guest`` to the clusters round-robin until ``n`` distinct seeds exist."""
    streams = [_placements(c, guest, rng) for c in clusters]
    out: list[XYZ] = []
    while streams and len(out) < n:
        for stream in list(streams):
            seed = next(stream, None)
            if seed is None:
                streams.remove(stream)
            elif _novel(seed, out):
                out.append(seed)
                if len(out) == n:
                    break
    return out


def seeds(host: XYZ, guests: Sequence[XYZ], *, n_seeds: int = 6) -> list[XYZ]:
    """Up to ``n_seeds`` pairwise distinct placements of host + guests (atoms in that order)."""
    rng = np.random.default_rng(RNG_SEED)
    clusters = [XYZ(list(host.symbols), np.asarray(host.coords, dtype=float))]
    for guest in guests:
        part = XYZ(list(guest.symbols), np.asarray(guest.coords, dtype=float))
        clusters = _grow(clusters, part, rng, n_seeds)
    return clusters


def _in_ring(bonded: frozenset[tuple[int, int]], i: int, j: int) -> bool:
    """True when i and j stay connected without the bond (i, j)."""
    seen, stack = {i}, [i]
    while stack:
        k = stack.pop()
        for p, q in bonded:
            if k in (p, q) and (p, q) != (i, j) and (m := p + q - k) not in seen:
                if m == j:
                    return True
                seen.add(m)
                stack.append(m)
    return False


def rotatable_bonds(symbols: Sequence[str], coords: np.ndarray) -> int:
    """Acyclic bonds between two heavy atoms that both have another covalent partner."""
    bonded = bonds(symbols, coords)
    degree = np.bincount([k for pair in bonded for k in pair], minlength=len(symbols))
    return sum(symbols[i] != "H" and symbols[j] != "H" and degree[i] > 1 and degree[j] > 1
               and not _in_ring(bonded, i, j) for i, j in bonded)


def is_small_rigid(symbols: Sequence[str], coords: np.ndarray) -> bool:
    """At most 3 heavy atoms and no rotatable bond: the conformer search is skipped."""
    return sum(s != "H" for s in symbols) <= 3 and rotatable_bonds(symbols, coords) == 0

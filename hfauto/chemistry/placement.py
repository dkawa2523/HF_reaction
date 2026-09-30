"""Placement seeds for complexes (design §8.2 conformers, chem 10, CH-19).

With a polar H (bonded to N/O/F/S/P or a halogen) on one side and an acceptor (N/P/O/S or a
halogen with a lone pair left) on the other, the H sits on the acceptor's lone-pair cone,
110–120° from the acceptor's bonds (the lone-pair axis of an sp3 acceptor), at azimuths
0/120/240°, and the donor fragment is tilted by ±30° about an RNG-drawn axis through that H.
Otherwise the guest is placed as a rigid body in a random orientation at van der Waals
contact + 0.5 Å. All directions are taken in molecule-fixed frames, so seeds do not depend
on the input orientation. Three or more fragments are stacked one at a time. CREST ``--nci``
samples from seed00; the seeds themselves are the output when CREST fails (no GFN2 minimum
complex or an SCC failure), and there the H-bond cone keeps complex states that rigid contact seeds lose in
the xTB screen.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterator, Sequence
from itertools import product

import numpy as np

from hfauto.chemistry.elements import vdw_radius
from hfauto.chemistry.geometry import rotation_about
from hfauto.chemistry.identity import permutation_invariant_rmsd
from hfauto.chemistry.topology import bonds
from hfauto.chemistry.xyz import XYZ

CONE_DEG = (110.0, 120.0)  # A···H direction measured from the acceptor's bonds
AZIMUTHS_DEG = (0.0, 120.0, 240.0)
TILT_DEG = 30.0
HBOND_SCALE = 0.75  # H···A = 0.75 × (r_vdw(H) + r_vdw(A))
CONTACT_PAD_A = 0.5  # rigid fallback: closest contact = Σ r_vdw + 0.5 Å
MIN_H_A = 1.2  # intermolecular pairs involving H
MIN_HEAVY_A = 2.2  # intermolecular heavy-atom pairs
MAX_TRIES = 50  # rejection sampling per placement
DUPLICATE_RMSD_A = 0.1
_TINY = 1e-3
_LABILE_PARTNERS = frozenset({"N", "O", "F", "S", "Cl", "Br", "I", "P"})
# An atom with a lone pair left: at most this many covalent partners.
_ACCEPTOR_MAX_BONDS = {"N": 3, "P": 3, "O": 2, "S": 2, "F": 1, "Cl": 1, "Br": 1, "I": 1}


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


def _frame(x: np.ndarray, origin: np.ndarray, axis: np.ndarray) -> np.ndarray:
    """Molecule-fixed orthonormal rows (axis, e1, e2); e1 points to the atoms' side."""
    rel = x - origin
    off = rel - np.outer(rel @ axis, axis)
    hint = off.sum(axis=0)
    if np.linalg.norm(hint) < 0.1:  # balanced about the axis: the atom farthest from it
        dist = np.linalg.norm(off, axis=1)
        hint = off[_first_max(dist)] if dist.max() > 0.1 else None
    e1 = _perpendicular(axis, hint)
    return np.array([axis, e1, np.cross(axis, e1)])


def _body_frame(x: np.ndarray) -> np.ndarray:
    """Frame at the centroid whose axis points to the (first) farthest atom."""
    c = x - x.mean(axis=0)
    r = np.linalg.norm(c, axis=1)
    return np.eye(3) if r.max() < 0.1 else _frame(c, np.zeros(3), _unit(c[_first_max(r)]))


def _random_rotation(rng: np.random.Generator) -> np.ndarray:
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q = q * np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return q


def _partners(symbols: Sequence[str], x: np.ndarray, atom: int) -> list[int]:
    return sorted(i + j - atom for i, j in bonds(symbols, x) if atom in (i, j))


def _labile_hydrogens(m: XYZ) -> list[int]:
    return sorted({h for pair in bonds(m.symbols, m.coords) for h, partner in (pair, pair[::-1])
                   if m.symbols[h] == "H" and m.symbols[partner] in _LABILE_PARTNERS})


def _acceptor_atoms(m: XYZ) -> list[int]:
    degree = Counter(i for pair in bonds(m.symbols, m.coords) for i in pair)
    return [i for i, s in enumerate(m.symbols)
            if s in _ACCEPTOR_MAX_BONDS and degree[i] <= _ACCEPTOR_MAX_BONDS[s]]


def _lone_pair(frame: np.ndarray, bond_dirs: np.ndarray, cone_deg: float, azimuth_deg: float
               ) -> np.ndarray:
    """Direction at ``azimuth_deg`` about frame[0] (the bond-sum axis) with the smallest polar
    angle that is ≥ ``cone_deg`` from every bond of the acceptor. With one bond this is the
    cone itself; where no direction gets that far (three bonds, sp3) the direction farthest
    from all bonds is used, i.e. the lone-pair axis."""
    phi = np.radians(azimuth_deg)
    side = np.cos(phi) * frame[1] + np.sin(phi) * frame[2]
    if not len(bond_dirs):
        theta = np.radians(cone_deg)
        return np.cos(theta) * frame[0] + np.sin(theta) * side
    thetas = np.radians(np.arange(0.0, 180.01, 0.25))
    dirs = np.outer(np.cos(thetas), frame[0]) + np.outer(np.sin(thetas), side)
    angles = np.degrees(np.arccos(np.clip(dirs @ bond_dirs.T, -1.0, 1.0))).min(axis=1)
    wide = np.flatnonzero(angles >= cone_deg)
    return dirs[wide[0] if wide.size else int(np.argmax(angles))]


def _clear(sym_a: Sequence[str], xa: np.ndarray, sym_b: Sequence[str], xb: np.ndarray) -> bool:
    d = np.linalg.norm(xa[:, None] - xb[None], axis=-1)
    heavy = np.outer([s != "H" for s in sym_a], [s != "H" for s in sym_b])
    return bool(np.all(d >= np.where(heavy, MIN_HEAVY_A, MIN_H_A)))


def _hbond(donor: XYZ, h: int, acceptor: XYZ, a: int, azimuth: float, sign: float,
           rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """(R, t) moving the donor rigidly (donor @ R.T + t) so that its H binds atom a."""
    xd, xa = donor.coords, acceptor.coords
    bond_dirs = np.array([_unit(xa[j] - xa[a]) for j in _partners(acceptor.symbols, xa, a)])
    total = bond_dirs.sum(axis=0) if len(bond_dirs) else np.zeros(3)
    if np.linalg.norm(total) > _TINY:
        axis = _unit(total)
    else:  # no bond, or bonds that cancel (linear X–A–X)
        axis = _perpendicular(bond_dirs[0], None) if len(bond_dirs) else np.array([0.0, 0, 1])
    acc = _frame(xa, xa[a], axis)
    d = _lone_pair(acc, bond_dirs, rng.uniform(*CONE_DEG), azimuth)
    heavy = min(_partners(donor.symbols, xd, h), key=lambda j: float(np.linalg.norm(xd[j] - xd[h])))
    g1 = _perpendicular(d, acc[1])
    target = np.array([d, g1, np.cross(d, g1)])
    alpha = rng.uniform(0.0, 2.0 * np.pi)
    tilt = rotation_about(np.cos(alpha) * g1 + np.sin(alpha) * target[2], sign * TILT_DEG)
    # X–H of the donor points along d, away from the acceptor, before the tilt
    rot = tilt @ target.T @ _frame(xd, xd[h], _unit(xd[heavy] - xd[h]))
    h_at = xa[a] + HBOND_SCALE * (vdw_radius("H") + vdw_radius(acceptor.symbols[a])) * d
    return rot, h_at - rot @ xd[h]


def _combined(cluster: XYZ, guest: XYZ, coords: np.ndarray) -> XYZ:
    return XYZ([*cluster.symbols, *guest.symbols], np.vstack([cluster.coords, coords]))


def _hbond_placements(cluster: XYZ, guest: XYZ, rng: np.random.Generator) -> Iterator[XYZ]:
    """Guest H → cluster acceptor pairs first, then cluster H → guest acceptor."""
    pairs = [(True, h, a) for h in _labile_hydrogens(guest) for a in _acceptor_atoms(cluster)]
    pairs += [(False, h, a) for h in _labile_hydrogens(cluster) for a in _acceptor_atoms(guest)]
    for sign, azimuth, (guest_donates, h, a) in product((1.0, -1.0), AZIMUTHS_DEG, pairs):
        for _ in range(MAX_TRIES):
            if guest_donates:
                rot, t = _hbond(guest, h, cluster, a, azimuth, sign, rng)
                moved = guest.coords @ rot.T + t
            else:  # move the guest by the inverse of the cluster's placement
                rot, t = _hbond(cluster, h, guest, a, azimuth, sign, rng)
                moved = (guest.coords - t) @ rot
            if _clear(cluster.symbols, cluster.coords, guest.symbols, moved):
                yield _combined(cluster, guest, moved)
                break


def _rigid_placements(cluster: XYZ, guest: XYZ, rng: np.random.Generator) -> Iterator[XYZ]:
    """Random orientation and direction in the cluster frame; contact at Σ r_vdw + 0.5 Å."""
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
        if _clear(cluster.symbols, cluster.coords, guest.symbols, moved):
            yield _combined(cluster, guest, moved)


def _placements(cluster: XYZ, guest: XYZ, rng: np.random.Generator) -> Iterator[XYZ]:
    found = False
    for seed in _hbond_placements(cluster, guest, rng):
        found = True
        yield seed
    if not found:
        yield from _rigid_placements(cluster, guest, rng)


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


def seeds(host: XYZ, guests: Sequence[XYZ], *, n_seeds: int = 6,
          rng_seed: int = 20260925) -> list[XYZ]:
    """Collision-free, pairwise distinct placements of host + guests (atoms in that order),
    most preferred first; empty when no collision-free seed exists."""
    rng = np.random.default_rng(rng_seed)
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

"""Generic reaction trials for single-ended discovery and relaxation discoveries (§8.2 explore).

A trial is a product-free drive (atom pairs to join / separate); the discovery engine finds the
product. Trials come in four tiers, in this order, capped per source minimum (chem 9):

1. ``polar_h``: a labile H moves to an acceptor; with two or more labile H, two-step relays
   (A ← H1–D1 ← H2–D2) come before single transfers.
2. ``h_shift``: 1,2- / 1,3-H shifts that break or form a C–H bond (H within 3 Å of the target).
3. ``heavy_bond``: at most two heavy-atom bond changes inside a fragment (3.5 Å, no aromatic
   ring bond).
4. ``association``: the closest heavy-atom pair of two fragments.

Linear molecules are bent by 10° and randomly displaced by 0.05 Å per atom (seeded), so that
NT2 does not start on a symmetry line (CH-30). Torsional isomerizations are not trials: they
are conformer pairs in ``hypotheses``.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Iterable, Sequence
from typing import Literal

import numpy as np

from hfauto.chemistry import topology
from hfauto.chemistry.vibrations import external_basis
from hfauto.core.hashing import sha256_text
from hfauto.core.records import DiscoveryRecord, MinimumRecord, ReactionTrial, SpeciesRecord

Pair = tuple[int, int]
Kind = Literal["polar_h", "h_shift", "heavy_bond", "association"]
Drive = tuple[Kind, tuple[Pair, ...], tuple[Pair, ...]]  # kind, associations, dissociations

RELAY_CONTACT_MAX_A = 3.2  # A···H1 and D1···H2 contacts of a two-step relay
H_SHIFT_MAX_A = 3.0
HEAVY_BOND_MAX_A = 3.5
LINEAR_BEND_DEG = 10.0
PERTURB_A = 0.05
RNG_SEED = 20260925
_AROMATIC_ELEMENTS = frozenset({"C", "N", "O", "S"})
_RING_SIZES = (5, 6)
_PLANAR_TOL_A = 0.1


def _pair(i: int, j: int) -> Pair:
    return (min(i, j), max(i, j))


def _neighbours(n_atoms: int, bonded: Iterable[Pair]) -> list[set[int]]:
    out: list[set[int]] = [set() for _ in range(n_atoms)]
    for i, j in bonded:
        out[i].add(j)
        out[j].add(i)
    return out


def _polar_h(symbols: Sequence[str], x: np.ndarray, nbr: list[set[int]]) -> list[Drive]:
    """Relays first (two or more labile H), then single transfers; nearest contacts first."""
    acceptors = set(topology.acceptor_atoms(symbols, x))
    donors = [(h, d) for h in topology.labile_hydrogens(symbols, x) for d in sorted(nbr[h])
              if symbols[d] != "H"]
    r = np.linalg.norm(x[:, None] - x[None], axis=-1)
    relays: list[tuple[float, Drive]] = []
    for (h1, d1), (h2, d2) in itertools.permutations(donors, 2):
        if h1 == h2 or d1 == d2 or d1 not in acceptors or r[d1, h2] > RELAY_CONTACT_MAX_A:
            continue
        for a in sorted(acceptors - {d1, d2} - nbr[h1]):
            if r[a, h1] <= RELAY_CONTACT_MAX_A:
                drive: Drive = ("polar_h", (_pair(a, h1), _pair(d1, h2)),
                                (_pair(d1, h1), _pair(d2, h2)))
                relays.append((r[a, h1] + r[d1, h2], drive))
    singles: list[tuple[float, Drive]] = [
        (r[a, h], ("polar_h", (_pair(a, h),), (_pair(d, h),)))
        for h, d in donors for a in sorted(acceptors - {d} - nbr[h])]
    return [drive for _, drive in sorted(relays) + sorted(singles)]


def _h_shifts(symbols: Sequence[str], x: np.ndarray, nbr: list[set[int]]) -> list[Drive]:
    """1,2- and 1,3-shifts of an H between heavy atoms, one of them carbon."""
    out: list[tuple[float, Drive]] = []
    for h in (i for i, s in enumerate(symbols) if s == "H"):
        for d in sorted(i for i in nbr[h] if symbols[i] != "H"):
            near = nbr[d] | {k for j in nbr[d] for k in nbr[j]}
            for t in sorted(near - {d} - nbr[h]):
                dist = float(np.linalg.norm(x[t] - x[h]))
                if symbols[t] != "H" and "C" in (symbols[d], symbols[t]) and dist <= H_SHIFT_MAX_A:
                    out.append((dist, ("h_shift", (_pair(t, h),), (_pair(d, h),))))
    return [drive for _, drive in sorted(out)]


def _rings(heavy_nbr: list[set[int]]) -> list[tuple[int, ...]]:
    """Simple cycles of 5 or 6 atoms (each once, as an ordered path)."""
    found: dict[frozenset[int], tuple[int, ...]] = {}

    def walk(path: tuple[int, ...]) -> None:
        if len(path) in _RING_SIZES and path[0] in heavy_nbr[path[-1]]:
            found.setdefault(frozenset(path), path)
        if len(path) < max(_RING_SIZES):
            for nxt in sorted(heavy_nbr[path[-1]]):
                if nxt > path[0] and nxt not in path:
                    walk((*path, nxt))

    for start in range(len(heavy_nbr)):
        walk((start,))
    return list(found.values())


def aromatic_bonds(symbols: Sequence[str], x: np.ndarray, nbr: list[set[int]]) -> set[Pair]:
    """Bonds of planar 5/6-rings of sp2-like C/N/O/S atoms (no RDKit needed)."""
    heavy_nbr = [{j for j in n if symbols[j] != "H"} if symbols[i] != "H" else set()
                 for i, n in enumerate(nbr)]
    out: set[Pair] = set()
    for ring in _rings(heavy_nbr):
        atoms = list(ring)
        centred = x[atoms] - x[atoms].mean(axis=0)
        planar = np.linalg.svd(centred, compute_uv=False)[-1] / math.sqrt(len(atoms))
        if (all(symbols[i] in _AROMATIC_ELEMENTS and len(nbr[i]) <= 3 for i in atoms)
                and planar < _PLANAR_TOL_A):
            out |= {_pair(a, b) for a, b in zip(atoms, (*atoms[1:], atoms[0]), strict=True)}
    return out


def _breakable(symbols: Sequence[str], x: np.ndarray, nbr: list[set[int]]) -> list[Pair]:
    """Heavy-atom bonds outside aromatic rings."""
    fixed = aromatic_bonds(symbols, x, nbr)
    heavy = [i for i, s in enumerate(symbols) if s != "H"]
    return sorted(_pair(i, j) for i in heavy for j in nbr[i]
                  if i < j and symbols[j] != "H" and _pair(i, j) not in fixed)


def _formations(symbols: Sequence[str], x: np.ndarray, nbr: list[set[int]],
                fragment: list[int]) -> list[Pair]:
    """Non-bonded heavy-atom pairs of one fragment within 3.5 Å, nearest first."""
    heavy = [i for i, s in enumerate(symbols) if s != "H"]
    forms = sorted((float(np.linalg.norm(x[i] - x[j])), (i, j))
                   for i, j in itertools.combinations(heavy, 2)
                   if j not in nbr[i] and fragment[i] == fragment[j])
    return [p for r, p in forms if r <= HEAVY_BOND_MAX_A]


def _heavy_bonds(symbols: Sequence[str], x: np.ndarray, nbr: list[set[int]],
                 fragment: list[int]) -> list[Drive]:
    """Formations, formation + breaking at a shared atom, breakings."""
    breakable = _breakable(symbols, x, nbr)
    forms = _formations(symbols, x, nbr, fragment)
    joins: list[Drive] = [("heavy_bond", (p,), ()) for p in forms]
    swaps: list[Drive] = [("heavy_bond", (p,), (b,))
                          for p in forms for b in breakable if set(p) & set(b)]
    cuts: list[Drive] = [("heavy_bond", (), (b,)) for b in breakable]
    return joins + swaps + cuts


def _associations(symbols: Sequence[str], x: np.ndarray,
                  groups: Sequence[Sequence[int]]) -> list[Drive]:
    out: list[tuple[float, Drive]] = []
    for fa, fb in itertools.combinations(groups, 2):
        pairs = [(i, j) for i in fa for j in fb if symbols[i] != "H" and symbols[j] != "H"]
        pairs = pairs or [(i, j) for i in fa for j in fb]
        dist, i, j = min((float(np.linalg.norm(x[i] - x[j])), i, j) for i, j in pairs)
        out.append((dist, ("association", (_pair(i, j),), ())))
    return [drive for _, drive in sorted(out)]


def is_linear(symbols: Sequence[str], coords: np.ndarray) -> bool:
    return len(symbols) >= 3 and external_basis(symbols, coords).shape[1] == 5


def perturb_linear(coords: np.ndarray, *, seed: int = RNG_SEED) -> np.ndarray:
    """Bend the atoms on one side of the middle atom by 10°, then displace every atom 0.05 Å."""
    x = np.asarray(coords, dtype=float).reshape(-1, 3)
    axis = np.linalg.svd(x - x.mean(axis=0))[2][0]
    along = (x - x.mean(axis=0)) @ axis
    pivot = int(np.argsort(along)[len(along) // 2])
    normal = np.cross(axis, np.eye(3)[int(np.argmin(np.abs(axis)))])
    normal /= np.linalg.norm(normal)
    angle = math.radians(LINEAR_BEND_DEG)
    k = np.array([[0, -normal[2], normal[1]], [normal[2], 0, -normal[0]],
                  [-normal[1], normal[0], 0]])
    rotation = np.eye(3) + math.sin(angle) * k + (1 - math.cos(angle)) * k @ k
    bent = x.copy()
    side = along < along[pivot]
    bent[side] = (x[side] - x[pivot]) @ rotation.T + x[pivot]
    step = np.random.default_rng(seed).normal(size=x.shape)
    return bent + PERTURB_A * step / np.linalg.norm(step, axis=1, keepdims=True)


def _trial_id(source_minimum: str, drive: Drive) -> str:
    kind, assoc, dissoc = drive
    return "trial_" + sha256_text(f"{source_minimum}|{kind}|{assoc}|{dissoc}")


def generate(source_minimum: str, symbols: Sequence[str], coords: np.ndarray, *,
             max_trials: int = 10) -> tuple[np.ndarray, list[ReactionTrial]]:
    """(start coordinates, trials): the start is perturbed for linear molecules.

    Every trial starts with NT2; the explore stage falls back to AFIR on the same drive.
    """
    x = np.asarray(coords, dtype=float).reshape(-1, 3)
    bonded = topology.bonds(symbols, x)
    nbr = _neighbours(len(symbols), bonded)
    groups = topology.fragments(symbols, x)
    fragment = [0] * len(symbols)
    for index, group in enumerate(groups):
        for atom in group:
            fragment[atom] = index
    drives = [*_polar_h(symbols, x, nbr), *_h_shifts(symbols, x, nbr),
              *_heavy_bonds(symbols, x, nbr, fragment), *_associations(symbols, x, groups)]
    unique: dict[tuple[frozenset[Pair], frozenset[Pair]], Drive] = {}
    for drive in drives:
        unique.setdefault((frozenset(drive[1]), frozenset(drive[2])), drive)
    perturbed = is_linear(symbols, x)
    trials = [
        ReactionTrial(trial_id=_trial_id(source_minimum, d), source_minimum=source_minimum,
                      kind=d[0], mechanism="nt2", associations=d[1], dissociations=d[2],
                      perturbed=perturbed)
        for d in list(unique.values())[:max_trials]
    ]
    return (perturb_linear(x) if perturbed else x), trials


def product_verdict(mechanism: str, *, ts_validated: bool, barrier_kj: float | None,
                    reaction_kj: float | None, barrier_max_kj: float = 150.0,
                    reaction_max_kj: float = 100.0) -> str | None:
    """None when a reported product is kept, else the negative reason (CH-28).

    An NT2 product needs its frequency-validated TS and IRC; a reaction energy that could not
    be evaluated is outside the window (fail-closed). AFIR products have no barrier.
    """
    if mechanism == "nt2" and not ts_validated:
        return "ts_not_validated"
    if reaction_kj is None or reaction_kj > reaction_max_kj:
        return "out_of_window"
    if barrier_kj is not None and barrier_kj > barrier_max_kj:
        return "out_of_window"
    return None


def relaxation_discoveries(species: Iterable[SpeciesRecord],
                           minima: Iterable[MinimumRecord]) -> list[DiscoveryRecord]:
    """A seed whose state label differs from its basin's collapsed without a barrier."""
    labels = {s.species_id: s.state_label for s in species}
    out: list[DiscoveryRecord] = []
    for minimum in minima:
        for species_id in dict.fromkeys((minimum.species_id, *minimum.members)):
            label = labels.get(species_id)
            if label is None or label == minimum.state_label:
                continue
            out.append(DiscoveryRecord(
                discovery_id=f"relax_{minimum.minimum_id}_{species_id}",
                source_minimum=minimum.minimum_id, mechanism="relaxation", outcome="negative",
                reason=f"collapsed_to:{minimum.state_label}"))
    return out

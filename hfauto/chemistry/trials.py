"""Reaction trials for single-ended discovery (§8.2 explore).

A trial is a product-free drive (bonds to form and to break); the discovery engine finds the
product. One element-independent enumerator uses formation candidates (a, b): not bonded,
r_ab ≤ Σr_vdW, at least three bonds apart (inside or across fragments).

- ``transfer`` (T1): form a–b, break b–c (c ≠ a): H transfer, 1,2- / 1,3-shift, SN2, H
  abstraction, ligand exchange. a–b may be two bonds apart only when c bridges them (1,2-shift).
- ``relay`` (T2): two chained H transfers, the second H moving into the atom the first one left
  (an exchange, the second H leaving the first acceptor, needs no contact of its own).
- ``formation`` (T3): form a–b only (addition, association, ring closure).
- ``dissociation`` (T4): break one bond; only for open-shell or charged sources (a closed-shell
  neutral homolysis is not a restricted Kohn–Sham reaction).

An atom that gains bonds may not end above ``elements.max_coordination`` bonded neighbours.
Drives with the same WL atom classes of their formed and broken bonds and the same formed-pair
distances (0.1 Å) are one class, one trial. A drive's value is r_ab/Σr_cov: a transfer's then
its a–b–c angle (closest to linear first), a relay's its farther contact then its nearer one, a
dissociation's the most stretched bond first. 1,2-shifts come after the through-space contacts:
their r_ab is set by the bond angle, not by how the source is arranged. A class's trial is its
drive of least value, then of least canonical rank (each atom's WL class and its sorted
distances to the atoms of every class). The classes of a template go by their trial's value,
then their description; the templates take turns T1 → T2 → T3 → T4 up to the cap, and classes
of one template and value that the cap would split are all left out. A trial's id is its class
description and its source. So the trials depend on the structure, not on the atom numbering.
Linear molecules are bent by 10° and displaced by 0.05 Å per atom (seeded), so that NT2 does
not start on a symmetry line. A torsion is studied only when declared (``hypotheses``).
"""

from __future__ import annotations

import itertools
from collections import Counter, defaultdict
from collections.abc import Sequence
from typing import Literal

import numpy as np
from scipy.sparse.csgraph import shortest_path

from hfauto.chemistry import topology
from hfauto.chemistry.elements import covalent_radius, max_coordination, vdw_radius
from hfauto.chemistry.geometry import rotation_about
from hfauto.chemistry.vibrations import external_basis
from hfauto.core.hashing import sha256_text
from hfauto.core.records import ReactionTrial

Pair = tuple[int, int]
Kind = Literal["transfer", "relay", "formation", "dissociation"]
Value = tuple[float, ...]  # the order within a template
Drive = tuple[Value, tuple[Pair, ...], tuple[Pair, ...]]  # value, form, break
Key = tuple[tuple[Pair, ...], tuple[Pair, ...], tuple[float, ...]]  # the class description
Print = tuple[int, tuple[tuple[float, ...], ...]]  # WL class, sorted distances per class

LINEAR_BEND_DEG = 10.0
PERTURB_A = 0.05
RNG_SEED = 20260925


def _pair(i: int, j: int) -> Pair:
    return (min(i, j), max(i, j))


class _Source:
    """Bond graph, contacts and free valences of one source structure."""

    def __init__(self, symbols: Sequence[str], x: np.ndarray) -> None:
        self.symbols, self.x = symbols, x
        self.bonded = topology.bonds(symbols, x)
        adjacency = np.zeros((len(symbols), len(symbols)), dtype=bool)
        for i, j in self.bonded:
            adjacency[i, j] = adjacency[j, i] = True
        self.nbr = [set(np.flatnonzero(row).tolist()) for row in adjacency]
        self.hops = shortest_path(adjacency, unweighted=True)  # inf across fragments
        self.r = np.linalg.norm(x[:, None] - x[None], axis=-1)
        cov, vdw = (np.array([f(s) for s in symbols]) for f in (covalent_radius, vdw_radius))
        self.ratio = self.r / (cov[:, None] + cov[None])
        self.contact = (self.r <= vdw[:, None] + vdw[None]) & (self.hops >= 2)
        self.room = [max_coordination(s) - len(self.nbr[i]) for i, s in enumerate(symbols)]
        self.classes = topology.wl_classes(symbols, self.bonded)
        members = [np.equal(self.classes, c) for c in sorted(set(self.classes))]
        self.prints: list[Print] = [
            (c, tuple(tuple(np.sort(self.r[i, m]).tolist()) for m in members))
            for i, c in enumerate(self.classes)]

    def fits(self, formed: Sequence[Pair], broken: Sequence[Pair]) -> bool:
        """No atom that gains bonds ends above its maximum coordination."""
        net = Counter(i for pair in formed for i in pair)
        net.subtract(i for pair in broken for i in pair)
        return all(gain <= self.room[i] for i, gain in net.items() if gain > 0)

    def _cos(self, a: int, b: int, c: int) -> float:
        u, v = self.x[a] - self.x[b], self.x[c] - self.x[b]
        return float(u @ v / (np.linalg.norm(u) * np.linalg.norm(v)))

    def steps(self) -> list[tuple[Value, float, int, int, int]]:
        """((1,2-shift, r_ab/Σr_cov), cos a–b–c, a, b, c): form a–b and break b–c, valence not
        checked."""
        return [((float(self.hops[a, b] == 2), float(self.ratio[a, b])), self._cos(a, b, c),
                 a, b, c)
                for a, b in _pairs(self.contact) for c in self.nbr[b] - {a}
                if self.hops[a, b] >= 3 or c in self.nbr[a]]

    def relays(self, h_steps: list[tuple[Value, float, int, int, int]]) -> list[Drive]:
        """h1 moves c -> a, then h2 moves d -> c. An exchange (d = a) needs no second contact:
        the first transfer already holds a and c together (NH3·HF double H exchange)."""
        chains = [(s1, a, h1, c, s2, h2, d) for s1, _, a, h1, c in h_steps
                  for s2, _, into, h2, d in h_steps if into == c and d != a]
        chains += [(s1, a, h1, c, (0.0, float(self.ratio[c, h2])), h2, a)
                   for s1, _, a, h1, c in h_steps for h2 in self.nbr[a]
                   if self.symbols[h2] == "H" and self.hops[c, h2] >= 3]
        return [(max(s1, s2) + min(s1, s2), (_pair(a, h1), _pair(c, h2)),
                 (_pair(h1, c), _pair(h2, d))) for s1, a, h1, c, s2, h2, d in chains]

    def templates(self, open_or_charged: bool) -> dict[Kind, list[Drive]]:
        """Each template's drives; drives that break the valence rule are left out."""
        steps = self.steps()
        formations = _pairs(np.triu(self.contact & (self.hops >= 3)))
        found: dict[Kind, list[Drive]] = {
            "transfer": [((*s, cos), (_pair(a, b),), (_pair(b, c),)) for s, cos, a, b, c in steps],
            "relay": self.relays([step for step in steps if self.symbols[step[3]] == "H"]),
            "formation": [((float(self.ratio[p]),), (p,), ()) for p in formations],
            "dissociation": [((-float(self.ratio[p]),), (), (p,))
                             for p in (self.bonded if open_or_charged else ())],
        }
        return {kind: [d for d in ds if self.fits(d[1], d[2])] for kind, ds in found.items()}

    def key(self, drive: Drive) -> Key:
        """WL class pairs of the formed and broken bonds and the formed-pair distances (0.1 Å)."""
        _, formed, broken = drive
        f, b = (tuple(sorted(_pair(self.classes[i], self.classes[j]) for i, j in bonds))
                for bonds in (formed, broken))
        return f, b, tuple(sorted(round(float(self.r[i, j]), 1) for i, j in formed))

    def rank(self, drive: Drive) -> tuple[tuple[tuple[Print, ...], ...], ...]:
        """Canonical rank: the formed and the broken bonds as sorted pairs of atom prints."""
        return tuple(tuple(sorted(tuple(sorted((self.prints[i], self.prints[j])))
                                  for i, j in bonds)) for bonds in drive[1:])


def _pairs(mask: np.ndarray) -> list[Pair]:
    return [(int(i), int(j)) for i, j in zip(*np.nonzero(mask), strict=True)]


def _classes(symbols: Sequence[str], x: np.ndarray, open_or_charged: bool,
             max_trials: int) -> list[tuple[Kind, Key, Drive]]:
    """(kind, class description, its trial's drive) up to the cap, the templates taking turns."""
    source = _Source(symbols, x)
    turns: list[list[tuple[Value, Key, Kind, Drive]]] = []
    for kind, template in source.templates(open_or_charged).items():
        drives: defaultdict[Key, list[Drive]] = defaultdict(list)
        for drive in template:
            drives[source.key(drive)].append(drive)
        best = {key: min(ds, key=lambda d: (d[0], source.rank(d))) for key, ds in drives.items()}
        turns.append(sorted(((d[0], key, kind, d) for key, d in best.items()),
                            key=lambda c: c[:2]))
    order = [c for turn in itertools.zip_longest(*turns) for c in turn if c is not None]
    cut = {(kind, value) for value, _, kind, _ in order[max_trials:]}
    return [(kind, key, drive) for value, key, kind, drive in order[:max_trials]
            if (kind, value) not in cut]


def is_linear(symbols: Sequence[str], coords: np.ndarray) -> bool:
    return len(symbols) >= 3 and external_basis(symbols, coords).shape[1] == 5


def perturb_linear(coords: np.ndarray) -> np.ndarray:
    """Bend the atoms on one side of the middle atom by 10°, then displace every atom 0.05 Å."""
    x = np.asarray(coords, dtype=float).reshape(-1, 3)
    axis = np.linalg.svd(x - x.mean(axis=0))[2][0]
    along = (x - x.mean(axis=0)) @ axis
    pivot = int(np.argsort(along)[len(along) // 2])
    normal = np.cross(axis, np.eye(3)[int(np.argmin(np.abs(axis)))])
    rotation = rotation_about(normal / np.linalg.norm(normal), LINEAR_BEND_DEG)
    bent = x.copy()
    side = along < along[pivot]
    bent[side] = (x[side] - x[pivot]) @ rotation.T + x[pivot]
    step = np.random.default_rng(RNG_SEED).normal(size=x.shape)
    return bent + PERTURB_A * step / np.linalg.norm(step, axis=1, keepdims=True)


def generate(source_minimum: str, symbols: Sequence[str], coords: np.ndarray, *,
             charge: int = 0, multiplicity: int = 1,
             max_trials: int = 10) -> tuple[np.ndarray, list[ReactionTrial]]:
    """(start coordinates, NT2 trials): a linear molecule starts bent."""
    x = np.asarray(coords, dtype=float).reshape(-1, 3)
    trials = [
        ReactionTrial(trial_id="trial_" + sha256_text(f"{source_minimum}|{kind}|{key}"),
                      source_minimum=source_minimum, kind=kind, associations=drive[1],
                      dissociations=drive[2])
        for kind, key, drive in _classes(symbols, x, multiplicity > 1 or charge != 0, max_trials)
    ]
    return (perturb_linear(x) if is_linear(symbols, x) else x), trials

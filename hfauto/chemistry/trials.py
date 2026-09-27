"""Reaction trials for single-ended discovery and relaxation discoveries (§8.2 explore).

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
The templates take turns T1 → T2 → T3 → T4 up to the cap, each ordered by r_ab/Σr_cov (a relay
by its farther contact, several c of one transfer by the a–b–c angle closest to linear, a
dissociation by the most stretched bond). 1,2-shifts come after the through-space contacts:
their r_ab is set by the bond angle, not by how the source is arranged. Drives with the same
WL atom classes and formed-pair distances (0.1 Å) are one trial. Linear molecules are bent by
10° and displaced by 0.05 Å per atom (seeded), so that NT2 does not start on a symmetry line
(CH-30). A torsion is studied only when declared (``hypotheses``).
"""

from __future__ import annotations

import itertools
import math
from collections import Counter
from collections.abc import Hashable, Iterable, Sequence
from typing import Literal

import numpy as np
from scipy.sparse.csgraph import shortest_path

from hfauto.chemistry import topology
from hfauto.chemistry.elements import covalent_radius, max_coordination, vdw_radius
from hfauto.chemistry.vibrations import external_basis
from hfauto.core.hashing import sha256_text
from hfauto.core.records import DiscoveryRecord, MinimumRecord, ReactionTrial, SpeciesRecord

Pair = tuple[int, int]
Kind = Literal["transfer", "relay", "formation", "dissociation"]
Order = tuple[bool, float]  # (1,2-shift, r/Σr_cov) within a template
Drive = tuple[Order, tuple[Pair, ...], tuple[Pair, ...]]  # order, form, break
Trial = tuple[Kind, tuple[Pair, ...], tuple[Pair, ...]]

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

    def fits(self, formed: Sequence[Pair], broken: Sequence[Pair]) -> bool:
        """No atom that gains bonds ends above its maximum coordination."""
        net = Counter(i for pair in formed for i in pair)
        net.subtract(i for pair in broken for i in pair)
        return all(gain <= self.room[i] for i, gain in net.items() if gain > 0)

    def _cos(self, a: int, b: int, c: int) -> float:
        u, v = self.x[a] - self.x[b], self.x[c] - self.x[b]
        return float(u @ v / (np.linalg.norm(u) * np.linalg.norm(v)))

    def steps(self) -> list[tuple[Order, int, int, int]]:
        """(order, a, b, c): form a–b and break b–c, valence not checked. The c of one (a, b)
        come by the a–b–c angle, closest to linear first."""
        found = sorted(((bool(self.hops[a, b] == 2), float(self.ratio[a, b])),
                        self._cos(a, b, c), a, b, c)
                       for a, b in _pairs(self.contact) for c in sorted(self.nbr[b] - {a})
                       if self.hops[a, b] >= 3 or c in self.nbr[a])
        return [(s, a, b, c) for s, _, a, b, c in found]

    def relays(self, h_steps: list[tuple[Order, int, int, int]]) -> list[Drive]:
        """h1 moves c -> a, then h2 moves d -> c. An exchange (d = a) needs no second contact:
        the first transfer already holds a and c together (NH3·HF double H exchange)."""
        chains = [(s1, a, h1, c, s2, h2, d) for s1, a, h1, c in h_steps
                  for s2, into, h2, d in h_steps if into == c and d != a]
        chains += [(s1, a, h1, c, (False, float(self.ratio[c, h2])), h2, a)
                   for s1, a, h1, c in h_steps for h2 in sorted(self.nbr[a])
                   if self.symbols[h2] == "H" and self.hops[c, h2] >= 3]
        return [(max(s1, s2), (_pair(a, h1), _pair(c, h2)), (_pair(h1, c), _pair(h2, d)))
                for s1, a, h1, c, s2, h2, d in chains]

    def templates(self, open_or_charged: bool) -> dict[Kind, list[Drive]]:
        """Each template's drives in order; drives that break the valence rule are left out."""
        steps = self.steps()
        relays = self.relays([step for step in steps if self.symbols[step[2]] == "H"])
        formations = _pairs(np.triu(self.contact & (self.hops >= 3)))
        cuts = sorted(self.bonded) if open_or_charged else []
        found: dict[Kind, list[Drive]] = {
            "transfer": [(s, (_pair(a, b),), (_pair(b, c),)) for s, a, b, c in steps],
            "relay": sorted(relays),
            "formation": sorted(((False, float(self.ratio[p])), (p,), ()) for p in formations),
            "dissociation": sorted(((False, -float(self.ratio[p])), (), (p,)) for p in cuts),
        }
        return {kind: [d for d in ds if self.fits(d[1], d[2])] for kind, ds in found.items()}

    def key(self, drive: Drive) -> Hashable:
        """WL class pairs of the formed and broken bonds and the formed-pair distances (0.1 Å)."""
        _, formed, broken = drive
        classes = [frozenset(tuple(sorted((self.classes[i], self.classes[j]))) for i, j in bonds)
                   for bonds in (formed, broken)]
        return (*classes, tuple(sorted(round(float(self.r[i, j]), 1) for i, j in formed)))


def _pairs(mask: np.ndarray) -> list[Pair]:
    return [(int(i), int(j)) for i, j in zip(*np.nonzero(mask), strict=True)]


def _drives(symbols: Sequence[str], x: np.ndarray, open_or_charged: bool) -> list[Trial]:
    """The distinct drives of one source, the templates taking turns."""
    source = _Source(symbols, x)
    seen: set[Hashable] = set()
    turns: list[list[Trial]] = []
    for kind, template in source.templates(open_or_charged).items():
        turns.append([])
        for drive in template:
            if (key := source.key(drive)) not in seen:
                seen.add(key)
                turns[-1].append((kind, drive[1], drive[2]))
    return [d for turn in itertools.zip_longest(*turns) for d in turn if d is not None]


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


def generate(source_minimum: str, symbols: Sequence[str], coords: np.ndarray, *,
             charge: int = 0, multiplicity: int = 1,
             max_trials: int = 10) -> tuple[np.ndarray, list[ReactionTrial]]:
    """(start coordinates, trials): the start is perturbed for linear molecules.

    Every trial starts with NT2; the explore stage falls back to AFIR on the same drive.
    """
    x = np.asarray(coords, dtype=float).reshape(-1, 3)
    perturbed = is_linear(symbols, x)
    trials = [
        ReactionTrial(
            trial_id="trial_" + sha256_text(f"{source_minimum}|{kind}|{form}|{cut}"),
            source_minimum=source_minimum, kind=kind, mechanism="nt2", associations=form,
            dissociations=cut, perturbed=perturbed)
        for kind, form, cut in _drives(symbols, x, multiplicity > 1 or charge != 0)[:max_trials]
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

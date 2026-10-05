"""Which structures the minima stage refines (design §8.2, minima(dft)).

The DFT entry of the low-level edges (``admit``) is decided once, on DFT//xTB single points:
an edge is admitted when the highest point on its way from a start state (its TS, else its
higher end) lies within the reaction window of that start, in order of that height and up to a
count per composition. The candidates are the screen minima and the end species no screen
basin holds, of the reacting compositions only (``reacting_candidates``). Per (composition,
state label) at most ``per_state`` structures within ``window_kcal`` of the group's lowest
energy, plus every always-kept structure (the declared endpoints and the ends of the admitted
edges); each group keeps its lowest structure (the Curtin-Hammett reference). ``rerank``
re-scores by single points only the groups with more than ``per_state`` candidates: a smaller
group is kept whole.
"""

from __future__ import annotations

import heapq
import math
from collections import Counter, defaultdict
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from typing import cast

from hfauto.chemistry.thermo import declared_monomers
from hfauto.chemistry.xyz import hill_formula
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.records import DiscoveryRecord, MinimumRecord, SpeciesRecord
from hfauto.core.system import SystemConfig

Point = float | str  # a DFT//xTB single point (Eh), or the kind of its failure
Verdict = tuple[str | None, float | None]  # (why not admitted, height above the start kcal/mol)


@dataclass(frozen=True)
class Edge:
    """A low-level edge (a product discovery) between two states (``source``, ``product``: any
    keys that equal for one state). ``points`` takes its single points when first asked: at its
    two end structures and, when it has a TS, the TS. ``first``: its source is a start state
    (explore's first generation). A state's energy is the lowest of its end structures' asked
    so far (its conformers interconvert below any edge: Curtin-Hammett), so an edge leaving a
    state is reached from whichever conformer another edge reached it at."""

    discovery_id: str
    composition_id: str
    source: str
    product: str
    points: Callable[[], tuple[Point, ...]]
    first: bool = True


@dataclass(frozen=True)
class Candidate:
    species_id: str
    composition_id: str
    state_label: str
    energy_hartree: float | None  # screen (or single-point) energy; None is unranked
    always: bool = False  # declared endpoints and admitted edge ends are refined regardless


class _Walk:
    """Best-first ways from each start state over the edges, asking an edge's single points
    only when its source is reached within the window: the height only grows along a way, so an
    edge beyond a source reached above the window is above it too (S10: 850 edges after the
    budget, most of them far above)."""

    def __init__(self, edges: Sequence[Edge], window: float) -> None:
        self.leaving: defaultdict[str, list[Edge]] = defaultdict(list)
        for e in edges:
            self.leaving[e.source].append(e)
        self.window = window
        self.asked: dict[str, tuple[Point, ...]] = {}
        self.energy: dict[str, float] = {}
        self.best: dict[str, tuple[float, float, int]] = {}  # id -> (height, span Eh, steps)
        self.failed: dict[str, str] = {}
        self.beyond: set[str] = set()  # states reached above the window only

    def points(self, e: Edge) -> tuple[Point, ...]:
        if e.discovery_id not in self.asked:
            self.asked[e.discovery_id] = pts = e.points()
            for state, p in zip((e.source, e.product), pts, strict=False):
                if not isinstance(p, str):
                    self.energy[state] = min(p, self.energy.get(state, math.inf))
            if any(isinstance(p, str) for p in pts):
                self.failed[e.discovery_id] = next(p for p in pts if isinstance(p, str))
        return self.asked[e.discovery_id]

    def run(self, start: str) -> None:
        """One start: states in order of (height, span, edges from the start), each expanded
        once while its height is within the window; an edge keeps its lowest such way over the
        starts."""
        base = self.energy[start]
        heap: list[tuple[float, float, int, str, float]] = [(0.0, 0.0, 0, start, base)]
        done: set[str] = set()
        while heap:
            height, span, steps, state, low = heapq.heappop(heap)
            if state in done:
                continue
            done.add(state)
            if height * HARTREE_TO_KCAL_MOL > self.window:
                self.beyond.add(state)
                continue
            for e in self.leaving[state]:
                pts = self.points(e)
                if e.discovery_id in self.failed:
                    continue
                top = max(cast(float, p) for p in pts[1:])  # the product and the TS
                way = (max(height, top - base), max(span, top - low), steps + 1)
                self.best[e.discovery_id] = min(self.best.get(e.discovery_id, way), way)
                heapq.heappush(heap, (*way, e.product, min(low, self.energy[e.product])))


def _downstream(leaving: Mapping[str, Sequence[Edge]], states: Iterable[str]) -> set[str]:
    """The states and every state their edges lead to."""
    seen, todo = set(states), list(states)
    while todo:
        for e in leaving.get(todo.pop(), ()):
            if e.product not in seen:
                seen.add(e.product)
                todo.append(e.product)
    return seen


def _heights(edges: Sequence[Edge], window: float
             ) -> tuple[dict[str, Verdict], dict[str, tuple[float, int]]]:
    """Discovery id -> (None, height kcal/mol), ``out_of_window`` (its source reached above
    the window only, or only beyond such a state: never asked) or ``not_evaluated:<why>``
    (``unreached``: no evaluated way, e.g. only past a failed single point); and, for each judged
    edge, the energetic span (kcal/mol) and the number of edges of the way of its height."""
    walk = _Walk(edges, window)
    starts = sorted({e.source for e in edges if e.first})
    for e in (e for s in starts for e in walk.leaving[s]):
        walk.points(e)  # a start's energy: its structures' lowest
    for s in starts:
        if s in walk.energy:
            walk.run(s)
    failed = {e.product for e in edges if e.discovery_id in walk.failed}
    beyond = _downstream(walk.leaving, walk.beyond) - _downstream(walk.leaving, failed)
    out: dict[str, Verdict] = {}
    order: dict[str, tuple[float, int]] = {}
    for e in edges:
        if e.discovery_id in walk.failed:
            out[e.discovery_id] = (f"not_evaluated:{walk.failed[e.discovery_id]}", None)
        elif e.discovery_id in walk.best:
            height, span, steps = walk.best[e.discovery_id]
            out[e.discovery_id] = (None, height * HARTREE_TO_KCAL_MOL)
            order[e.discovery_id] = (span * HARTREE_TO_KCAL_MOL, steps)
        elif e.source in beyond:
            out[e.discovery_id] = ("out_of_window", None)
        else:
            out[e.discovery_id] = ("not_evaluated:unreached", None)
    return out, order


def admit(edges: Sequence[Edge], *, window_kcal: float, per_composition: int
          ) -> dict[str, Verdict]:
    """Discovery id -> (None when admitted, else why not; its height above its start).

    The height of an edge is its highest point (its TS, else its higher end; its product too,
    so the reaction energy is bounded with it) on its lowest way from a start state, above that
    start, for the start that gives the lowest height. ``not_evaluated:<failure>``: a single
    point failed, or ``unreached``: no evaluated way leads to its source. In order of height,
    then of energetic span on that way (S5: every edge below a high start has height 0, and a
    deep well passed on the way must be climbed out of again), then of its number of edges (an
    edge before those it leads to), an edge above ``window_kcal`` is
    ``out_of_window`` (also, with no height, one whose source lies beyond the window: its
    single points are never taken) and one beyond ``per_composition`` pairs of ends
    ``over_cap``."""
    out, order = _heights(edges, window_kcal)
    pairs: defaultdict[str, set[tuple[str, str]]] = defaultdict(set)
    judged = {e.discovery_id: h for e in edges if (h := out[e.discovery_id][1]) is not None}
    for e in sorted((e for e in edges if e.discovery_id in judged),
                    key=lambda e: (judged[e.discovery_id], *order[e.discovery_id],
                                   e.discovery_id)):
        kept, pair, height = pairs[e.composition_id], (e.source, e.product), judged[e.discovery_id]
        if height > window_kcal:
            out[e.discovery_id] = ("out_of_window", height)
        elif pair not in kept and len(kept) >= per_composition:
            out[e.discovery_id] = ("over_cap", height)
        else:
            kept.add(pair)
    return out


def _state(c: Candidate) -> tuple[str, str]:
    return c.composition_id, c.state_label


def _reacting(candidates: Iterable[Candidate], species: Mapping[str, SpeciesRecord],
              system: SystemConfig) -> set[str]:
    """Compositions of the always-kept candidates, plus the monomers of those complexes: the
    association and separated references of the thermo stage (thermo.declared_monomers)."""
    reacting = {c.composition_id for c in candidates if c.always}
    monomers = declared_monomers(species.values(), system.compositions)
    formulas = {(hill_formula(s.geometry.symbols), s.charge) for s in species.values()
                if s.composition_id in reacting}
    return reacting | {composition for f in formulas for parts in monomers.get(f, ())
                       for composition, _ in parts}


def reacting_candidates(
    species: Mapping[str, SpeciesRecord], admitted: Iterable[DiscoveryRecord],
    screen: Sequence[MinimumRecord], system: SystemConfig,
) -> tuple[list[Candidate], list[str]]:
    """The candidates of the reacting compositions, and the other compositions, sorted.

    The species of each of ``screen``, the screen minima of ``species``, stands for its basin,
    always kept when the basin holds a declared endpoint or an end of an ``admitted`` edge;
    such a species that no screen basin holds stands for itself, always kept. One candidate
    per species id."""
    keep = {s.id for s in system.species if s.role == "endpoint"}
    keep |= {s for d in admitted for s in (d.source_species, d.product_species) if s}
    pool = {m.species_id: Candidate(m.species_id, m.composition_id, m.state_label,
                                    m.energy_hartree,
                                    always=bool(keep & {m.species_id, *m.members}))
            for m in screen}
    held = {s for m in screen for s in (m.species_id, *m.members)}
    for s in (species[sid] for sid in sorted((keep & species.keys()) - held)):
        pool[s.species_id] = Candidate(s.species_id, s.composition_id, s.state_label,
                                       s.energy_hartree, always=True)
    reacting = _reacting(pool.values(), species, system)
    return ([c for c in pool.values() if c.composition_id in reacting],
            sorted({c.composition_id for c in pool.values()} - reacting))


def select_for_refinement(
    candidates: Sequence[Candidate], *, per_state: int = 3, window_kcal: float = 6.0
) -> list[Candidate]:
    """The selected candidates in input order; unranked ones only when always kept."""
    chosen = {c.species_id for c in candidates if c.always}
    states: defaultdict[tuple[str, str], list[Candidate]] = defaultdict(list)
    for c in candidates:
        states[_state(c)].append(c)
    for group in states.values():
        ranked = sorted(((c.energy_hartree, c) for c in group if c.energy_hartree is not None),
                        key=lambda pair: (pair[0], pair[1].species_id))
        chosen.update(
            c.species_id for energy, c in ranked[:max(per_state, 1)]
            if (energy - ranked[0][0]) * HARTREE_TO_KCAL_MOL <= window_kcal
        )
    return [c for c in candidates if c.species_id in chosen]


def rerank(
    candidates: Sequence[Candidate],
    single_points: Callable[[list[Candidate]], Mapping[str, float]],
    keep_per_state: int,
    window_kcal: float,
) -> list[Candidate]:
    """Re-score with ``single_points`` (species id -> energy) only the candidates, not always
    kept, of the groups that hold more than ``keep_per_state`` of them (the only ones whose
    selection a single point can change), and keep, per state, the ``keep_per_state`` lowest
    within ``window_kcal`` of the lowest (the better estimate cuts the window again); every
    other candidate stays.

    A crowded candidate without a single point drops out.
    """
    count = Counter(_state(c) for c in candidates if not c.always)
    busy = [c for c in candidates if not c.always and count[_state(c)] > keep_per_state]
    sp_energies = single_points(busy)
    rescored = [replace(c, energy_hartree=sp_energies.get(c.species_id)) for c in busy]
    kept = select_for_refinement(rescored, per_state=keep_per_state, window_kcal=window_kcal)
    dropped = {c.species_id for c in busy} - {c.species_id for c in kept}
    return [c for c in candidates if c.species_id not in dropped]

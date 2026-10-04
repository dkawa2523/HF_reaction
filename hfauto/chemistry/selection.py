"""Which structures the minima stage refines (design §8.2, minima(dft)).

The candidates are the screen minima, the discovery ends no screen basin holds and the
relaxation seeds, of the reacting compositions only (``reacting_candidates``). Per (composition,
state label) at most ``per_state`` structures within ``window_kcal`` of the group's lowest
energy, plus every always-kept structure (discovery endpoints, relaxation seeds); each group
keeps its lowest structure (the Curtin-Hammett reference). ``rerank`` re-scores by single points
only the groups with more than ``per_state`` candidates: a smaller group is kept whole.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, replace

from hfauto.chemistry.thermo import declared_monomers
from hfauto.chemistry.xyz import hill_formula
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.records import DiscoveryRecord, MinimumRecord, SpeciesRecord
from hfauto.core.system import SystemConfig


@dataclass(frozen=True)
class Candidate:
    species_id: str
    composition_id: str
    state_label: str
    energy_hartree: float | None  # screen (or single-point) energy; None is unranked
    always: bool = False  # discovery endpoints are refined regardless of energy


def _state(c: Candidate) -> tuple[str, str]:
    return c.composition_id, c.state_label


def _reacting(candidates: Iterable[Candidate], species: Mapping[str, SpeciesRecord],
              system: SystemConfig) -> set[str]:
    """Compositions of the declared endpoints and of the always-kept candidates, plus the
    monomers of those complexes: the association and separated references of the thermo stage
    (thermo.declared_monomers)."""
    ends = {species[s.id].composition_id for s in system.species
            if s.role == "endpoint" and s.id in species}
    reacting = ends | {c.composition_id for c in candidates if c.always}
    monomers = declared_monomers(species.values(), system.compositions)
    formulas = {(hill_formula(s.geometry.symbols), s.charge) for s in species.values()
                if s.composition_id in reacting}
    return reacting | {composition for f in formulas for parts in monomers.get(f, ())
                       for composition, _ in parts}


def reacting_candidates(
    species: Mapping[str, SpeciesRecord], discoveries: Iterable[DiscoveryRecord],
    screen: Sequence[MinimumRecord], seeds: Iterable[SpeciesRecord], system: SystemConfig,
) -> tuple[list[Candidate], list[str]]:
    """The candidates of the reacting compositions, and the other compositions, sorted.

    The species of each of ``screen``, the screen minima of ``species``, stands for its basin,
    always kept when the basin holds a discovery source or product; a product that skipped
    screen and each of ``seeds``, the relaxation seeds as their own species
    (hypotheses.seed_species_id), stand for themselves, always kept. One candidate per species
    id: a later one replaces it."""
    found = [d for d in discoveries if d.outcome == "product"]
    always = {d.source_minimum for d in found}
    always |= {d.product_species for d in found if d.product_species is not None}
    pool = {m.species_id: Candidate(m.species_id, m.composition_id, m.state_label,
                                    m.energy_hartree,
                                    always=bool(always & {m.minimum_id, m.species_id, *m.members}))
            for m in screen}
    held = {s for m in screen for s in (m.species_id, *m.members)}
    for s in [*(species[sid] for sid in sorted((always & species.keys()) - held)), *seeds]:
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

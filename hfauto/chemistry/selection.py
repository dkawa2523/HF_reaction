"""Which structures the minima stage refines (design §8.2, minima(dft)).

Per (composition, state label) at most ``per_state`` structures within ``window_kcal`` of the
group's lowest energy, plus every always-kept structure (discovery endpoints); each group keeps
its lowest structure (the Curtin-Hammett reference). ``rerank`` re-scores by single points only
the groups with more than ``per_state`` candidates (``crowded``): a smaller group is kept whole.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace

from hfauto.core.constants import HARTREE_TO_KCAL_MOL


@dataclass(frozen=True)
class Candidate:
    species_id: str
    composition_id: str
    state_label: str
    energy_hartree: float | None  # screen (or single-point) energy; None is unranked
    always: bool = False  # discovery endpoints are refined regardless of energy


def _state(c: Candidate) -> tuple[str, str]:
    return c.composition_id, c.state_label


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


def crowded(candidates: Sequence[Candidate], per_state: int) -> list[Candidate]:
    """The candidates, not always kept, of the groups that hold more than ``per_state`` of
    them: the only ones whose selection a single point can change."""
    count = Counter(_state(c) for c in candidates if not c.always)
    return [c for c in candidates if not c.always and count[_state(c)] > per_state]


def rerank(
    candidates: Sequence[Candidate],
    sp_energies: Mapping[str, float],
    keep_per_state: int,
    window_kcal: float,
) -> list[Candidate]:
    """Re-score the crowded groups with single-point energies and keep, per state, the
    ``keep_per_state`` lowest within ``window_kcal`` of the lowest (the better estimate cuts
    the window again); every other candidate stays.

    A crowded candidate without a single point drops out.
    """
    busy = crowded(candidates, keep_per_state)
    rescored = [replace(c, energy_hartree=sp_energies.get(c.species_id)) for c in busy]
    kept = select_for_refinement(rescored, per_state=keep_per_state, window_kcal=window_kcal)
    dropped = {c.species_id for c in busy} - {c.species_id for c in kept}
    return [c for c in candidates if c.species_id not in dropped]

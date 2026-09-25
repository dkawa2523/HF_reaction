"""Which structures the minima stage refines (design §8.2, minima(dft)).

``window``: per (composition, state label) at most ``per_state`` structures within
``window_kcal`` of the group's lowest energy, plus every always-kept structure (discovery
endpoints). Each group keeps its lowest structure, so the lowest structure of every input
monomer (needed for association thermochemistry) is always included. ``all``: everything.
"""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import Literal

from hfauto.core.constants import HARTREE_TO_KCAL_MOL


@dataclass(frozen=True)
class Candidate:
    species_id: str
    composition_id: str
    state_label: str
    energy_hartree: float | None  # screen (or single-point) energy; None is unranked
    always: bool = False  # discovery endpoints are refined regardless of energy


def _lowest_first(group: Sequence[Candidate]) -> list[tuple[float, Candidate]]:
    ranked = [(c.energy_hartree, c) for c in group if c.energy_hartree is not None]
    return sorted(ranked, key=lambda pair: (pair[0], pair[1].species_id))


def select_for_refinement(
    candidates: Sequence[Candidate],
    *,
    include: Literal["window", "all"] = "window",
    per_state: int = 3,
    window_kcal: float = 6.0,
) -> list[Candidate]:
    """The selected candidates in input order; unranked ones only when always kept."""
    if include == "all":
        return list(candidates)
    chosen = {c.species_id for c in candidates if c.always}
    states: defaultdict[tuple[str, str], list[Candidate]] = defaultdict(list)
    for c in candidates:
        states[(c.composition_id, c.state_label)].append(c)
    for group in states.values():
        ranked = _lowest_first(group)
        chosen.update(
            c.species_id for energy, c in ranked[:max(per_state, 1)]
            if (energy - ranked[0][0]) * HARTREE_TO_KCAL_MOL <= window_kcal
        )
    return [c for c in candidates if c.species_id in chosen]


def rerank(
    candidates: Sequence[Candidate], sp_energies: Mapping[str, float], keep_per_state: int
) -> list[Candidate]:
    """Re-score with single-point energies and keep the ``keep_per_state`` lowest per state.

    A candidate without a single point drops out unless it is always kept.
    """
    rescored = [replace(c, energy_hartree=sp_energies.get(c.species_id)) for c in candidates]
    return select_for_refinement(rescored, per_state=keep_per_state, window_kcal=math.inf)

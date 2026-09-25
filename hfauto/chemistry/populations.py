"""Boltzmann populations for a finite, explicitly defined minimum ensemble."""

from __future__ import annotations

from math import exp, isfinite
from typing import Any

from hfauto.core.constants import HARTREE_TO_KCAL_MOL

R_KCAL_MOL_K = 0.00198720425864083


def conditional_boltzmann_populations(
    entries: list[dict[str, Any]],
    *,
    temperature_K: float,
) -> list[dict[str, Any]]:
    """Return populations normalized over the supplied finite ensemble.

    Every entry must contain ``G_standard_hartree``. ``degeneracy`` defaults
    to one and represents a declared statistical degeneracy, not spin
    multiplicity. The output deliberately calls the result conditional: this
    function cannot establish that conformer or basin discovery was exhaustive.
    """

    temperature = float(temperature_K)
    if not isfinite(temperature) or temperature <= 0.0:
        raise ValueError("temperature_K must be finite and positive")
    if not entries:
        return []

    normalized: list[tuple[dict[str, Any], float, float]] = []
    for entry in entries:
        energy = float(entry["G_standard_hartree"])
        degeneracy = float(entry.get("degeneracy", 1.0))
        if not isfinite(energy):
            raise ValueError("G_standard_hartree must be finite")
        if not isfinite(degeneracy) or degeneracy <= 0.0:
            raise ValueError("degeneracy must be finite and positive")
        normalized.append((dict(entry), energy, degeneracy))

    minimum = min(energy for _entry, energy, _degeneracy in normalized)
    weighted: list[tuple[dict[str, Any], float, float]] = []
    for entry, energy, degeneracy in normalized:
        delta_kcal = (energy - minimum) * HARTREE_TO_KCAL_MOL
        weight = degeneracy * exp(-delta_kcal / (R_KCAL_MOL_K * temperature))
        weighted.append((entry, delta_kcal, weight))
    partition = sum(weight for _entry, _delta, weight in weighted)
    if not isfinite(partition) or partition <= 0.0:
        raise ValueError("Boltzmann partition function is not positive and finite")

    rows: list[dict[str, Any]] = []
    for entry, delta_kcal, weight in weighted:
        rows.append(
            {
                **entry,
                "T_K": temperature,
                "delta_G_ensemble_kcal_mol": delta_kcal,
                "boltzmann_weight": weight,
                "conditional_population": weight / partition,
                "population_is_conditional": True,
            }
        )
    return rows

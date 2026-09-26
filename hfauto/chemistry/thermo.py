"""Thermochemistry arithmetic (design §8.2 thermo): the frequencies GoodVibes gets, composite
G, standard states, ensembles, association and reaction deltas.

Pure functions. Energies are in Hartree unless a name ends in ``_kcal``; free energies of
single species are gas-phase 1 atm values (GoodVibes' default reference) and are moved to
other standard states only through ``standard_state_shift``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Literal

from hfauto.core.constants import HARTREE_TO_KCAL_MOL, R_KCAL_MOL_K
from hfauto.core.hashing import fingerprint_dict
from hfauto.core.method import ThermoSettings
from hfauto.core.records import CaseOutcome, ReactionRecord, SpeciesThermo

StandardState = Literal["1atm", "1bar", "1M"]

_L_ATM_PER_MOL_K = 0.082057366080960  # gas constant in L atm / (mol K), CODATA 2018
_ATM_PER_BAR = 1.0 / 1.01325
_TS_OUTCOMES = frozenset(
    {CaseOutcome.ELEMENTARY_STEP, CaseOutcome.DEGENERATE, CaseOutcome.REASSIGNED}
)


def thermo_frequencies(freqs_cm1: Sequence[float], *, saddle: bool) -> tuple[float, ...]:
    """The modes that enter ZPE, U and S, sorted: a minimum keeps every mode as |nu|; a saddle
    drops its lowest mode (the reaction coordinate) and keeps the others as |nu|."""
    modes = sorted(freqs_cm1)[1:] if saddle else freqs_cm1
    return tuple(sorted(abs(nu) for nu in modes))


def settings_sha(settings: ThermoSettings) -> str:
    """ThermoResult.settings_sha / SpeciesThermo.settings_sha of settings as passed."""
    return fingerprint_dict(settings.model_dump(mode="json"))


def composite(E_sp: float, G_gv: float, E_gv: float) -> float:
    """G = E_SP + (G_GV - E_GV): the thermal correction moved onto the energy layer."""
    return E_sp + (G_gv - E_gv)


def standard_state_shift(dn: float, T: float, to: StandardState) -> float:
    """kcal/mol added to a 1 atm reaction free energy with dn = n(products) - n(reactants)."""
    if to == "1atm":
        return 0.0
    pressure_atm = _ATM_PER_BAR if to == "1bar" else _L_ATM_PER_MOL_K * T  # 1 mol/L as a gas
    return dn * R_KCAL_MOL_K * T * math.log(pressure_atm)


def _rt_hartree(T: float) -> float:
    if not (math.isfinite(T) and T > 0.0):
        raise ValueError("temperature must be finite and positive")
    return R_KCAL_MOL_K * T / HARTREE_TO_KCAL_MOL


def boltzmann_populations(G_hartree: Sequence[float], T: float) -> list[float]:
    """Normalized Boltzmann weights of a finite ensemble (conditional on its completeness)."""
    if not G_hartree:
        return []
    rt, low = _rt_hartree(T), min(G_hartree)
    weights = [math.exp(-(g - low) / rt) for g in G_hartree]
    total = sum(weights)
    return [w / total for w in weights]


def ensemble_G(G_hartree: Sequence[float], T: float) -> float:
    """-RT ln sum exp(-G_i / RT): the free energy of the ensemble (<= its lowest member)."""
    if not G_hartree:
        raise ValueError("empty ensemble")
    rt, low = _rt_hartree(T), min(G_hartree)
    return low - rt * math.log(sum(math.exp(-(g - low) / rt) for g in G_hartree))


def association(
    G_complex: float | None,
    G_ts: float | None,
    monomer_G: Sequence[float],
    counts: Sequence[int],
    T: float,
    state: StandardState,
) -> tuple[float | None, float | None]:
    """(dG_assoc, dG_act_vs_separated) in kcal/mol against n_i separated monomers, with the
    standard-state correction for dn = 1 - sum(n_i)."""
    separated = sum(g * n for g, n in zip(monomer_G, counts, strict=True))
    shift = standard_state_shift(1 - sum(counts), T, state)

    def delta(G: float | None) -> float | None:
        return None if G is None else (G - separated) * HARTREE_TO_KCAL_MOL + shift

    return delta(G_complex), delta(G_ts)


def _kcal(a: float | None, b: float | None) -> float | None:
    return None if a is None or b is None else (a - b) * HARTREE_TO_KCAL_MOL


def participants(reaction: ReactionRecord) -> tuple[str, ...]:
    """Thermo subjects of a reaction: both minima, then the TS (SaddleClaim.freq_calc) for
    elementary / degenerate / reassigned outcomes with a saddle claim."""
    if reaction.saddle is None or reaction.outcome not in _TS_OUTCOMES:
        return reaction.minima
    return (*reaction.minima, reaction.saddle.freq_calc)


def reaction_delta(
    reaction: ReactionRecord, species_thermo: Mapping[str, SpeciesThermo]
) -> tuple[float | None, float | None, float | None]:
    """(dG_act, dG_rxn, dzpe_act) in kcal/mol at 1 atm from the thermo of each subject.

    dG_act and dzpe_act exist only when ``participants`` includes the TS; a degenerate
    rearrangement has dG_rxn = 0 by symmetry."""

    def value(subject: str, field: str) -> float | None:
        thermo = species_thermo.get(subject)
        return None if thermo is None else getattr(thermo, field)

    reactant, product, *ts = participants(reaction)
    dG_rxn = _kcal(value(product, "G_hartree"), value(reactant, "G_hartree"))
    if reaction.degenerate and dG_rxn is not None:
        dG_rxn = 0.0
    if not ts:
        return None, dG_rxn, None
    return (
        _kcal(value(ts[0], "G_hartree"), value(reactant, "G_hartree")),
        dG_rxn,
        _kcal(value(ts[0], "zpe_hartree"), value(reactant, "zpe_hartree")),
    )

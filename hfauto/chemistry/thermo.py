"""Thermochemistry (design §8.2 thermo): GoodVibes in this process, the modes it gets at the
point group's structure (chemistry.symmetry), the optical isomer term, standard states,
association, reaction deltas and the ranking quantity dG_eff.

Energies are in Hartree unless a name ends in ``_kcal``; free energies of single species are
gas-phase 1 atm values (GoodVibes' default reference) and are moved to other standard states
only through ``standard_state_shift``.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from typing import NamedTuple

import numpy as np

from hfauto.chemistry.elements import atomic_number, mass
from hfauto.chemistry.symmetry import Symmetry
from hfauto.chemistry.vibrations import projected_frequencies, rotational_constants_ghz
from hfauto.chemistry.xyz import hill_formula
from hfauto.core.constants import HARTREE_TO_KCAL_MOL, R_KCAL_MOL_K
from hfauto.core.hashing import fingerprint_dict
from hfauto.core.method import ThermoSettings
from hfauto.core.records import (
    CONNECTED_OUTCOMES,
    ReactionRecord,
    SpeciesRecord,
    StandardState,
)
from hfauto.core.system import CompositionInput

State = tuple[str, str]  # (composition_id, state_label)
Monomers = dict[tuple[str, int], list[tuple[State, int]]]  # complex (Hill, charge) -> parts

GOODVIBES_VERSION = "4.3.0"  # species_thermo is validated against it (tests/golden, G06)
_K_PER_GHZ = 0.0479924307  # h / k_B in K per GHz (CODATA 2018)
_L_ATM_PER_MOL_K = 0.082057366080960  # gas constant in L atm / (mol K), CODATA 2018
_ATM_PER_BAR = 1.0 / 1.01325


def thermo_frequencies(freqs_cm1: Sequence[float], *, saddle: bool) -> tuple[float, ...]:
    """The modes that enter ZPE, U and S, sorted: a minimum keeps every mode as |nu|; a saddle
    drops its lowest mode (the reaction coordinate) and keeps the others as |nu|."""
    modes = sorted(freqs_cm1)[1:] if saddle else freqs_cm1
    return tuple(sorted(abs(nu) for nu in modes))


def thermal_modes(symbols: Sequence[str], sym: Symmetry, hessian: np.ndarray, *, saddle: bool
                  ) -> tuple[float, ...]:
    """thermo_frequencies of the freq Hessian projected at the symmetrized structure with the
    point group's external count: 3N - 5 modes if linear, else 3N - 6 (an atom: none)."""
    freqs, _, _ = projected_frequencies(hessian, symbols, sym.coords, linear=sym.linear)
    return thermo_frequencies(freqs.tolist(), saddle=saddle)


class Thermal(NamedTuple):
    """Thermal terms of one species above its electronic energy (Hartree)."""

    G: float  # H - T S, S with the quasi-RRHO vibrational entropy of the settings
    H: float
    zpe: float


def species_thermo(symbols: Sequence[str], sym: Symmetry, modes_cm1: Sequence[float], *,
                   multiplicity: int, settings: ThermoSettings, T: float) -> Thermal:
    """GoodVibes (``compute_thermo``) in this process for an ideal gas at 1 atm: QH=False (H is
    RRHO), the vibrational entropy of ``settings.qs`` with ``settings.cutoff_cm1``, vib_scale
    on the modes and the ZPE, sigma and linearity of the point group and S_el = R ln(2S+1).

    The modes are ``thermal_modes``; the rotational temperatures come from the symmetrized
    structure. calc_bbe reads ``zero_point_corr`` only as a gate (None: nothing computed) and as
    the monatomic test (== 0.0), and computes nothing from an empty rotemp, so an atom gets 0.0
    and a dummy [0.0]. ``file`` stays "" (no such path), so calc_bbe never parses a file.
    """
    from goodvibes.api import compute_thermo
    from goodvibes.io import QCData

    symbols = list(symbols)
    atom = len(symbols) == 1
    constants = () if atom else rotational_constants_ghz(symbols, sym.coords, sym.linear)
    qcdata = QCData(
        scf_energy=0.0, multiplicity=multiplicity, atom_types=symbols,
        atom_nums=[atomic_number(s) for s in symbols], cartesians=sym.coords.tolist(),
        frequency_wn=list(modes_cm1), linear_mol=sym.linear,
        molecular_mass=sum(mass(s) for s in symbols),
        rotemp=[b * _K_PER_GHZ for b in constants if b > 0.0] or [0.0],
        zero_point_corr=0.0 if atom else 1.0, symmno=sym.sigma,
    )
    r = compute_thermo(qcdata=qcdata, QS=settings.qs, s_freq_cutoff=settings.cutoff_cm1,
                       temperature=T, freq_scale_factor=settings.vib_scale,
                       zpe_scale_factor=settings.vib_scale, symm=False)
    if r.qh_gibbs_free_energy is None or r.enthalpy is None or r.zpe is None:
        raise ValueError("GoodVibes computed no thermal terms")
    return Thermal(G=r.qh_gibbs_free_energy, H=r.enthalpy, zpe=r.zpe)


def monomer_states(species: Iterable[SpeciesRecord],
                   compositions: Iterable[CompositionInput]) -> Monomers:
    """(Hill formula, charge) of each composition of several molecules -> (state, count) of
    its parts: the separated reference."""
    by_id = {s.species_id: s for s in species}
    out: Monomers = {}
    for comp in compositions:
        parts = [(by_id[i], n) for i, n in comp.components.items() if i in by_id]
        if len(parts) == len(comp.components) and sum(comp.components.values()) > 1:
            symbols = [x for s, n in parts for x in s.geometry.symbols * n]
            key = (hill_formula(symbols), sum(s.charge * n for s, n in parts))
            out[key] = [((s.composition_id, s.state_label), n) for s, n in parts]
    return out


def settings_sha(settings: ThermoSettings) -> str:
    """SpeciesThermo.settings_sha of settings as passed."""
    return fingerprint_dict(settings.model_dump(mode="json"))


def standard_state_shift(dn: float, T: float, to: StandardState) -> float:
    """kcal/mol added to a 1 atm reaction free energy with dn = n(products) - n(reactants)."""
    if to == "1atm":
        return 0.0
    pressure_atm = _ATM_PER_BAR if to == "1bar" else _L_ATM_PER_MOL_K * T  # 1 mol/L as a gas
    return dn * R_KCAL_MOL_K * T * math.log(pressure_atm)


def chiral_G(T: float) -> float:
    """-RT ln 2 (Hartree): the G term of a chiral structure (Symmetry.m = 2: its point group has
    no improper operation), whose mirror image is the same basin (Fernandez-Ramos et al.,
    Theor. Chem. Acc. 118, 813)."""
    if not (math.isfinite(T) and T > 0.0):
        raise ValueError("temperature must be finite and positive")
    return -R_KCAL_MOL_K * T * math.log(2.0) / HARTREE_TO_KCAL_MOL


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


def participants(reaction: ReactionRecord) -> tuple[str, ...]:
    """Thermo subjects of a reaction: both minima, then the TS (SaddleClaim.freq_calc) for
    elementary / degenerate / reassigned outcomes with a saddle claim."""
    if reaction.saddle is None or reaction.outcome not in CONNECTED_OUTCOMES.values():
        return reaction.minima
    return (*reaction.minima, reaction.saddle.freq_calc)


def effective_barrier(G_ts: float | None, G_R: float, G_P: float) -> float:
    """dG_eff = max(G_TS, G_R, G_P) - G_R, in the units of the arguments: the highest free
    energy met between the reactant and product states. A TS below either state is no
    bottleneck (microscopic reversibility; Truhlar, Garrett, Klippenstein, J. Phys. Chem. 100,
    12771 (1996)); without a TS (barrierless, or a barrier submerged by ZPE) it is
    max(dG_rxn, 0)."""
    return max(G_R, G_P, G_R if G_ts is None else G_ts) - G_R

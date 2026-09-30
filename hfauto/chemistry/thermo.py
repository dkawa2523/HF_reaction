"""Thermochemistry (design §8.2 thermo): GoodVibes in this process, the frequencies it gets,
the optical isomer term, standard states, association, reaction deltas and the ranking
quantity dG_eff.

Energies are in Hartree unless a name ends in ``_kcal``; free energies of single species are
gas-phase 1 atm values (GoodVibes' default reference) and are moved to other standard states
only through ``standard_state_shift``.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from typing import NamedTuple

from hfauto.chemistry.elements import atomic_number, mass
from hfauto.chemistry.vibrations import rotational_constants_ghz
from hfauto.chemistry.xyz import XYZ, hill_formula
from hfauto.core.constants import HARTREE_TO_KCAL_MOL, R_KCAL_MOL_K
from hfauto.core.hashing import fingerprint_dict
from hfauto.core.method import ThermoSettings
from hfauto.core.records import (
    CONNECTED_OUTCOMES,
    ReactionRecord,
    SpeciesRecord,
    SpeciesThermo,
    StandardState,
)
from hfauto.core.system import CompositionInput

State = tuple[str, str]  # (composition_id, state_label)
Monomers = dict[tuple[str, int], list[tuple[State, int]]]  # complex (Hill, charge) -> parts

GOODVIBES_VERSION = "4.3.0"  # species_thermo is validated against it (tests/golden, G06)
_K_PER_GHZ = 0.0479924307  # h / k_B in K per GHz (CODATA 2018)
_L_ATM_PER_MOL_K = 0.082057366080960  # gas constant in L atm / (mol K), CODATA 2018
_ATM_PER_BAR = 1.0 / 1.01325
# libmsym's relative equivalence threshold. Its default 5e-4 split the C3v I-...CH3I complex into
# Cs when an optimisation left the I- 0.03 deg off the axis (sigma 3 -> 1, dG 0.65 kcal/mol run to
# run); 2e-3 gives the same point group as the default for every other validation structure.
_SYMMETRY_EQUIVALENCE = 2.0e-3


def symmetry_number(symbols: Sequence[str], coords: Sequence[Sequence[float]]) -> int:
    """External rotational symmetry number of the point group libmsym finds (pymsym, the
    detector GoodVibes uses) with ``_SYMMETRY_EQUIVALENCE``; 1 for an atom or when libmsym finds
    no consistent group (pymsym returns C1 then as well)."""
    import pymsym
    from pymsym.high_level import SYMMNO_BY_POINT_GROUP

    if len(symbols) == 1:
        return 1
    elements = [pymsym.Element(name=s, coordinates=[float(v) for v in c])
                for s, c in zip(symbols, coords, strict=True)]
    try:
        with pymsym.Context(elements=elements) as ctx:
            ctx.set_thresholds(equivalence=_SYMMETRY_EQUIVALENCE)
            group = ctx.find_symmetry()
    except Exception:  # libmsym raises when no group fits
        return 1
    return SYMMNO_BY_POINT_GROUP[{"D0h": "Dinfh", "C0v": "Cinfv"}.get(group, group)]


def thermo_frequencies(freqs_cm1: Sequence[float], *, saddle: bool) -> tuple[float, ...]:
    """The modes that enter ZPE, U and S, sorted: a minimum keeps every mode as |nu|; a saddle
    drops its lowest mode (the reaction coordinate) and keeps the others as |nu|."""
    modes = sorted(freqs_cm1)[1:] if saddle else freqs_cm1
    return tuple(sorted(abs(nu) for nu in modes))


class Thermal(NamedTuple):
    """Thermal terms of one species above its electronic energy (Hartree)."""

    G: float  # H - T S, S with the quasi-RRHO vibrational entropy of the settings
    H: float
    zpe: float


def species_thermo(xyz: XYZ, frequencies_cm1: Sequence[float], *, saddle: bool,
                   multiplicity: int, settings: ThermoSettings, T: float) -> Thermal:
    """GoodVibes (``compute_thermo``) in this process for an ideal gas at 1 atm: QH=False (H is
    RRHO), the vibrational entropy of ``settings.qs`` with ``settings.cutoff_cm1``, vib_scale
    on the modes and the ZPE, sigma from ``symmetry_number`` and S_el = R ln(2S+1).

    The modes are ``thermo_frequencies``; the mass and the rotational temperatures come from
    the geometry, and a vanishing moment of inertia marks a linear molecule. calc_bbe reads
    ``zero_point_corr`` only as a gate (None: nothing computed) and as the monatomic test
    (== 0.0), and computes nothing from an empty rotemp, so an atom gets 0.0 and a dummy [0.0].
    ``file`` stays "" (no such path), so calc_bbe never parses a file.
    """
    from goodvibes.api import compute_thermo
    from goodvibes.io import QCData

    symbols = list(xyz.symbols)
    constants = [b for b in rotational_constants_ghz(symbols, xyz.coords) if b > 0.0]
    qcdata = QCData(
        scf_energy=0.0, multiplicity=multiplicity, atom_types=symbols,
        atom_nums=[atomic_number(s) for s in symbols], cartesians=xyz.coords.tolist(),
        frequency_wn=list(thermo_frequencies(frequencies_cm1, saddle=saddle)),
        linear_mol=len(constants) == 2, molecular_mass=sum(mass(s) for s in symbols),
        rotemp=[b * _K_PER_GHZ for b in constants] or [0.0],
        zero_point_corr=0.0 if len(symbols) == 1 else 1.0,
        symmno=symmetry_number(symbols, xyz.coords.tolist()),
    )
    r = compute_thermo(qcdata=qcdata, QS=settings.qs, s_freq_cutoff=settings.cutoff_cm1,
                       temperature=T, freq_scale_factor=settings.vib_scale,
                       zpe_scale_factor=settings.vib_scale, symm=False)
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
    """-RT ln 2 (Hartree): the G term of a chiral structure, whose mirror image is the same
    basin (optical isomer number m = 2; Fernandez-Ramos et al., Theor. Chem. Acc. 118, 813)."""
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


def _kcal(a: float | None, b: float | None) -> float | None:
    return None if a is None or b is None else (a - b) * HARTREE_TO_KCAL_MOL


def participants(reaction: ReactionRecord) -> tuple[str, ...]:
    """Thermo subjects of a reaction: both minima, then the TS (SaddleClaim.freq_calc) for
    elementary / degenerate / reassigned outcomes with a saddle claim."""
    if reaction.saddle is None or reaction.outcome not in CONNECTED_OUTCOMES.values():
        return reaction.minima
    return (*reaction.minima, reaction.saddle.freq_calc)


def reaction_delta(
    reaction: ReactionRecord, species_thermo: Mapping[str, SpeciesThermo]
) -> tuple[float | None, float | None]:
    """(dG_act, dG_rxn) in kcal/mol at 1 atm from the thermo of each subject; dG_act exists only
    when ``participants`` includes the TS."""

    def value(subject: str) -> float | None:
        thermo = species_thermo.get(subject)
        return None if thermo is None else thermo.G_hartree

    reactant, product, *ts = participants(reaction)
    G_R = value(reactant)
    return (_kcal(value(ts[0]), G_R) if ts else None), _kcal(value(product), G_R)


def effective_barrier(G_ts: float | None, G_R: float, G_P: float) -> float:
    """dG_eff = max(G_TS, G_R, G_P) - G_R, in the units of the arguments: the highest free
    energy met between the reactant and product states. A TS below either state is no
    bottleneck (microscopic reversibility; Truhlar, Garrett, Klippenstein, J. Phys. Chem. 100,
    12771 (1996)); without a TS (barrierless, or a barrier submerged by ZPE) it is
    max(dG_rxn, 0)."""
    return max(G_R, G_P, G_R if G_ts is None else G_ts) - G_R

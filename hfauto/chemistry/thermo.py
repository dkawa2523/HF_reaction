"""Thermochemistry (design §8.2 thermo): GoodVibes in this process for the nuclear terms, the
electronic term from the levels of the ground term, the optical isomer term, standard states, and
the points a reaction reads with the ranking quantity over them.

The modes are those of the structure itself; the point group (chemistry.symmetry) gives sigma,
m, linearity and the moments of inertia of its symmetrised structure. GoodVibes sees a singlet:
the electronic term is ``electronic`` alone.

Energies are in Hartree unless a name ends in ``_kcal``; free energies of single species are
gas-phase 1 atm values (GoodVibes' default reference) and are moved to other standard states
only through ``standard_state_shift``.

The points a reaction reads (``reaction_points``) are one definition for the sp stage, the
thermo stage and the method panel. A reaction is a chain of points with a common zero,
[R_sep?, R, TS, P, P_sep?]: a separated point exists only where the end's fragments are the
declared monomers of a composition (``separated_states``); an association's reactant side is its
separated monomers (``ReactionRecord.monomers``), and a split child's ends lie inside its
parent's chain, so it has none.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from typing import NamedTuple

import numpy as np

from hfauto.chemistry import topology
from hfauto.chemistry.electronic_state import Levels
from hfauto.chemistry.elements import atomic_number, mass
from hfauto.chemistry.symmetry import Symmetry
from hfauto.chemistry.vibrations import projected_frequencies, rotational_constants_ghz
from hfauto.chemistry.xyz import XYZ, hill_formula
from hfauto.core.constants import CM1_TO_HARTREE, HARTREE_TO_KCAL_MOL, R_KCAL_MOL_K
from hfauto.core.evidence import Evidence
from hfauto.core.hashing import fingerprint_dict
from hfauto.core.manifest import Artifact
from hfauto.core.method import ThermoSettings
from hfauto.core.records import (
    CONNECTED_OUTCOMES,
    ReactionRecord,
    SpeciesRecord,
    StandardState,
)
from hfauto.core.system import CompositionInput

State = tuple[str, str]  # (composition_id, state_label)
Separated = tuple[State, ...]  # declared monomer states, one per molecule
Monomers = dict[tuple[str, int], list[Separated]]  # (Hill, charge) -> each composition's parts
Sides = tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]  # reactant, product, TS ids

GOODVIBES_VERSION = "4.3.0"  # species_thermo is validated against it (tests/golden, G06)
_K_PER_GHZ = 0.0479924307  # h / k_B in K per GHz (CODATA 2018)
_L_ATM_PER_MOL_K = 0.082057366080960  # gas constant in L atm / (mol K), CODATA 2018
_ATM_PER_BAR = 1.0 / 1.01325


def thermo_frequencies(freqs_cm1: Sequence[float], *, saddle: bool) -> tuple[float, ...]:
    """The modes that enter ZPE, U and S, sorted: a minimum keeps every mode as |nu|; a saddle
    drops its lowest mode (the reaction coordinate) and keeps the others as |nu|."""
    modes = sorted(freqs_cm1)[1:] if saddle else freqs_cm1
    return tuple(sorted(abs(nu) for nu in modes))


def thermal_modes(symbols: Sequence[str], coords: np.ndarray, hessian: np.ndarray, *,
                  linear: bool, saddle: bool) -> tuple[float, ...]:
    """thermo_frequencies of the freq Hessian projected at the structure itself with the point
    group's external count: 3N - 5 modes if ``linear``, else 3N - 6 (an atom: none)."""
    freqs, _, _ = projected_frequencies(hessian, symbols, coords, linear=linear)
    return thermo_frequencies(freqs.tolist(), saddle=saddle)


class Thermal(NamedTuple):
    """Nuclear thermal terms of one species above its electronic energy (Hartree)."""

    G: float  # H - T S, S with the quasi-RRHO vibrational entropy of the settings
    H: float  # with the quasi-harmonic vibrational energy at the same cutoff
    zpe: float


def species_thermo(symbols: Sequence[str], sym: Symmetry, modes_cm1: Sequence[float], *,
                   settings: ThermoSettings, T: float) -> Thermal:
    """GoodVibes (``compute_thermo``) in this process for an ideal gas at 1 atm and multiplicity
    1: the vibrational entropy of ``settings.qs`` with ``settings.cutoff_cm1``, QH=True (Head-
    Gordon's damped vibrational energy) with that cutoff, vib_scale on the modes and the ZPE,
    sigma, linearity and the moments of inertia of the point group. With qs = grimme the
    vibrational H and S damp towards the free rotor at one cutoff; truhlar's S stays RRQHO.

    The modes are ``thermal_modes``. calc_bbe reads ``zero_point_corr`` only as a gate (None:
    nothing computed) and as the monatomic test (== 0.0), and computes nothing from an empty
    rotemp, so an atom gets 0.0 and a dummy [0.0]. ``file`` stays "" (no such path), so
    calc_bbe never parses a file.
    """
    from goodvibes.api import compute_thermo
    from goodvibes.io import QCData

    symbols = list(symbols)
    atom = len(symbols) == 1
    constants = () if atom else rotational_constants_ghz(symbols, sym.coords, sym.linear)
    qcdata = QCData(
        scf_energy=0.0, multiplicity=1, atom_types=symbols,
        atom_nums=[atomic_number(s) for s in symbols], cartesians=sym.coords.tolist(),
        frequency_wn=list(modes_cm1), linear_mol=sym.linear,
        molecular_mass=sum(mass(s) for s in symbols),
        rotemp=[b * _K_PER_GHZ for b in constants if b > 0.0] or [0.0],
        zero_point_corr=0.0 if atom else 1.0, symmno=sym.sigma,
    )
    r = compute_thermo(qcdata=qcdata, QS=settings.qs, QH=True, s_freq_cutoff=settings.cutoff_cm1,
                       h_freq_cutoff=settings.cutoff_cm1, temperature=T,
                       freq_scale_factor=settings.vib_scale, zpe_scale_factor=settings.vib_scale,
                       symm=False)
    if r.qh_gibbs_free_energy is None or r.qh_enthalpy is None or r.zpe is None:
        raise ValueError("GoodVibes computed no thermal terms")
    return Thermal(G=r.qh_gibbs_free_energy, H=r.qh_enthalpy, zpe=r.zpe)


def electronic(levels: Levels, T: float) -> float:
    """G_el = -RT ln q_el (Hartree) of the levels (degeneracy, cm⁻¹) of a ground term, q_el =
    Σ g exp(-(E - E_0)/kT) above the lowest level E_0; a spin multiplet alone, ((2S + 1, 0),),
    gives -RT ln(2S + 1)."""
    if not (math.isfinite(T) and T > 0.0):
        raise ValueError("temperature must be finite and positive")
    kT = R_KCAL_MOL_K * T / HARTREE_TO_KCAL_MOL
    low = min(E for _, E in levels)
    return -kT * math.log(sum(g * math.exp(-(E - low) * CM1_TO_HARTREE / kT) for g, E in levels))


def spin_orbit(levels: Levels) -> float:
    """E_SO = E_0 - Σ g E / Σ g (Hartree): the lowest level against the term's mean, which a
    non-relativistic energy describes."""
    mean = sum(g * E for g, E in levels) / sum(g for g, _ in levels)
    return (min(E for _, E in levels) - mean) * CM1_TO_HARTREE


def declared_monomers(species: Iterable[SpeciesRecord],
                      compositions: Iterable[CompositionInput]) -> Monomers:
    """(Hill formula, charge) -> the parts of every declared composition of several molecules
    with that formula and charge: its monomers' states, one per molecule."""
    by_id = {s.species_id: s for s in species}
    out: Monomers = {}
    for comp in compositions:
        parts = [by_id[i] for i, n in comp.components.items() if i in by_id for _ in range(n)]
        if len(parts) == sum(comp.components.values()) > 1:
            key = (hill_formula([x for s in parts for x in s.geometry.symbols]),
                   sum(s.charge for s in parts))
            out.setdefault(key, []).append(tuple((s.composition_id, s.state_label) for s in parts))
    return out


def separated_states(xyz: XYZ, charge: int, monomers: Monomers) -> Separated:
    """The declared monomers whose state labels are the fragments of ``xyz`` (as a multiset,
    their charges summing to ``charge``), one per molecule; () when no composition's are."""
    labels = sorted(topology.state_label([xyz.symbols[i] for i in group], xyz.coords[list(group)])
                    for group in topology.fragments(xyz.symbols, xyz.coords))
    return next((parts for parts in monomers.get((hill_formula(xyz.symbols), charge), ())
                 if sorted(label for _, label in parts) == labels), ())


def single_points(calculations: Iterable[Artifact]
                  ) -> dict[tuple[str, str], tuple[str, Evidence]]:
    """(subject, Level.full_key) -> (calculation id, Evidence) of every sp calculation, its
    subjects being its parents. Two sp calculations of one subject at one Level are an error:
    no order decides between them."""
    out: dict[tuple[str, str], tuple[str, Evidence]] = {}
    for a in calculations:
        if isinstance(ev := a.payload, Evidence) and ev.task == "sp":
            for subject in a.parents:
                if (key := (subject, ev.level.full_key())) in out:
                    raise ValueError(f"two sp calculations of {subject} at one level: "
                                     f"{out[key][0]}, {a.artifact_id}")
                out[key] = (a.artifact_id, ev)
    return out


def settings_sha(settings: ThermoSettings) -> str:
    """SpeciesThermo.settings_sha of settings as passed."""
    return fingerprint_dict(settings.model_dump(mode="json"))


def standard_state_shift(n: float, T: float, to: StandardState) -> float:
    """kcal/mol added to the 1 atm free energy of n molecules (to a reaction free energy with
    n = n(products) - n(reactants))."""
    if to == "1atm":
        return 0.0
    pressure_atm = _ATM_PER_BAR if to == "1bar" else _L_ATM_PER_MOL_K * T  # 1 mol/L as a gas
    return n * R_KCAL_MOL_K * T * math.log(pressure_atm)


def chiral_G(T: float) -> float:
    """-RT ln 2 (Hartree): the G term of a chiral structure (Symmetry.m = 2: its point group has
    no improper operation), whose mirror image is the same basin (Fernandez-Ramos et al.,
    Theor. Chem. Acc. 118, 813)."""
    if not (math.isfinite(T) and T > 0.0):
        raise ValueError("temperature must be finite and positive")
    return -R_KCAL_MOL_K * T * math.log(2.0) / HARTREE_TO_KCAL_MOL


def participants(reaction: ReactionRecord) -> Sides:
    """A reaction's own points: its reactant (an association: its separated monomers, each by
    count), product and, for an elementary / degenerate / reassigned outcome with a saddle
    claim, the TS (SaddleClaim.freq_calc). dE, dG_act, dG_rxn and the method panel read them."""
    ts = (reaction.saddle.freq_calc,) if (
        reaction.saddle is not None and reaction.outcome in CONNECTED_OUTCOMES.values()) else ()
    return reaction.monomers or reaction.minima[:1], reaction.minima[1:], ts


class Point(NamedTuple):
    """A well (``states``, one per molecule, each read on the LOT of the reaction end
    ``subject``) or the saddle (no states; ``subject`` is the TS). Its standard-state shift is
    that of its molecules: len(states), a saddle one."""

    states: Separated
    subject: str


class ReactionPoints(NamedTuple):
    """The points a reaction reads. ``own``: participants. ``chain``: dG_eff of a connected
    outcome, [R_sep?, R, TS, P, P_sep?] (an association's collapsed complex is no R).
    ``separated`` (R_sep) and ``complex`` (R where R_sep exists): dG_assoc and
    dG_act_vs_separated; for a barrierless association they are auxiliary points."""

    own: Sides
    chain: tuple[Point, ...]
    separated: Point | None
    complex: Point | None

    def reads(self, state_of: Mapping[str, State]) -> tuple[list[State | str], list[State]]:
        """(main, aux): the states and TS subjects the own points and the chain read, and the
        states only the auxiliary points (dG_assoc) read."""
        reactant, product, ts = self.own
        main = [*(state_of[i] for i in (*reactant, *product) if i in state_of), *ts,
                *(k for pt in self.chain for k in pt.states or (pt.subject,))]
        aux = [s for pt in (self.separated, self.complex) if pt for s in pt.states]
        return main, [s for s in aux if s not in main]


def reaction_points(reaction: ReactionRecord, state_of: Mapping[str, State],
                    separated: Mapping[str, Separated]) -> ReactionPoints:
    """``state_of``: dft minimum id -> its state; ``separated``: reaction end -> its
    separated_states. A reaction whose ends are no dft minima reads its own points only."""
    own = participants(reaction)
    r, p = reaction.minima
    if any(m not in state_of for m in (r, p, *reaction.monomers)):
        return ReactionPoints(own, (), None, None)
    split = reaction.source == "split"
    r_sep = tuple(state_of[m] for m in reaction.monomers) or (
        () if split else separated.get(r, ()))
    p_sep = () if split else separated.get(p, ())
    sep = Point(r_sep, r) if r_sep else None
    rp = None if reaction.monomers and r == p else Point((state_of[r],), r)
    found = (sep, rp, *(Point((), t) for t in own[2]), Point((state_of[p],), p),
             Point(p_sep, p) if p_sep else None)
    chain = tuple(pt for pt in found if pt is not None) if own[2] else ()
    return ReactionPoints(own, chain, sep, rp if sep else None)


def chain_barrier(points: Sequence[tuple[float, bool]]) -> float:
    """dG_eff = max_k (G_k - min over wells j <= k of G_j) over a chain of (G, is_well): the
    highest rise above the lowest well met before it. A TS below its reactant is no
    bottleneck (microscopic reversibility; Truhlar, Garrett, Klippenstein, J. Phys. Chem. 100,
    12771 (1996)), and for [R, TS, P] this is max(G_TS, G_R, G_P) - G_R; it is the serial-path
    form of the energetic span (Kozuch, Shaik, Acc. Chem. Res. 44, 101 (2011))."""
    low, rise = math.inf, 0.0
    for G, well in points:
        low = min(low, G) if well else low
        rise = max(rise, G - low)
    return rise

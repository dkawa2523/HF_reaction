"""thermo stage (design §4.1 #7, §8.2 thermo): GoodVibes in this process
(chemistry.thermo.species_thermo), composite G, sensitivity band and the reaction values over
the points of thermo.reaction_points.

Subjects: every dft minimum (monomers included) and SaddleClaim freq calc, at the main
settings and the qs x cutoff variants and at every configured temperature. Each subject's
point group (chemistry.symmetry, from the freq geometry and Hessian; open shell: multiplicity
above 1) gives sigma, linearity, the moments of inertia and m: with m = 2 G gains -RT ln 2, as
the mirror image is the same basin or saddle, counted once. The modes are those of the freq
structure (chemistry.thermo.thermal_modes: a TS drops its reaction coordinate). Species values
are 1 atm; reactions get one record per standard state. E is the energy layer: with
energy_method the one sp at that method whose parents name the subject (thermo.single_points;
without one G = None, energy_layer_missing), else the freq itself. A subject is
spin_contaminated unless its energy passes gates.energy_spin_ok (the method panel's definition
too).

The electronic term (thermo.electronic, thermo.spin_orbit) comes from the levels of the ground
term: a minimum of a declared species' state with electronic_levels (a linear radical with
Λ > 0) takes them, an isolated atom its NIST term (electronic_state.atom_levels), anything else
the spin multiplet alone. E_SO joins E (dE, G and the submerged test); in a complex or a TS
the spin-orbit splitting counts as quenched. G = E + nuclear terms + G_el (+ -RT ln 2).

A state is a composition and a topology.state_label. Its G (_state_G) is the lowest G of its
spin-clean minima on one freq and energy LOT (fast conformer equilibria make one state); it fails
closed (None) when one of them lacks its G, e.g. its energy layer. A well of the chain sums its
states' G on the LOT of its reaction end and every point gets the standard-state shift of its
molecule count. dG_eff (thermo.chain_barrier) is defined for a connected outcome only, its
zero the lower of R and R_sep (``reference``); the TS drops out, with the note
submerged_barrier, when its forward or reverse dE0 = dE + dZPE on the own points is <= 0. A
barrierless outcome has no dG_eff: its row gives dG_rxn. dG_act, dG_rxn and dE refer to the own
points (thermo.participants), dG_assoc and dG_act_vs_separated to R_sep. The blockers of a
ReactionThermo are the only source of thermo_unavailable (the row's value, dG_eff or else
dG_rxn, is missing), mixed_level_of_theory (the own points differ in LOT, charge and
multiplicity aside), spin_contaminated and not_stationary (an own point is).
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import ClassVar, cast

import numpy as np

from hfauto.chemistry import symmetry
from hfauto.chemistry import thermo as th
from hfauto.chemistry.electronic_state import Levels, atom_levels
from hfauto.chemistry.gates import energy_spin_ok, same_pes
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Evidence, Level
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec, ThermoSettings, level_mismatches
from hfauto.core.records import (
    ArtifactType,
    MinimumRecord,
    ReactionRecord,
    ReactionThermo,
    SpeciesRecord,
    SpeciesThermo,
    StandardState,
)
from hfauto.stages.spec import StageConfig, StageRuntime, StageSpec

_QS, _CUTOFFS_CM1 = ("grimme", "truhlar"), (50.0, 100.0, 150.0)
Table = dict[str, SpeciesThermo]  # subject -> thermo at one temperature and one settings


class ThermoConfig(StageConfig):
    settings: ThermoSettings = ThermoSettings()
    energy_method: str | None = None  # sp layer of the composite G
    temperatures_K: tuple[float, ...] = (298.15,)
    standard_states: tuple[StandardState, ...] = ("1atm",)


@dataclass(frozen=True)
class _Subject:
    id: str  # minimum_id or SaddleClaim.freq_calc
    freq_calc: str
    freq: Evidence
    energy_calc: str | None  # the sp at energy_method whose parents name this subject
    energy: Evidence  # that sp, else the freq itself
    state: th.State | None  # None for a TS
    notes: tuple[str, ...]  # of the MinimumRecord / SaddleClaim
    separated: th.Separated  # a minimum's thermo.separated_states
    contaminated: bool  # the freq fails spin_ok, or the energy layer is on another spin state
    layer_missing: bool  # energy_method set but no such sp
    sym: symmetry.Symmetry  # of the freq geometry
    modes: tuple[float, ...]  # th.thermal_modes
    levels: Levels  # of the electronic term

    @property
    def E(self) -> float:
        """The energy layer with the spin-orbit lowering of the ground level."""
        return self.energy.energy_hartree + th.spin_orbit(self.levels)


def _variants(settings: ThermoSettings) -> list[ThermoSettings]:
    grid = [settings] + [settings.model_copy(update={"qs": qs, "cutoff_cm1": c})
                         for qs in _QS for c in _CUTOFFS_CM1 if settings.sensitivity]
    return list({th.settings_sha(s): s for s in grid}.values())  # main first, no repeats


def _subjects(inputs: Manifest, method: MethodSpec | None, rt: StageRuntime
              ) -> dict[str, _Subject]:
    layer = {} if method is None else {  # subject -> (calc id, Evidence) of its sp at method
        subject: sp for (subject, _), sp in
        th.single_points(inputs.of(ArtifactType.CALCULATION)).items()
        if not level_mismatches(method, sp[1].level, version_pin=sp[1].level.version)}
    species = inputs.records(ArtifactType.SPECIES, SpeciesRecord)
    monomers = th.declared_monomers(species, rt.system.compositions)
    declared = {s.id: s.electronic_levels for s in rt.system.species if s.electronic_levels}
    levels_of = {(s.composition_id, s.state_label): declared[s.species_id]
                 for s in species if s.source == "input" and s.species_id in declared}

    def make(sid: str, calc: str, state: th.State | None, notes: tuple[str, ...]) -> _Subject:
        freq = inputs.evidence(calc)
        assert freq.hessian is not None  # Evidence guarantee (3) of a freq
        xyz, hessian = rt.load_xyz(freq.final), np.load(rt.resolve(freq.hessian))
        charge, mult = freq.level.charge, freq.level.multiplicity
        sym = symmetry.analyze(xyz.symbols, xyz.coords, hessian, open_shell=mult > 1)
        sp_calc, sp = layer.get(sid, (None, freq))
        levels = None if state is None else levels_of.get(state) or (
            atom_levels(xyz.symbols[0], charge, mult) if len(xyz.symbols) == 1 else None)
        return _Subject(sid, calc, freq, sp_calc, sp, state, notes,
                        () if state is None else th.separated_states(xyz, charge, monomers),
                        contaminated=not energy_spin_ok(sp, freq, rt.policy),
                        layer_missing=method is not None and sp_calc is None, sym=sym,
                        modes=th.thermal_modes(xyz.symbols, xyz.coords, hessian,
                                               linear=sym.linear, saddle=state is None),
                        levels=levels or ((mult, 0.0),))

    out = {m.minimum_id: make(m.minimum_id, m.freq_calc, (m.composition_id, m.state_label),
                              m.notes)
           for m in inputs.records(ArtifactType.MINIMUM, MinimumRecord) if m.tier == "dft"}
    for r in inputs.records(ArtifactType.REACTION, ReactionRecord):
        if r.saddle is not None:
            calc = r.saddle.freq_calc
            out[calc] = make(calc, calc, None, r.saddle.notes)
    return out


def _species(sub: _Subject, settings: ThermoSettings, T: float) -> SpeciesThermo:
    """G = E + nuclear thermal terms of the freq + G_el (- RT ln 2 if m = 2)."""
    row = SpeciesThermo(
        subject=sub.id, freq_calc=sub.freq_calc, energy_calc=sub.energy_calc, T_K=T,
        G_hartree=None, H_hartree=None, zpe_hartree=None, settings_sha=th.settings_sha(settings),
        point_group=sub.sym.point_group, sigma=sub.sym.sigma, m=sub.sym.m, notes=sub.notes)
    if sub.layer_missing:  # G stays None
        return row.model_copy(update={"notes": (*sub.notes, "energy_layer_missing")})
    c = th.species_thermo(sub.freq.final.symbols, sub.sym, sub.modes, settings=settings, T=T)
    mirror = th.chiral_G(T) if sub.sym.m == 2 else 0.0
    return row.model_copy(update={
        "G_hartree": sub.E + c.G + th.electronic(sub.levels, T) + mirror,
        "H_hartree": sub.E + c.H, "zpe_hartree": c.zpe})


def _members(state: th.State, ref: _Subject, subjects: dict[str, _Subject]) -> list[_Subject]:
    """The spin-clean minima of ``state`` on the LOT of ``ref`` (charge and multiplicity aside):
    its freq LOT and, for a minimum with its energy layer, its energy LOT."""
    return [s for s in subjects.values()
            if s.state == state and not s.contaminated
            and same_pes(ref.freq.level, s.freq.level, state=False)
            and (s.layer_missing or same_pes(ref.energy.level, s.energy.level, state=False))]


def _state_G(state: th.State, ref: _Subject, subjects: dict[str, _Subject], table: Table
             ) -> float | None:
    """The G of ``state``: the lowest G of its members; None without one or when a member
    lacks its G (energy_layer_missing): a missing conformer may be the lowest."""
    G = [table[s.id].G_hartree for s in _members(state, ref, subjects)]
    return None if not G or None in G else min(cast(list[float], G))


def _sum(values: Iterable[float | None]) -> float | None:
    """The sum of one side's values; None for an empty side or a missing value."""
    found = list(values)
    return None if not found or None in found else sum(cast(list[float], found))


def _kcal(G: float | None, n: int, T: float, state: StandardState) -> float | None:
    """G (Hartree) of n molecules in kcal/mol at the standard state."""
    return None if G is None else G * HARTREE_TO_KCAL_MOL + th.standard_state_shift(n, T, state)


def _point(pt: th.Point | None, subjects: dict[str, _Subject], table: Table, T: float,
           state: StandardState) -> float | None:
    """G of a chain point in kcal/mol: a well sums its states' G on its end's LOT."""
    if pt is None or pt.subject not in subjects:
        return None
    end = subjects[pt.subject]
    G = _sum(_state_G(s, end, subjects, table) for s in pt.states) if pt.states else (
        table[end.id].G_hartree)
    return _kcal(G, len(pt.states) or 1, T, state)


def _barrier(chain: Sequence[th.Point], subjects: dict[str, _Subject], table: Table, T: float,
             state: StandardState) -> float | None:
    G = [_point(pt, subjects, table, T, state) for pt in chain]
    return None if None in G else th.chain_barrier(
        [(cast(float, g), bool(pt.states)) for g, pt in zip(G, chain, strict=True)])


def _label(level: Level) -> str:  # e.g. pbe0-d3bj/def2-svpd
    method = f"{level.method}-{level.dispersion}" if level.dispersion else level.method
    return f"{method}/{level.basis}" if level.basis else method


def _energy_level(subs: Sequence[_Subject]) -> str | None:
    """The energy layer all own points share; None when it is missing or mixed."""
    labels = {None if s.layer_missing else _label(s.energy.level) for s in subs}
    return labels.pop() if len(labels) == 1 else None


def _diff(a: float | None, b: float | None) -> float | None:
    return None if a is None or b is None else a - b


def _E(subjects: dict[str, _Subject], ids: Sequence[str]) -> float | None:
    """The energy layer summed over ``ids`` in kcal/mol; None without one (no mixed-LOT dE)."""
    E = _sum(None if (s := subjects.get(i)) is None or s.layer_missing else s.E for i in ids)
    return None if E is None else E * HARTREE_TO_KCAL_MOL


def _submerged(sides: th.Sides, subjects: dict[str, _Subject], table: Table) -> bool:
    """Forward or reverse dE0 = dE + dZPE <= 0 on the energy layer: the saddle is no
    bottleneck and its TST barrier has no meaning."""
    def e0(ids: Sequence[str]) -> float | None:
        return _sum(subjects[n].E + st.zpe_hartree
                    if (st := table.get(n)) is not None and st.zpe_hartree is not None else None
                    for n in ids)

    r, p, ts = (e0(ids) for ids in sides)
    return r is not None and p is not None and ts is not None and ts <= max(r, p)


def _blockers(subs: Sequence[_Subject], value: float | None) -> tuple[str, ...]:
    """thermo_unavailable: the row's value is missing; mixed_level_of_theory: the own points
    differ in LOT (charge and multiplicity aside); spin_contaminated or not_stationary: one of
    them is (the driver's note when the trust-region certification failed after its relax)."""
    hits = {"thermo_unavailable": value is None,
            "mixed_level_of_theory": not (same_pes(*[s.freq.level for s in subs], state=False)
                                          and same_pes(*[s.energy.level for s in subs],
                                                       state=False)),
            "spin_contaminated": any(s.contaminated for s in subs),
            "not_stationary": any("not_stationary" in s.notes for s in subs)}
    return tuple(name for name, hit in hits.items() if hit)


def _reaction(rx: ReactionRecord, points: th.ReactionPoints, subjects: dict[str, _Subject],
              T: float, state: StandardState, tables: Sequence[Table]) -> ReactionThermo:
    """tables: the thermo of every settings variant at T (index 0 = main settings)."""
    table, (reactant, product, ts) = tables[0], points.own
    subs = [subjects[p] for p in dict.fromkeys((*reactant, *product, *ts)) if p in subjects]

    def own(ids: Sequence[str]) -> float | None:  # own minima's G with the shift of their count
        return _kcal(_sum(table[i].G_hartree if i in table else None for i in ids), len(ids),
                     T, state)

    submerged = _submerged(points.own, subjects, table)
    chain = [pt for pt in points.chain if pt.states or not submerged]
    band = [_barrier(chain, subjects, t, T, state) for t in tables] if chain else [None]
    G_R, G_sep, G_cpx = own(reactant), *(_point(p, subjects, table, T, state)
                                          for p in (points.separated, points.complex))
    dG_rxn = _diff(own(product), G_R)
    reference = None if band[0] is None or G_sep is None else (
        "complex" if G_cpx is not None and G_cpx < G_sep else "separated")
    return ReactionThermo(
        reaction_id=rx.reaction_id, T_K=T, standard_state=state,
        dE_act_kcal=_diff(_E(subjects, ts), _E(subjects, reactant)),
        dE_rxn_kcal=_diff(_E(subjects, product), _E(subjects, reactant)),
        dG_act_kcal=_diff(own(ts), G_R), dG_rxn_kcal=dG_rxn,
        dG_assoc_kcal=_diff(G_cpx, G_sep), dG_act_vs_separated_kcal=_diff(own(ts), G_sep),
        band_kcal=None if None in band else (min(cast(list[float], band)),
                                             max(cast(list[float], band))),
        blockers=_blockers(subs, band[0] if points.chain else dG_rxn), dG_eff_kcal=band[0],
        reference=reference, energy_level=_energy_level(subs),
        notes=("submerged_barrier",) if submerged else (),
    )


class ThermoStage:
    spec: ClassVar[StageSpec] = StageSpec(
        "thermo", ThermoConfig, consumes=(ArtifactType.MINIMUM, ArtifactType.CALCULATION),
        produces=(ArtifactType.SPECIES_THERMO, ArtifactType.REACTION_THERMO))

    def run(self, inputs: Manifest, config: StageConfig, rt: StageRuntime) -> list[Artifact]:
        cfg = cast(ThermoConfig, config)
        method = None if cfg.energy_method is None else rt.method(cfg.energy_method)
        subjects = _subjects(inputs, method, rt)
        variants = _variants(cfg.settings)  # index 0 = main settings
        state_of = {s.id: s.state for s in subjects.values() if s.state is not None}
        separated = {s.id: s.separated for s in subjects.values()}
        reactions = [(rx, th.reaction_points(rx, state_of, separated))
                     for rx in inputs.records(ArtifactType.REACTION, ReactionRecord)]
        out: list[Artifact] = []
        for T in cfg.temperatures_K:
            tables = [{sid: _species(sub, s, T) for sid, sub in subjects.items()}
                      for s in variants]
            out += [Artifact(
                artifact_id=f"species_thermo_{st.subject}_{T:g}K",
                type=ArtifactType.SPECIES_THERMO, payload=st,
                parents=tuple(dict.fromkeys(
                    filter(None, (st.subject, st.freq_calc, st.energy_calc)))),
            ) for st in tables[0].values()]
            out += [Artifact(
                artifact_id=f"reaction_thermo_{rx.reaction_id}_{T:g}K_{state}",
                type=ArtifactType.REACTION_THERMO, parents=(rx.reaction_id,),
                payload=_reaction(rx, points, subjects, T, state, tables),
            ) for rx, points in reactions for state in cfg.standard_states]
        return out

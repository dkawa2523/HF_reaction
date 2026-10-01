"""thermo stage (design §4.1 #7, §8.2 thermo): GoodVibes in this process
(chemistry.thermo.species_thermo), composite G, sensitivity band, association
thermochemistry and the ranking quantity dG_eff.

Subjects: every dft minimum (monomers included) and SaddleClaim freq calc, at the main
settings and the qs x cutoff variants and at every configured temperature. Each subject's
point group (chemistry.symmetry, from the freq geometry and Hessian) gives sigma, linearity,
the structure the modes (chemistry.thermo.thermal_modes: a TS drops its reaction coordinate)
and moments are evaluated at, and m: with m = 2 G gains -RT ln 2, as the mirror image is the
same basin or saddle, counted once. Species values are 1 atm; reactions get one record per
standard state. G is the energy layer plus the thermal terms of the freq: with energy_method
the sp whose parents name the subject (the later one in the view wins; without it G = None,
energy_layer_missing), else the freq itself. The energy layer passes spin_ok like the freq.
The blockers of a ReactionThermo are the only source of thermo_unavailable,
mixed_level_of_theory and spin_contaminated.

A state is a composition and a topology.state_label. Its G (_state_G) is the lowest G of its
spin-clean minima on one freq and energy LOT (fast conformer equilibria make one state). dG_eff
(thermo.effective_barrier) refers to the reactant's and the product's state, dG_assoc to the
monomers' states on the complex's LOT; the TS drops out, with the note submerged_barrier, when
its forward or reverse dE0 = dE + dZPE is <= 0.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import ClassVar, cast

import numpy as np

from hfauto.chemistry import symmetry
from hfauto.chemistry import thermo as th
from hfauto.chemistry.gates import same_pes, spin_ok
from hfauto.chemistry.xyz import hill_formula
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Evidence, Level
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec, ThermoSettings, level_mismatches
from hfauto.core.records import (
    ArtifactType,
    CaseOutcome,
    MinimumRecord,
    ReactionRecord,
    ReactionThermo,
    SpeciesRecord,
    SpeciesThermo,
    StandardState,
)
from hfauto.stages.spec import StageConfig, StageRuntime, StageSpec

_QS, _CUTOFFS_CM1 = ("grimme", "truhlar"), (50.0, 100.0, 150.0)
State = th.State
Monomers = th.Monomers
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
    state: State | None  # None for a TS
    notes: tuple[str, ...]  # of the MinimumRecord / SaddleClaim
    formula: tuple[str, int]  # (Hill formula, charge): finds a complex's monomers
    contaminated: bool  # spin_contaminated in notes, or the energy layer fails spin_ok
    layer_missing: bool  # energy_method set but no such sp
    sym: symmetry.Symmetry  # of the freq geometry
    modes: tuple[float, ...]  # th.thermal_modes


def _variants(settings: ThermoSettings) -> list[ThermoSettings]:
    grid = [settings] + [settings.model_copy(update={"qs": qs, "cutoff_cm1": c})
                         for qs in _QS for c in _CUTOFFS_CM1 if settings.sensitivity]
    return list({th.settings_sha(s): s for s in grid}.values())  # main first, no repeats


def _subjects(inputs: Manifest, method: MethodSpec | None, rt: StageRuntime
              ) -> dict[str, _Subject]:
    layer = {} if method is None else {  # subject -> the sp at method naming it
        parent: (a.artifact_id, a.payload)
        for a in inputs.of(ArtifactType.CALCULATION)
        if isinstance(a.payload, Evidence) and a.payload.task == "sp"
        and not level_mismatches(method, a.payload.level, version_pin=a.payload.level.version)
        for parent in a.parents
    }

    def make(sid: str, calc: str, state: State | None, notes: tuple[str, ...]) -> _Subject:
        freq = inputs.evidence(calc)
        if freq.hessian is None:
            raise ValueError(f"{calc}: a freq without a Hessian")
        xyz, hessian = rt.load_xyz(freq.final), np.load(rt.resolve(freq.hessian))
        sym = symmetry.analyze(xyz.symbols, xyz.coords, hessian)
        sp_calc, sp = layer.get(sid, (None, freq))
        return _Subject(sid, calc, freq, sp_calc, sp, state, notes,
                        (hill_formula(freq.final.symbols), freq.level.charge),
                        contaminated="spin_contaminated" in notes or not spin_ok(sp, rt.policy),
                        layer_missing=method is not None and sp_calc is None, sym=sym,
                        modes=th.thermal_modes(xyz.symbols, sym, hessian, saddle=state is None))

    out = {m.minimum_id: make(m.minimum_id, m.freq_calc, (m.composition_id, m.state_label),
                              m.notes)
           for m in inputs.records(ArtifactType.MINIMUM, MinimumRecord) if m.tier == "dft"}
    for r in inputs.records(ArtifactType.REACTION, ReactionRecord):
        if r.saddle is not None:
            calc = r.saddle.freq_calc
            out[calc] = make(calc, calc, None, r.saddle.notes)
    return out


def _species(sub: _Subject, settings: ThermoSettings, T: float) -> SpeciesThermo:
    """G = energy layer + thermal terms of the freq (- RT ln 2 if m = 2)."""
    row = SpeciesThermo(
        subject=sub.id, freq_calc=sub.freq_calc, energy_calc=sub.energy_calc, T_K=T,
        G_hartree=None, H_hartree=None, zpe_hartree=None, settings_sha=th.settings_sha(settings),
        notes=sub.notes)
    if sub.layer_missing:
        return row.model_copy(update={
            "notes": (*sub.notes, "thermo_unavailable", "energy_layer_missing")})
    c = th.species_thermo(sub.freq.final.symbols, sub.sym, sub.modes,
                          multiplicity=sub.freq.level.multiplicity, settings=settings, T=T)
    E, mirror = sub.energy.energy_hartree, th.chiral_G(T) if sub.sym.m == 2 else 0.0
    return row.model_copy(update={
        "G_hartree": E + c.G + mirror, "H_hartree": E + c.H, "zpe_hartree": c.zpe})


def _same_level(subs: Sequence[_Subject], *, state: bool = True) -> bool:
    return bool(same_pes(*[s.freq.level for s in subs], state=state)) and bool(
        same_pes(*[s.energy.level for s in subs], state=state))


def _state_G(state: State | None, ref: _Subject, subjects: dict[str, _Subject], table: Table
             ) -> float | None:
    """The G of ``state``: the lowest G of its spin-clean minima on the freq and energy LOT of
    ``ref`` (charge and multiplicity aside); None without one."""
    G = [g for s in subjects.values()
         if s.state == state and not s.contaminated and _same_level([ref, s], state=False)
         and (g := table[s.id].G_hartree) is not None]
    return min(G, default=None)


def _association(names: tuple[str, ...], subjects: dict[str, _Subject], monomers: Monomers,
                 table: Table, T: float, state: StandardState, with_ts: bool
                 ) -> tuple[float | None, float | None]:
    """dG_assoc and dG_act_vs_separated of the complex against the state G of each monomer on
    its LOT; a monomer state without one fails closed."""
    complex_ = subjects.get(names[0])
    parts = monomers.get(complex_.formula) if complex_ is not None else None
    if complex_ is None or not parts:
        return None, None
    G = [_state_G(part, complex_, subjects, table) for part, _ in parts]
    if None in G:
        return None, None
    G_ts = table[names[2]].G_hartree if with_ts else None
    return th.association(table[names[0]].G_hartree, G_ts, cast(list[float], G),
                          [n for _, n in parts], T, state)


def _blockers(names: tuple[str, ...], subs: Sequence[_Subject], table: Table
              ) -> tuple[str, ...]:
    hits = {
        "thermo_unavailable": any(p not in table or table[p].G_hartree is None for p in names),
        "mixed_level_of_theory": not _same_level(subs),
        "spin_contaminated": any(s.contaminated for s in subs),
    }
    return tuple(name for name, hit in hits.items() if hit)


def _label(level: Level) -> str:  # e.g. pbe0-d3bj/def2-svpd
    method = f"{level.method}-{level.dispersion}" if level.dispersion else level.method
    return f"{method}/{level.basis}" if level.basis else method


def _energy_level(subs: Sequence[_Subject]) -> str | None:
    """The energy layer all participants share; None when it is missing or mixed."""
    labels = {None if s.layer_missing else _label(s.energy.level) for s in subs}
    return labels.pop() if len(labels) == 1 else None


def _kcal(subjects: dict[str, _Subject], a: str, b: str) -> float | None:
    sa, sb = subjects.get(a), subjects.get(b)
    if sa is None or sb is None or sa.layer_missing or sb.layer_missing:  # no mixed-LOT dE
        return None
    return (sa.energy.energy_hartree - sb.energy.energy_hartree) * HARTREE_TO_KCAL_MOL


def _submerged(names: tuple[str, ...], subjects: dict[str, _Subject], table: Table) -> bool:
    """Forward or reverse dE0 = dE + dZPE <= 0 on the energy layer: the saddle is no
    bottleneck and its TST barrier has no meaning."""
    e0 = [subjects[n].energy.energy_hartree + st.zpe_hartree for n in names
          if (st := table.get(n)) is not None and st.zpe_hartree is not None]
    return len(names) == len(e0) == 3 and e0[2] <= max(e0[:2])


def _effective(rx: ReactionRecord, names: tuple[str, ...], subjects: dict[str, _Subject],
               table: Table, shifts: tuple[float, float], submerged: bool) -> float | None:
    """dG_eff in kcal/mol from one settings variant: defined with a TS participant or for a
    barrierless outcome; a submerged TS drops out."""
    with_ts = len(names) == 3 and not submerged
    if len(names) < 3 and rx.outcome is not CaseOutcome.BARRIERLESS:
        return None
    ends = (subjects.get(name) for name in names[:2])
    G_R, G_P = (None if s is None else _state_G(s.state, s, subjects, table) for s in ends)
    G_ts = table[names[2]].G_hartree if with_ts else None
    if G_R is None or G_P is None or (with_ts and G_ts is None):
        return None
    ts = None if G_ts is None else (G_ts - G_R) * HARTREE_TO_KCAL_MOL + shifts[0]
    return th.effective_barrier(ts, 0.0, (G_P - G_R) * HARTREE_TO_KCAL_MOL + shifts[1])


def _reaction(rx: ReactionRecord, subjects: dict[str, _Subject], monomers: Monomers, T: float,
              state: StandardState, tables: Sequence[Table]) -> ReactionThermo:
    """tables: the thermo of every settings variant at T (index 0 = main settings)."""
    names, table = th.participants(rx), tables[0]
    subs = [subjects[p] for p in names if p in subjects]
    dG_act, dG_rxn = th.reaction_delta(rx, table)
    n_r, n_p = (sum(t.coefficient for t in terms) for terms in (rx.reactants, rx.products))
    shifts = (th.standard_state_shift(1 - n_r, T, state),
              th.standard_state_shift(n_p - n_r, T, state))
    submerged = _submerged(names, subjects, table)
    band = [_effective(rx, names, subjects, t, shifts, submerged) for t in tables]
    assoc, vs_separated = _association(names, subjects, monomers, table, T, state,
                                       dG_act is not None)
    return ReactionThermo(
        reaction_id=rx.reaction_id, T_K=T, standard_state=state,
        dE_act_kcal=_kcal(subjects, names[2], names[0]) if len(names) == 3 else None,
        dE_rxn_kcal=_kcal(subjects, names[1], names[0]),
        dG_act_kcal=None if dG_act is None else dG_act + shifts[0],
        dG_rxn_kcal=None if dG_rxn is None else dG_rxn + shifts[1],
        dG_assoc_kcal=assoc, dG_act_vs_separated_kcal=vs_separated,
        band_kcal=None if None in band else (min(cast(list[float], band)),
                                             max(cast(list[float], band))),
        blockers=_blockers(names, subs, table), dG_eff_kcal=band[0],
        energy_level=_energy_level(subs), notes=("submerged_barrier",) if submerged else (),
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
        monomers = th.monomer_states(inputs.records(ArtifactType.SPECIES, SpeciesRecord),
                                     rt.system.compositions)
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
                payload=_reaction(rx, subjects, monomers, T, state, tables),
            ) for rx in inputs.records(ArtifactType.REACTION, ReactionRecord)
                for state in cfg.standard_states]
        return out

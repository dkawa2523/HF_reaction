"""thermo stage (design §4.1 #7, §8.2 thermo): GoodVibes, consistency gate, composite G,
sensitivity band, association thermochemistry and basin populations.

Subjects: every dft minimum (monomers included) and SaddleClaim freq calc; one engine call
each covers the main settings and the qs x cutoff variants at the conditions' temperatures.
Modes follow chemistry.thermo.thermo_frequencies (a TS drops its reaction coordinate). A
failed call or gate gives G = None (thermo_unavailable), never retried. Species values are
1 atm; reactions get one record per standard state. With energy_method, G = E_SP + (G_GV -
E_GV) from the sp on the same geometry; a subject without that sp is energy_layer_missing.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from typing import ClassVar, cast

from hfauto.backends.protocols import Capability, ThermoEngine
from hfauto.chemistry import thermo as th
from hfauto.chemistry.gates import Gate, Policy, same_pes, thermo_consistent
from hfauto.chemistry.xyz import hill_formula
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Evidence, Failure
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec, ThermoSettings, level_mismatches
from hfauto.core.records import (
    ArtifactType,
    MinimumRecord,
    ReactionRecord,
    ReactionThermo,
    SpeciesRecord,
    SpeciesThermo,
)
from hfauto.stages.spec import StageConfig, StageRuntime, StageSpec

_QS, _CUTOFFS_CM1 = ("grimme", "truhlar"), (50.0, 100.0, 150.0)
Key = tuple[str, int, int]  # (Hill formula, charge, multiplicity) of a structure
Monomers = dict[tuple[str, int], list[tuple[Key, int]]]  # complex (Hill, charge) -> parts
Table = dict[str, SpeciesThermo]  # subject -> thermo at one temperature and one settings


class ThermoConfig(StageConfig):
    engine: str  # thermo engine name; no default
    settings: ThermoSettings = ThermoSettings()
    energy_method: str | None = None  # sp layer of the composite G


@dataclass(frozen=True)
class _Subject:
    id: str  # minimum_id or SaddleClaim.freq_calc
    freq_calc: str
    freq: Evidence
    energy_calc: str | None  # sp at energy_method on the same geometry
    energy: Evidence  # that sp, else the freq itself
    composition: str | None  # None for a TS
    notes: tuple[str, ...]  # of the MinimumRecord / SaddleClaim
    key: Key
    layer_missing: bool  # energy_method set but no sp on this geometry


def _variants(settings: ThermoSettings) -> list[ThermoSettings]:
    grid = [settings] + [settings.model_copy(update={"qs": qs, "cutoff_cm1": c})
                         for qs in _QS for c in _CUTOFFS_CM1 if settings.sensitivity]
    return list({th.settings_sha(s): s for s in grid}.values())  # main first, no repeats


def _subjects(inputs: Manifest, method: MethodSpec | None) -> dict[str, _Subject]:
    layer = {} if method is None else {
        a.payload.start.fingerprint: (a.artifact_id, a.payload)
        for a in inputs.of(ArtifactType.CALCULATION)
        if isinstance(a.payload, Evidence) and a.payload.task == "sp"
        and not level_mismatches(method, a.payload.level, version_pin=a.payload.level.version)
    }

    def make(sid: str, calc: str, composition: str | None, notes: tuple[str, ...]) -> _Subject:
        freq = inputs.evidence(calc)
        sp_calc, sp = layer.get(freq.final.fingerprint, (None, freq))
        key = (hill_formula(freq.final.symbols), freq.level.charge, freq.level.multiplicity)
        return _Subject(sid, calc, freq, sp_calc, sp, composition, notes, key,
                        layer_missing=method is not None and sp_calc is None)

    out = {m.minimum_id: make(m.minimum_id, m.freq_calc, m.composition_id, m.notes)
           for m in inputs.records(ArtifactType.MINIMUM, MinimumRecord) if m.tier == "dft"}
    for r in inputs.records(ArtifactType.REACTION, ReactionRecord):
        if r.saddle is not None:
            out[r.saddle.freq_calc] = make(r.saddle.freq_calc, r.saddle.freq_calc, None,
                                           r.saddle.notes)
    return out


def _species(sub: _Subject, engine: ThermoEngine, variants: Sequence[ThermoSettings],
             temperatures: Sequence[float], policy: Policy) -> dict[float, list[SpeciesThermo]]:
    """Per temperature, the thermo of every variant (index 0 = main settings)."""
    shas, saddle = [th.settings_sha(s) for s in variants], sub.composition is None  # TS
    out = [] if sub.layer_missing else engine.thermo(
        sub.freq, variants, temperatures_K=temperatures, saddle=saddle)
    found = {} if isinstance(out, Failure) else {(r.settings_sha, r.T_K): r for r in out}
    missing = (str(out.kind),) if isinstance(out, Failure) else (
        "energy_layer_missing" if sub.layer_missing else "no_thermo_result",)
    modes = th.thermo_frequencies(sub.freq.frequencies_cm1 or (), saddle=saddle)
    rows: dict[float, list[SpeciesThermo]] = {}
    for T in temperatures:
        results = [found.get((sha, T)) for sha in shas]
        r0 = results[0]
        gate = Gate(False, missing) if r0 is None else thermo_consistent(
            sub.freq, frequencies_cm1=modes, gv_zpe_hartree=r0.zpe_hartree,
            gv_energy_hartree=r0.E_hartree, gv_n_real=r0.n_real, scale=variants[0].vib_scale,
            policy=policy)
        base = SpeciesThermo(
            subject=sub.id, freq_calc=sub.freq_calc, energy_calc=sub.energy_calc, T_K=T,
            G_hartree=None, H_hartree=None, zpe_hartree=None, settings_sha=shas[0],
            notes=(*sub.notes, "thermo_unavailable", *gate.reasons))
        rows[T] = [base if not gate or r is None else base.model_copy(update={
            "G_hartree": th.composite(sub.energy.energy_hartree, r.G_hartree, r.E_hartree),
            "H_hartree": th.composite(sub.energy.energy_hartree, r.H_hartree, r.E_hartree),
            "zpe_hartree": r.zpe_hartree,
            "notes": tuple(dict.fromkeys((*sub.notes, *r.notes))),
        }) for r in results]
    return rows


def _populated(table: Table, subjects: dict[str, _Subject], T: float) -> Table:
    """Boltzmann weights of the minima within each composition and level of theory."""
    groups: dict[tuple, list[SpeciesThermo]] = defaultdict(list)
    for sid, st in table.items():
        s = subjects[sid]
        if s.composition is not None and st.G_hartree is not None:
            groups[s.composition, s.freq.level.full_key(), s.energy.level.full_key()].append(st)
    return table | {
        st.subject: st.model_copy(update={"population": w}) for sts in groups.values()
        for st, w in zip(sts, th.boltzmann_populations([cast(float, x.G_hartree) for x in sts], T),
                         strict=True)}


def _monomers(inputs: Manifest, rt: StageRuntime) -> Monomers:
    """(Hill, charge) of each system composition -> (monomer key, count) of its components."""
    species = {s.species_id: s for s in inputs.records(ArtifactType.SPECIES, SpeciesRecord)}
    out: Monomers = {}
    for comp in rt.system.compositions:
        parts = [(species[i], n) for i, n in comp.components.items() if i in species]
        if len(parts) == len(comp.components) and sum(comp.components.values()) > 1:
            symbols = [x for s, n in parts for x in s.geometry.symbols * n]
            key = (hill_formula(symbols), sum(s.charge * n for s, n in parts))
            out[key] = [((hill_formula(s.geometry.symbols), s.charge, s.multiplicity), n)
                        for s, n in parts]
    return out


def _same_level(subs: Sequence[_Subject], *, state: bool = True) -> bool:
    return bool(same_pes(*[s.freq.level for s in subs], state=state)) and bool(
        same_pes(*[s.energy.level for s in subs], state=state))


def _association(names: tuple[str, ...], subjects: dict[str, _Subject], monomers: Monomers,
                 table: Table, T: float, state: th.StandardState, with_ts: bool
                 ) -> tuple[float | None, float | None]:
    """dG_assoc and dG_act_vs_separated against the monomers' ensemble G (same LOT apart from
    charge and multiplicity; a monomer without G fails closed)."""
    complex_ = subjects.get(names[0])
    parts = monomers.get(complex_.key[:2]) if complex_ is not None else None
    if complex_ is None or not parts:
        return None, None
    energies, used = [], [complex_]
    for key, _ in parts:
        subs = [s for s in subjects.values() if s.composition is not None and s.key == key]
        G = [table[s.id].G_hartree for s in subs]
        if not G or None in G:
            return None, None
        energies.append(th.ensemble_G(cast(list[float], G), T))
        used += subs
    if not _same_level(used, state=False):
        return None, None
    G_ts = table[names[2]].G_hartree if with_ts else None
    return th.association(table[names[0]].G_hartree, G_ts, energies, [n for _, n in parts],
                          T, state)


def _blockers(names: tuple[str, ...], subjects: dict[str, _Subject], table: Table
              ) -> tuple[str, ...]:
    subs = [subjects[p] for p in names if p in subjects]
    hits = {
        "thermo_unavailable": any(p not in table or table[p].G_hartree is None for p in names),
        "mixed_level_of_theory": not _same_level(subs),
        "spin_contaminated": any("spin_contaminated" in s.notes for s in subs),
    }
    return tuple(name for name, hit in hits.items() if hit)


def _kcal(subjects: dict[str, _Subject], a: str, b: str) -> float | None:
    sa, sb = subjects.get(a), subjects.get(b)
    return None if sa is None or sb is None else (
        sa.energy.energy_hartree - sb.energy.energy_hartree) * HARTREE_TO_KCAL_MOL


def _reaction(rx: ReactionRecord, subjects: dict[str, _Subject], monomers: Monomers, T: float,
              state: th.StandardState, tables: Sequence[Table]) -> ReactionThermo:
    """tables: the thermo of every settings variant at T (index 0 = main settings)."""
    names, table = th.participants(rx), tables[0]
    dG_act, dG_rxn, dzpe = th.reaction_delta(rx, table)
    n_r, n_p = (sum(t.coefficient for t in terms) for terms in (rx.reactants, rx.products))
    shift_act, shift_rxn = (th.standard_state_shift(dn, T, state) for dn in (1 - n_r, n_p - n_r))
    band = [th.reaction_delta(rx, t)[0] for t in tables] if dG_act is not None else [None]
    dE_rxn = _kcal(subjects, names[1], names[0])
    if rx.degenerate and dE_rxn is not None:
        dE_rxn = 0.0
    assoc, vs_separated = _association(names, subjects, monomers, table, T, state,
                                       dG_act is not None)
    return ReactionThermo(
        reaction_id=rx.reaction_id, T_K=T, standard_state=state,
        dE_act_kcal=_kcal(subjects, names[2], names[0]) if len(names) == 3 else None,
        dE_rxn_kcal=dE_rxn, dzpe_act_kcal=dzpe,
        dG_act_kcal=None if dG_act is None else dG_act + shift_act,
        dG_rxn_kcal=None if dG_rxn is None else dG_rxn + shift_rxn,
        dG_assoc_kcal=assoc, dG_act_vs_separated_kcal=vs_separated,
        band_kcal=None if None in band else (
            min(cast(list[float], band)) + shift_act, max(cast(list[float], band)) + shift_act),
        blockers=_blockers(names, subjects, table),
    )


class ThermoStage:
    spec: ClassVar[StageSpec] = StageSpec(
        "thermo", ThermoConfig, consumes=(ArtifactType.MINIMUM, ArtifactType.CALCULATION),
        produces=(ArtifactType.SPECIES_THERMO, ArtifactType.REACTION_THERMO))

    def run(self, inputs: Manifest, config: StageConfig, rt: StageRuntime) -> list[Artifact]:
        cfg = cast(ThermoConfig, config)
        engine = cast(ThermoEngine, rt.engine(Capability.THERMO, cfg.engine))
        method = None if cfg.energy_method is None else rt.method(cfg.energy_method)
        subjects = _subjects(inputs, method)
        temperatures, variants = rt.conditions.temperatures_K, _variants(cfg.settings)
        rows = rt.thread_map(lambda s: _species(s, engine, variants, temperatures, rt.policy),
                             list(subjects.values()), threads_per_item=1)
        monomers = _monomers(inputs, rt)
        out: list[Artifact] = []
        for T in temperatures:
            tables = [{s.id: r[T][i] for s, r in zip(subjects.values(), rows, strict=True)}
                      for i in range(len(variants))]
            tables[0] = _populated(tables[0], subjects, T)
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
                for state in rt.conditions.standard_states]
        return out

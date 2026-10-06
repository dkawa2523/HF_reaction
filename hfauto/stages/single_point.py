"""sp stage (design §4.1 #6, §8.2 sp): fixed-geometry single points on the points the reactions
read.

Energy layers for composite free energies and method panels (WFT methods included). The targets
are the points of thermo.reaction_points of every reaction with a connected or barrierless
outcome: every dft minimum of a state a point reads (fast conformer equilibria: the state's G
is its lowest minimum) and the TS. Each target is computed on the final geometry of its freq
calculation, starting from that freq's SCF (scf_guess: design §7.1); the calculation's parents
are its subjects (a minimum_id, or SaddleClaim.freq_calc for a TS; several when subjects share
one freq: one job). The jobs are independent: they run at once (StageRuntime.thread_map),
each contained on its own. A failed single point of an auxiliary point only (R_sep and the
precursor complex of a barrierless step: dG_assoc) is left out, so it blocks no value and fails
no stage. CCSD(T) runs only where the freq passes spin_ok: a broken-symmetry point has no
single-reference state to correlate.
"""

from __future__ import annotations

from typing import ClassVar, cast

from hfauto.backends.protocols import Capability, QMEngine
from hfauto.chemistry import thermo as th
from hfauto.chemistry.gates import RANKABLE_OUTCOMES, spin_ok
from hfauto.chemistry.xyz import Molecule
from hfauto.core.evidence import Evidence, Failure
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.records import (
    ArtifactType,
    CaseOutcome,
    MinimumRecord,
    ReactionRecord,
    SpeciesRecord,
)
from hfauto.drivers.minimum import calc_id
from hfauto.stages.spec import StageConfig, StageRuntime, StageSpec

_VALUED = RANKABLE_OUTCOMES | {CaseOutcome.BARRIERLESS}  # outcomes the report gives values of


class SinglePointConfig(StageConfig):
    engine: str  # QM engine name; no default
    methods: list[str]


def _subjects(inputs: Manifest, rt: StageRuntime) -> dict[str, tuple[Evidence, bool]]:
    """Subject -> (its freq Evidence, False for an auxiliary point only)."""
    minima = {m.minimum_id: m for m in inputs.records(ArtifactType.MINIMUM, MinimumRecord)
              if m.tier == "dft"}
    freq = {i: inputs.evidence(m.freq_calc) for i, m in minima.items()}
    state_of = {i: (m.composition_id, m.state_label) for i, m in minima.items()}
    reactions = [r for r in inputs.records(ArtifactType.REACTION, ReactionRecord)
                 if r.outcome in _VALUED]
    monomers = th.declared_monomers(inputs.records(ArtifactType.SPECIES, SpeciesRecord),
                                    rt.system.compositions)
    separated = {m: th.separated_states(rt.load_xyz(freq[m].final), freq[m].level.charge,
                                        monomers)
                 for m in {m for r in reactions for m in r.minima if m in freq}}
    read: dict[th.State | str, bool] = {}  # a state or a TS subject -> read by a main value
    for rx in reactions:
        main, aux = th.reaction_points(rx, state_of, separated).reads(state_of)
        read |= {s: read.get(s, False) for s in aux} | dict.fromkeys(main, True)
    return {i: (freq[i], read[s]) for i, s in state_of.items() if s in read} | {
        t: (inputs.evidence(t), True) for t in read if isinstance(t, str)}


def _artifact(subjects: list[str], method_id: str, result: Evidence | Failure) -> Artifact:
    if isinstance(result, Failure):
        return Artifact(artifact_id=f"sp_{subjects[0]}_{method_id}",
                        type=ArtifactType.CALCULATION, parents=tuple(subjects), status="failed",
                        failure=result)
    return Artifact(artifact_id=calc_id(result), type=ArtifactType.CALCULATION,
                    parents=tuple(subjects), payload=result)


class SinglePointStage:
    spec: ClassVar[StageSpec] = StageSpec(
        name="sp",
        config=SinglePointConfig,
        consumes=(ArtifactType.MINIMUM, ArtifactType.REACTION),
        produces=(ArtifactType.CALCULATION,),
    )

    def run(self, inputs: Manifest, config: StageConfig, rt: StageRuntime) -> list[Artifact]:
        cfg = cast(SinglePointConfig, config)
        engine = cast(QMEngine, rt.engine(Capability.QM, cfg.engine))
        subjects = _subjects(inputs, rt)
        jobs: dict[tuple[str, str], list[str]] = {}  # (method, freq) -> the subjects sharing it
        for method_id in cfg.methods:
            for subject, (freq, _) in subjects.items():
                if rt.method(method_id).kind != "wft" or spin_ok(freq, rt.policy):
                    jobs.setdefault((method_id, freq.job_key), []).append(subject)

        def energy(job: tuple[tuple[str, str], list[str]]) -> Evidence | Failure:
            (method_id, _), names = job
            freq = subjects[names[0]][0]

            def sp() -> Evidence | Failure:
                mol = Molecule(rt.load_xyz(freq.final), freq.level.charge,
                               freq.level.multiplicity)
                return engine.energy(mol, rt.method(method_id), scf_guess=freq)

            return rt.contain(names[0], sp)

        todo = list(jobs.items())
        return [_artifact(names, method_id, result)
                for ((method_id, _), names), result in zip(todo, rt.thread_map(energy, todo),
                                                           strict=True)
                if not isinstance(result, Failure) or any(subjects[n][1] for n in names)]

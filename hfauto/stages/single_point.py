"""sp stage (design §4.1 #6, §8.2 sp): fixed-geometry single points on the points the reactions
read.

Energy layers for composite free energies and method panels (WFT methods included). The targets
are the points of thermo.reaction_points of every reaction with a connected or barrierless
outcome: every dft minimum of a state a point reads (fast conformer equilibria: the state's G
is its lowest minimum) and the TS. Each target is computed on the final geometry of its freq
calculation, starting from that freq's SCF (scf_guess: design §7.1); the calculation's parents
are its subjects (a minimum_id, or SaddleClaim.freq_calc for a TS; several when subjects share
one freq). A failed single point of an auxiliary point only (R_sep and the precursor complex of
a barrierless step: dG_assoc) is left out, so it blocks no value and fails no stage. CCSD(T)
runs only where the freq passes spin_ok: a broken-symmetry point has no single-reference state
to correlate.
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


def _artifact(subject: str, method_id: str, result: Evidence | Failure) -> Artifact:
    if isinstance(result, Failure):
        return Artifact(artifact_id=f"sp_{subject}_{method_id}", type=ArtifactType.CALCULATION,
                        parents=(subject,), status="failed", failure=result)
    return Artifact(artifact_id=calc_id(result), type=ArtifactType.CALCULATION,
                    parents=(subject,), payload=result)


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
        out: dict[str, Artifact] = {}
        for method_id in cfg.methods:  # DFT and WFT jobs run serially
            method = rt.method(method_id)
            for subject, (freq, main) in subjects.items():
                if method.kind == "wft" and not spin_ok(freq, rt.policy):
                    continue
                mol = Molecule(rt.load_xyz(freq.final), freq.level.charge,
                               freq.level.multiplicity)
                result = engine.energy(mol, method, scf_guess=freq)
                if isinstance(result, Failure) and not main:
                    continue  # an auxiliary value stays empty
                artifact = _artifact(subject, method_id, result)
                if artifact.artifact_id in out:  # one freq shared by several subjects
                    parents = out[artifact.artifact_id].parents + artifact.parents
                    artifact = artifact.model_copy(update={"parents": parents})
                out[artifact.artifact_id] = artifact
        return list(out.values())

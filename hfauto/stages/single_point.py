"""sp stage (design §4.1 #6, §8.2 sp): fixed-geometry single points on the points the ranking
reads.

Energy layers for composite free energies and method panels (WFT methods included). One
rule picks the targets: for every rankable reaction (gates' rankable outcomes), its TS and
every dft minimum of its reactant and product states (the lowest-G reference of dG_eff) and
of the monomer states of its composition (the separated reference). Each target is computed
on the final geometry of its freq calculation; the calculation's parents are its subjects (a
minimum_id, or SaddleClaim.freq_calc for a TS; several when subjects share one geometry).
"""

from __future__ import annotations

from typing import ClassVar, cast

from hfauto.backends.protocols import Capability, QMEngine
from hfauto.chemistry.gates import RANKABLE_OUTCOMES
from hfauto.chemistry.thermo import State, monomer_states, participants
from hfauto.chemistry.xyz import Molecule, hill_formula
from hfauto.core.evidence import Evidence, Failure
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.records import ArtifactType, MinimumRecord, ReactionRecord, SpeciesRecord
from hfauto.core.system import SystemConfig
from hfauto.drivers.minimum import calc_id
from hfauto.stages.spec import StageConfig, StageRuntime, StageSpec


class SinglePointConfig(StageConfig):
    engine: str  # QM engine name; no default
    methods: list[str]


def _targets(inputs: Manifest, system: SystemConfig) -> dict[str, Evidence]:
    """Subject id -> freq Evidence whose final geometry is computed (see the module doc)."""
    minima = {m.minimum_id: m for m in inputs.records(ArtifactType.MINIMUM, MinimumRecord)
              if m.tier == "dft"}
    monomers = monomer_states(inputs.records(ArtifactType.SPECIES, SpeciesRecord),
                              system.compositions)
    states: dict[State, None] = {}  # ordered set
    saddles: dict[str, Evidence] = {}
    for rx in inputs.records(ArtifactType.REACTION, ReactionRecord):
        ends = [minima[i] for i in rx.minima if i in minima]
        if rx.outcome not in RANKABLE_OUTCOMES or not ends:
            continue
        freq = inputs.evidence(ends[0].freq_calc)
        formula = (hill_formula(freq.final.symbols), freq.level.charge)
        states |= dict.fromkeys([*((m.composition_id, m.state_label) for m in ends),
                                 *(state for state, _ in monomers.get(formula, []))])
        if len(names := participants(rx)) == 3:
            saddles[names[2]] = inputs.evidence(names[2])
    return {mid: inputs.evidence(m.freq_calc) for mid, m in minima.items()
            if (m.composition_id, m.state_label) in states} | saddles


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
        subjects = _targets(inputs, rt.system)
        out: dict[str, Artifact] = {}
        for method_id in cfg.methods:  # DFT and WFT jobs run serially
            method = rt.method(method_id)
            for subject, freq in subjects.items():
                mol = Molecule(rt.load_xyz(freq.final), freq.level.charge,
                               freq.level.multiplicity)
                artifact = _artifact(subject, method_id, engine.energy(mol, method))
                if artifact.artifact_id in out:  # one geometry shared by several subjects
                    parents = out[artifact.artifact_id].parents + artifact.parents
                    artifact = artifact.model_copy(update={"parents": parents})
                out[artifact.artifact_id] = artifact
        return list(out.values())

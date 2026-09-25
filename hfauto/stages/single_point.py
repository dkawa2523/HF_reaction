"""sp stage (design §4.1 #6, §8.2 sp): fixed-geometry single points on stationary points.

Energy layers for composite free energies and method panels (WFT methods included). Each
target is computed on the final geometry of its freq calculation, so the thermo stage can
match a single point to its subject by geometry fingerprint; the calculation's parent is
the subject (a minimum_id, or SaddleClaim.freq_calc for a TS).
"""

from __future__ import annotations

from typing import ClassVar, Literal, cast

from hfauto.backends.protocols import Capability, QMEngine
from hfauto.chemistry.xyz import Molecule
from hfauto.core.evidence import Evidence, Failure
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.records import ArtifactType, MinimumRecord, ReactionRecord
from hfauto.drivers.minimum import calc_id
from hfauto.stages.spec import StageConfig, StageRuntime, StageSpec


class SinglePointConfig(StageConfig):
    engine: str  # QM engine name; no default
    methods: list[str]
    # all_minima: every dft minimum (monomers included, for association thermochemistry)
    # in addition to the reaction stationary points.
    targets: Literal["reaction_stationary_points", "all_minima"] = "reaction_stationary_points"


def _targets(inputs: Manifest, which: str) -> dict[str, Evidence]:
    """Subject id -> freq Evidence whose final geometry is computed."""
    minima = {
        m.minimum_id: m for m in inputs.records(ArtifactType.MINIMUM, MinimumRecord)
        if m.tier == "dft"
    }
    out: dict[str, Evidence] = {}
    if which == "all_minima":
        out |= {mid: inputs.evidence(m.freq_calc) for mid, m in minima.items()}
    for reaction in inputs.records(ArtifactType.REACTION, ReactionRecord):
        out |= {mid: inputs.evidence(minima[mid].freq_calc)
                for mid in reaction.minima if mid in minima}
        if reaction.saddle is not None:
            out[reaction.saddle.freq_calc] = inputs.evidence(reaction.saddle.freq_calc)
    return out


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
        consumes=(ArtifactType.MINIMUM,),
        produces=(ArtifactType.CALCULATION,),
    )

    def run(self, inputs: Manifest, config: StageConfig, rt: StageRuntime) -> list[Artifact]:
        cfg = cast(SinglePointConfig, config)
        engine = cast(QMEngine, rt.engine(Capability.QM, cfg.engine))
        subjects = _targets(inputs, cfg.targets)
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

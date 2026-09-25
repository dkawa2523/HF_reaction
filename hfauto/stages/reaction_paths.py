"""reaction-paths stage (design §4.1, §8.2): hypotheses → ReactionCaseDriver, one case at a time.

Engines have no defaults here; the pipeline YAML names them (``engines`` and ``screen``).
"""

from __future__ import annotations

from collections import deque
from dataclasses import fields, replace
from typing import Any, ClassVar, cast

from pydantic import Field, field_validator

from hfauto.backends.protocols import Capability, PathEngine, QMEngine, SaddleRefiner
from hfauto.chemistry.gates import Policy
from hfauto.chemistry.hypotheses import select
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.records import (
    ArtifactType,
    DiscoveryRecord,
    MinimumRecord,
    ReactionRecord,
    SpeciesRecord,
)
from hfauto.drivers.minimum import Registry
from hfauto.drivers.reaction_case.driver import CaseRuntime, drive_case
from hfauto.drivers.reaction_case.state import CasePolicy
from hfauto.stages.spec import StageConfig, StageRuntime, StageSpec

_POLICY_KEYS = frozenset(f.name for f in fields(CasePolicy)) - {"gates"} | {"walltime_h"}


class EngineNames(StageConfig):
    qm: str
    saddle: str
    path: str


class ScreenConfig(StageConfig):
    method: str
    qm: str
    path: str


class ReactionPathsConfig(StageConfig):
    method: str
    engines: EngineNames
    screen: ScreenConfig | None = None  # None: no barrier pre-check (CasePolicy.screen off)
    policy: dict[str, Any] = Field(default_factory=dict)  # CasePolicy overrides, or walltime_h
    reaction_ids: list[str] | None = None  # debugging: only these hypotheses

    @field_validator("policy")
    @classmethod
    def _known_keys(cls, value: dict[str, Any]) -> dict[str, Any]:
        unknown = sorted(set(value) - _POLICY_KEYS)
        if unknown:
            raise ValueError(f"unknown policy keys {unknown}; allowed: {sorted(_POLICY_KEYS)}")
        return value

    def case_policy(self, gates: Policy) -> CasePolicy:
        values = dict(self.policy)
        if "walltime_h" in values:
            values["walltime_s"] = float(values.pop("walltime_h")) * 3600.0
        if "qrc_bounds_A" in values:
            values["qrc_bounds_A"] = tuple(values["qrc_bounds_A"])
        policy = CasePolicy(**values, gates=gates)
        return replace(policy, screen=policy.screen and self.screen is not None)


def _case_runtime(config: ReactionPathsConfig, rt: StageRuntime, inputs: Manifest,
                  species: dict[str, SpeciesRecord]) -> CaseRuntime:
    """Engines, methods and the DFT minima registry of this stage."""
    dft = [m for m in inputs.records(ArtifactType.MINIMUM, MinimumRecord) if m.tier == "dft"]
    minima = {m.minimum_id: (m, inputs.evidence(m.opt_calc).final) for m in dft}
    run_dir = rt.stage_dir.parent  # RunLayout: <run>/<stage_id>/
    screen = config.screen
    return CaseRuntime(
        qm=cast(QMEngine, rt.engine(Capability.QM, config.engines.qm)),
        saddle=cast(SaddleRefiner, rt.engine(Capability.SADDLE, config.engines.saddle)),
        path=cast(PathEngine, rt.engine(Capability.PATH, config.engines.path)),
        screen_qm=cast(QMEngine, rt.engine(Capability.QM, screen.qm)) if screen else None,
        screen_path=cast(PathEngine, rt.engine(Capability.PATH, screen.path)) if screen else None,
        method=rt.method(config.method),
        screen_method=rt.method(screen.method) if screen else None,
        registry=Registry(minima.values(), rt.load_xyz),
        load_xyz=rt.load_xyz,
        file_ref=rt.file_ref,
        case_dir=rt.stage_dir / "cases",
        deadline=rt.deadline,
        resolve=lambda ref: run_dir / ref.path,
        minima=minima,
        species=species,
    )


class ReactionPathsStage:
    spec: ClassVar[StageSpec] = StageSpec(
        name="reaction-paths",
        config=ReactionPathsConfig,
        consumes=(ArtifactType.MINIMUM, ArtifactType.SPECIES),
        produces=(ArtifactType.REACTION, ArtifactType.CALCULATION, ArtifactType.MINIMUM,
                  ArtifactType.SPECIES),
    )

    def run(self, inputs: Manifest, config: StageConfig, rt: StageRuntime) -> list[Artifact]:
        if not isinstance(config, ReactionPathsConfig):
            raise TypeError(f"expected ReactionPathsConfig, got {type(config).__name__}")
        policy = config.case_policy(rt.policy)
        species = {s.species_id: s for s in inputs.records(ArtifactType.SPECIES, SpeciesRecord)}
        case_rt = _case_runtime(config, rt, inputs, species)
        cases = select(
            inputs.records(ArtifactType.MINIMUM, MinimumRecord), list(species.values()),
            inputs.records(ArtifactType.DISCOVERY, DiscoveryRecord), rt.system.reactions,
            rt.load_xyz, window_kcal=rt.policy.reaction_window_kcal,
        )
        if config.reaction_ids is not None:
            cases = [c for c in cases if c.reaction_id in config.reaction_ids]
        queue: deque[tuple[ReactionRecord, int]] = deque((case, 0) for case in cases)
        artifacts: dict[str, Artifact] = {}  # a job shared by cases is emitted once
        while queue:  # serial (§7.1); split children follow up to max_split_depth
            case, depth = queue.popleft()
            result = drive_case(case, case_rt, policy)
            artifacts.update((a.artifact_id, a) for a in result.artifacts)
            if depth < policy.max_split_depth:
                queue.extend((child, depth + 1) for child in result.children)
        return list(artifacts.values())

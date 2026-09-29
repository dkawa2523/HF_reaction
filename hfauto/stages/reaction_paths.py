"""reaction-paths stage (design §4.1, §8.2): hypotheses → ReactionCaseDriver, one case at a time.

Engines have no defaults here; the pipeline YAML names them (``engines`` and ``screen``). A
hypothesis has ``policy.walltime_h`` in all: its split children run right after it, on its
deadline, and see its calculations (a TS it validated for a child).
"""

from __future__ import annotations

import json
import os
from collections.abc import Iterable
from typing import ClassVar, cast

from hfauto.backends.protocols import Capability, PathEngine, QMEngine, SaddleRefiner
from hfauto.chemistry.hypotheses import select
from hfauto.core.evidence import Evidence
from hfauto.core.ids import path_token
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import Deadline
from hfauto.core.records import (
    ArtifactType,
    CaseOutcome,
    DiscoveryRecord,
    MinimumRecord,
    ReactionRecord,
    SpeciesRecord,
)
from hfauto.drivers.minimum import Registry
from hfauto.drivers.reaction_case.driver import CaseResult, CaseRuntime, drive_case
from hfauto.drivers.reaction_case.state import CaseRules, ReactionPathsPolicy
from hfauto.stages.spec import StageConfig, StageRuntime, StageSpec


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
    screen: ScreenConfig | None = None  # None: no barrier pre-check (no SCREEN row)
    policy: ReactionPathsPolicy = ReactionPathsPolicy()
    reaction_ids: list[str] | None = None  # debugging: only these hypotheses


def _case_runtime(config: ReactionPathsConfig, rt: StageRuntime, inputs: Manifest,
                  species: dict[str, SpeciesRecord]) -> CaseRuntime:
    """Engines, methods, the DFT minima registry and the calculations of this stage."""
    dft = [m for m in inputs.records(ArtifactType.MINIMUM, MinimumRecord) if m.tier == "dft"]
    minima = {m.minimum_id: (m, inputs.evidence(m.opt_calc).final) for m in dft}
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
        resolve=rt.resolve,
        map=rt.thread_map,
        minima=minima,
        species=species,
        calcs=_calcs(inputs.of(ArtifactType.CALCULATION)),
    )


def _calcs(artifacts: Iterable[Artifact]) -> dict[str, Evidence]:
    return {a.artifact_id: a.payload for a in artifacts if isinstance(a.payload, Evidence)}


def _drive(case: ReactionRecord, rt: CaseRuntime, rules: CaseRules, deadline: Deadline
           ) -> CaseResult:
    """drive_case; an exception leaves only this case UNRESOLVED (``error:<type>``)."""
    try:
        return drive_case(case, rt, rules, deadline)
    except Exception as exc:
        if os.environ.get("HFAUTO_STRICT") == "1":
            raise
        reason = f"error:{type(exc).__name__}"
        folder = rt.case_dir / path_token(case.reaction_id)
        folder.mkdir(parents=True, exist_ok=True)
        entry = {"action": "error", "reason": reason, "detail": str(exc)[:500]}
        with (folder / "log.jsonl").open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(entry, sort_keys=True) + "\n")
        record = case.model_copy(update={"outcome": CaseOutcome.UNRESOLVED, "reasons": (reason,),
                                         "log": f"cases/{folder.name}/log.jsonl"})
        artifact = Artifact(artifact_id=record.reaction_id, type=ArtifactType.REACTION,
                            payload=record, parents=tuple(m for m in record.minima if m))
        return CaseResult(record, (), (artifact,))


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
        budget = config.policy
        rules = CaseRules(gates=rt.policy, budget=budget, screen=config.screen is not None)
        species = {s.species_id: s for s in inputs.records(ArtifactType.SPECIES, SpeciesRecord)}
        case_rt = _case_runtime(config, rt, inputs, species)
        minima = [(m, inputs.evidence(m.opt_calc).final)
                  for m in inputs.records(ArtifactType.MINIMUM, MinimumRecord)]
        cases = select(
            minima, list(species.values()),
            inputs.records(ArtifactType.DISCOVERY, DiscoveryRecord), rt.system.reactions,
            rt.load_xyz, window_kcal=rt.policy.reaction_window_kcal,
        )
        if config.reaction_ids is not None:
            cases = [c for c in cases if c.reaction_id in config.reaction_ids]
        # serial (§7.1), a stack: split children (up to max_split_depth) run next, on the
        # deadline of their hypothesis; a hypothesis' own deadline starts when it does
        stack: list[tuple[ReactionRecord, int, Deadline | None]] = [
            (case, 0, None) for case in reversed(cases)]
        artifacts: dict[str, Artifact] = {}  # a job shared by cases is emitted once
        while stack:
            case, depth, deadline = stack.pop()
            deadline = deadline or rt.deadline(3600.0 * budget.walltime_h)
            result = _drive(case, case_rt, rules, deadline)
            artifacts.update((a.artifact_id, a) for a in result.artifacts)
            case_rt.calcs.update(_calcs(result.artifacts))
            if depth < budget.max_split_depth:
                stack += [(child, depth + 1, deadline) for child in reversed(result.children)]
        return list(artifacts.values())

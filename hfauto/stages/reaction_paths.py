"""reaction-paths stage (design §4.1, §8.2): hypotheses → ReactionCaseDriver, one case at a time.

The monomer states of the system's compositions (thermo.monomer_states) make an association's
separated reactant side (chemistry.hypotheses). Engines have no defaults here; the pipeline YAML
names them (``engines`` and ``screen``). A hypothesis' split children run right after it, each
with its own ``policy`` budget, and see its calculations (a TS it validated for a child). A split
child is not driven when a case of its case key (hypotheses.pair_key), run or queued, concludes:
it takes that conclusion, ``same_as:<reaction_id>``; a case that ends UNRESOLVED gives none, and
the child is driven. Nor is a child deeper than ``max_split_depth``: UNRESOLVED ``split_depth``.
Neither has a job or a log.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from typing import ClassVar, cast

from hfauto.backends.protocols import Capability, PathEngine, QMEngine, SaddleRefiner
from hfauto.chemistry.classification import undriven
from hfauto.chemistry.hypotheses import CaseKey, pair_key, select
from hfauto.chemistry.thermo import monomer_states
from hfauto.core.evidence import Evidence
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.records import (
    ArtifactType,
    CaseOutcome,
    DiscoveryRecord,
    MinimumRecord,
    ReactionRecord,
    SpeciesRecord,
)
from hfauto.drivers.minimum import Registry
from hfauto.drivers.reaction_case.driver import (
    CaseResult,
    CaseRuntime,
    drive_case,
    reaction_artifact,
)
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
    dft = [(m, inputs.evidence(m.opt_calc).final)
           for m in inputs.records(ArtifactType.MINIMUM, MinimumRecord) if m.tier == "dft"]
    screen = config.screen
    return CaseRuntime(
        qm=cast(QMEngine, rt.engine(Capability.QM, config.engines.qm)),
        saddle=cast(SaddleRefiner, rt.engine(Capability.SADDLE, config.engines.saddle)),
        path=cast(PathEngine, rt.engine(Capability.PATH, config.engines.path)),
        screen_qm=cast(QMEngine, rt.engine(Capability.QM, screen.qm)) if screen else None,
        screen_path=cast(PathEngine, rt.engine(Capability.PATH, screen.path)) if screen else None,
        method=rt.method(config.method),
        screen_method=rt.method(screen.method) if screen else None,
        registry=Registry(dft, rt.load_xyz),
        load_xyz=rt.load_xyz,
        file_ref=rt.file_ref,
        case_dir=rt.stage_dir / "cases",
        resolve=rt.resolve,
        map=rt.thread_map,
        species=species,
        calcs=_calcs(inputs.of(ArtifactType.CALCULATION)),
    )


def _calcs(artifacts: Iterable[Artifact]) -> dict[str, Evidence]:
    return {a.artifact_id: a.payload for a in artifacts if isinstance(a.payload, Evidence)}


Item = tuple[ReactionRecord, int]  # a case and its split depth


@dataclass
class _Book:
    """The cases of the stage and what they emit (a job shared by cases is emitted once). A
    split child whose case key a case has takes that case's conclusion, waiting for a queued
    one; a failure to conclude (UNRESOLVED) is no result to share, so the child is then driven,
    as a child of a new key is, up to ``max_depth`` (deeper: UNRESOLVED)."""

    rt: CaseRuntime
    max_depth: int
    first: dict[CaseKey, str] = field(default_factory=dict)  # key -> its first case to drive
    done: dict[str, ReactionRecord] = field(default_factory=dict)
    waiting: list[Item] = field(default_factory=list)  # children of a key still queued
    artifacts: dict[str, Artifact] = field(default_factory=dict)

    def key(self, case: ReactionRecord) -> CaseKey | None:
        """The case key of a case between two DFT minima; None for a blocked one."""
        a, b = (self.rt.registry.minima.get(m) for m in case.minima)
        return None if a is None or b is None else pair_key(a[0], b[0])

    def queue(self, cases: Sequence[ReactionRecord]) -> None:
        """Cases to drive: each takes its key unless an earlier one has it."""
        for case in cases:
            if (key := self.key(case)) is not None:
                self.first.setdefault(key, case.reaction_id)

    def emit(self, result: CaseResult, depth: int) -> list[Item]:
        """A driven case's artifacts; the split children of ``result`` to drive now."""
        self.done[result.reaction.reaction_id] = result.reaction
        self.artifacts.update((a.artifact_id, a) for a in result.artifacts)
        self.rt.calcs.update(_calcs(result.artifacts))
        return [i for child in result.children for i in self._child((child, depth + 1))]

    def waited(self) -> list[Item]:
        """The waiting children, once every queued case has been driven."""
        items, self.waiting = self.waiting, []
        return [i for item in items for i in self._child(item)]

    def _child(self, item: Item) -> list[Item]:
        child, depth = item
        key = self.key(child)
        case_id = None if key is None else self.first.get(key)
        if case_id is not None:
            like = self.done.get(case_id)
            if like is None:
                self.waiting.append(item)
                return []
            if like.outcome is not CaseOutcome.UNRESOLVED:
                return self._record(undriven(child, f"same_as:{case_id}", like))
        if depth > self.max_depth:
            return self._record(undriven(child, "split_depth"))
        self.queue([child])
        return [item]

    def _record(self, record: ReactionRecord) -> list[Item]:
        self.artifacts[record.reaction_id] = reaction_artifact(record)
        return []


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
        rules = CaseRules(gates=rt.policy, budget=config.policy, screen=config.screen is not None)
        species = {s.species_id: s for s in inputs.records(ArtifactType.SPECIES, SpeciesRecord)}
        case_rt = _case_runtime(config, rt, inputs, species)
        opts = [(m, inputs.evidence(m.opt_calc))
                for m in inputs.records(ArtifactType.MINIMUM, MinimumRecord)]
        cases = select(
            [(m, opt.final) for m, opt in opts], list(species.values()),
            inputs.records(ArtifactType.DISCOVERY, DiscoveryRecord), rt.system.reactions,
            rt.load_xyz, window_kcal=rt.policy.reaction_window_kcal,
            monomers=monomer_states(species.values(), rt.system.compositions),
            levels={m.minimum_id: opt.level for m, opt in opts if m.tier == "dft"},
        )
        if config.reaction_ids is not None:
            cases = [c for c in cases if c.reaction_id in config.reaction_ids]
        # serial (§7.1), a stack: split children run next (one waiting for a queued case of
        # its key: once the stack is empty)
        book = _Book(case_rt, config.policy.max_split_depth)
        book.queue(cases)
        stack: list[Item] = [(case, 0) for case in reversed(cases)]
        while stack:
            case, depth = stack.pop()
            stack += reversed(book.emit(drive_case(case, case_rt, rules), depth))
            stack = stack or list(reversed(book.waited()))
        return list(book.artifacts.values())

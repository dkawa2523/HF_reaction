"""Stage execution (design §7.4): ``Runtime``, ``execute_stage`` and ``run_pipeline``."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable, Collection, Mapping, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import TYPE_CHECKING, Literal, TypeVar

from hfauto.chemistry.gates import Policy
from hfauto.chemistry.xyz import XYZ, read_xyz
from hfauto.core.evidence import FileRef, Geometry
from hfauto.core.hashing import sha256_file
from hfauto.core.manifest import Manifest, save_manifest
from hfauto.core.method import Deadline, MethodSpec
from hfauto.core.system import Conditions, SystemConfig
from hfauto.pipeline.config import ResolvedConfig, SiteConfig, StageEntry, method_ids
from hfauto.pipeline.layout import JobCounts, RunLayout, now
from hfauto.stages import catalog

if TYPE_CHECKING:
    from hfauto.backends.protocols import Capability, Engine
    from hfauto.execution.jobs import JobRunner, JobStats

ItemT = TypeVar("ItemT")
ResultT = TypeVar("ResultT")


@dataclass
class Runtime:
    """``StageRuntime`` implementation; engines are cached per (capability, name) per run."""

    run_id: str
    run_dir: Path
    system: SystemConfig
    conditions: Conditions
    policy: Policy
    site: SiteConfig
    methods: Mapping[str, MethodSpec]
    jobs: JobRunner
    create: Callable[..., Engine]  # backends.engines.create
    stage_id: str = ""
    stage_dir: Path = Path()
    engines: dict[tuple[Capability, str], Engine] = field(default_factory=dict)

    def bind(self, stage_id: str, stage_dir: Path) -> Runtime:
        """The same runtime (shared engine cache and JobRunner) for one stage."""
        return replace(self, stage_id=stage_id, stage_dir=stage_dir)

    def engine(self, capability: Capability, name: str) -> Engine:
        key = (capability, name)
        if key not in self.engines:
            if name not in self.site.engines:
                raise KeyError(f"engine {name!r} is not configured in site {self.site.site!r}")
            site = self.site.engines[name]
            self.engines[key] = self.create(capability, name, jobs=self.jobs, site=site)
        return self.engines[key]

    def method(self, method_id: str) -> MethodSpec:
        if method_id not in self.methods:
            raise KeyError(f"method {method_id!r} is not referenced by the pipeline")
        return self.methods[method_id]

    def load_xyz(self, geometry: Geometry) -> XYZ:
        path = self.run_dir / geometry.file.path
        if sha256_file(path) != geometry.file.sha256:
            raise ValueError(f"{geometry.file.path}: sha256 differs from the recorded FileRef")
        return read_xyz(path)

    def file_ref(self, path: Path) -> FileRef:
        resolved = Path(path).resolve()
        rel = resolved.relative_to(self.run_dir)  # ValueError outside the run directory
        return FileRef(path=rel.as_posix(), sha256=sha256_file(resolved))

    def thread_map(
        self, fn: Callable[[ItemT], ResultT], items: Sequence[ItemT], *, threads_per_item: int
    ) -> list[ResultT]:
        from hfauto.execution.jobs import thread_map

        workers = max(1, self.site.cores // max(1, threads_per_item))
        return thread_map(fn, items, workers=workers)

    def deadline(self, seconds: float) -> Deadline:
        return Deadline.after(seconds)


def build_runtime(
    resolved: ResolvedConfig, layout: RunLayout, *, retry_failed: Collection[str] = ()
) -> Runtime:
    from hfauto.backends.engines import create
    from hfauto.execution.jobs import JobRunner
    from hfauto.execution.jobstore import JobStore

    store = JobStore(layout.jobs_dir, retry_failed=retry_failed)
    return Runtime(
        run_id=layout.run_id,
        run_dir=layout.run_dir,
        system=resolved.system,
        conditions=resolved.pipeline.conditions,
        policy=resolved.policy(),
        site=resolved.site,
        methods=resolved.methods,
        jobs=JobRunner(store, cores=resolved.site.cores),
        create=create,
    )


def config_sha(entry: StageEntry, resolved: ResolvedConfig) -> str:
    """What a stage's result depends on besides its inputs (xyz files by content)."""
    system = resolved.system.model_dump(mode="json")
    for species, raw in zip(resolved.system.species, system["species"], strict=True):
        raw["xyz"] = sha256_file(species.xyz) if species.xyz is not None else None
    methods = {m: resolved.methods[m].model_dump(mode="json") for m in method_ids(entry.settings())}
    data = {
        "entry": entry.model_dump(mode="json"),
        "conditions": resolved.pipeline.conditions.model_dump(mode="json"),
        "gates": resolved.pipeline.gates,
        "system": system,
        "methods": methods,
    }
    text = json.dumps(data, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _job_counts(before: JobStats, after: JobStats) -> JobCounts:
    failures = {
        kind: n - before.failures_by_kind.get(kind, 0)
        for kind, n in after.failures_by_kind.items()  # counters only grow
    }
    return JobCounts(
        hits=after.hits - before.hits,
        misses=after.misses - before.misses,
        failures_by_kind={kind: n for kind, n in failures.items() if n},
    )


def execute_stage(
    entry: StageEntry, resolved: ResolvedConfig, layout: RunLayout, runtime: Runtime
) -> Manifest:
    """The only way a stage runs: validate, check consumes (a missing type names up to 3 failed
    upstream artifacts and why), run, check produces, save, record."""
    pipeline_id = resolved.pipeline.pipeline_id
    layout.begin(entry.id, pipeline_id)
    before = runtime.jobs.stats()
    try:
        stage_cls = catalog.get(entry.stage)
        spec = stage_cls.spec
        config = spec.config.model_validate(entry.settings())
        inputs = layout.view(entry.id)
        missing = [t.value for t in spec.consumes if not inputs.of(t)]
        if missing:
            failed = sorted({f"{a.artifact_id}: {a.failure.reason}" for t in spec.consumes
                             for a in inputs.of(t, ok_only=False) if a.failure is not None})[:3]
            why = f"; upstream failures: {', '.join(failed)}" if failed else ""
            raise ValueError(f"stage {entry.id!r} ({spec.name}) has no input of type {missing}{why}")
        stage_dir = layout.stage_dir(entry.id)
        stage_dir.mkdir(parents=True, exist_ok=True)
        artifacts = stage_cls().run(inputs, config, runtime.bind(entry.id, stage_dir))
        unexpected = sorted({a.type.value for a in artifacts} - {t.value for t in spec.produces})
        if unexpected:
            raise ValueError(f"stage {entry.id!r} ({spec.name}) produced {unexpected}")
        manifest = Manifest(
            run_id=layout.run_id, stage_id=entry.id, created_at=now(), artifacts=list(artifacts)
        )
        save_manifest(manifest, layout.manifest_path(entry.id))
    except BaseException:
        jobs = _job_counts(before, runtime.jobs.stats())
        layout.update(entry.id, status="failed", finished=now(), jobs=jobs)
        raise
    layout.update(
        entry.id,
        status="done",
        finished=now(),
        input_sha=layout.input_sha(entry.id),
        config_sha=config_sha(entry, resolved),
        n_ok=sum(a.status == "success" for a in artifacts),
        n_failed=sum(a.status == "failed" for a in artifacts),
        jobs=_job_counts(before, runtime.jobs.stats()),
    )
    return manifest


@dataclass(frozen=True)
class PlannedStage:
    stage_id: str
    stage: str
    action: Literal["run", "skip"]


def _select(entries: list[StageEntry], start: str | None, stop: str | None) -> list[StageEntry]:
    ids = [e.id for e in entries]
    for name, value in (("--from", start), ("--to", stop)):
        if value is not None and value not in ids:
            raise ValueError(f"unknown {name} stage {value!r}; stages: {ids}")
    first = ids.index(start) if start is not None else 0
    last = ids.index(stop) if stop is not None else len(ids) - 1
    if first > last:
        raise ValueError(f"--from {start!r} comes after --to {stop!r}")
    return entries[first : last + 1]


def _fresh(
    entry: StageEntry, resolved: ResolvedConfig, layout: RunLayout, retry_failed: Collection[str]
) -> bool:
    """Done with unchanged inputs and config, and no job failure of a kind to retry."""
    state = layout.state(entry.id)
    return (
        state is not None
        and state.status == "done"
        and not set(state.jobs.failures_by_kind) & set(retry_failed)
        and state.input_sha == layout.input_sha(entry.id)
        and state.config_sha == config_sha(entry, resolved)
    )


def plan(
    resolved: ResolvedConfig,
    layout: RunLayout,
    *,
    start: str | None = None,
    stop: str | None = None,
    retry_failed: Collection[str] = (),
) -> list[PlannedStage]:
    """What ``run_pipeline`` does; every stage after a re-run stage runs again."""
    pipeline_id = resolved.pipeline.pipeline_id
    for entry in resolved.pipeline.stages:
        state = layout.state(entry.id)
        if state is not None and state.pipeline_id != pipeline_id:
            raise ValueError(f"stage id {entry.id!r} already belongs to {state.pipeline_id!r}")
    rerun = False
    planned: list[PlannedStage] = []
    for entry in _select(resolved.pipeline.stages, start, stop):
        catalog.get(entry.stage).spec.config.model_validate(entry.settings())
        rerun = rerun or entry.id == start or not _fresh(entry, resolved, layout, retry_failed)
        planned.append(PlannedStage(entry.id, entry.stage, "run" if rerun else "skip"))
    return planned


def run_pipeline(
    resolved: ResolvedConfig,
    run_dir: Path,
    *,
    start: str | None = None,
    stop: str | None = None,
    dry_run: bool = False,
    retry_failed: Collection[str] = (),
) -> RunLayout | list[PlannedStage]:
    """Run the selected stages while holding the SiteLock; ``dry_run`` only returns the plan.

    ``start`` (--from) marks that stage and everything executed after it stale; done stages
    whose input_sha and config_sha are unchanged (and whose jobs had no failure of a
    ``retry_failed`` kind) are skipped (resume).
    """
    layout = RunLayout(run_dir)
    if dry_run:
        return plan(resolved, layout, start=start, stop=stop, retry_failed=retry_failed)
    from hfauto.execution.lock import SiteLock

    entries = {entry.id: entry for entry in resolved.pipeline.stages}
    with SiteLock(Path(resolved.site.scratch_root), layout.run_dir):
        layout.run_dir.mkdir(parents=True, exist_ok=True)
        for entry in entries.values():  # stage id collisions fail before any run
            layout.ensure(entry.id, resolved.pipeline.pipeline_id)
        layout.write_resolved_config(resolved)
        if start is not None:
            layout.mark_stale(start)
        steps = plan(resolved, layout, start=start, stop=stop, retry_failed=retry_failed)
        runtime = build_runtime(resolved, layout, retry_failed=retry_failed)
        for step in steps:
            if step.action == "run":
                execute_stage(entries[step.stage_id], resolved, layout, runtime)
    return layout

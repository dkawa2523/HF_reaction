"""Run directory layout, run_state.json and stage input views (design §7.4).

``<run>/<stage_id>/{manifest.json,diagnostics.json,cases/}``, ``<run>/jobs/``,
``<run>/run_state.json`` and ``<run>/resolved_config.yaml``. run_state is the list of
stages in the order they were (last) executed; a stage's input is the union of the
manifests of the done stages before it (so another pipeline can append to a run).
Pending entries have not run yet: their position carries no meaning.
"""

from __future__ import annotations

import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import BaseModel, ConfigDict, TypeAdapter

from hfauto.core.files import write_atomic
from hfauto.core.hashing import sha256_file
from hfauto.core.manifest import Manifest, load_manifest
from hfauto.execution.jobs import JobStats
from hfauto.pipeline.config import ResolvedConfig

Status = Literal["pending", "running", "done", "failed", "stale"]


class StageState(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")
    stage_id: str
    pipeline_id: str
    status: Status = "pending"
    input_sha: str | None = None
    config_sha: str | None = None
    started: str | None = None
    finished: str | None = None
    n_ok: int = 0
    n_failed: int = 0
    jobs: JobStats = JobStats()  # of this stage's last execution


_STATES = TypeAdapter(list[StageState])


def now() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds")


def _index(states: list[StageState], stage_id: str) -> int:
    for i, state in enumerate(states):
        if state.stage_id == stage_id:
            return i
    raise KeyError(f"stage {stage_id!r} is not in run_state")


class RunLayout:
    def __init__(self, run_dir: Path) -> None:
        self.run_dir = Path(run_dir).resolve()

    @property
    def run_id(self) -> str:
        return self.run_dir.name

    @property
    def jobs_dir(self) -> Path:
        return self.run_dir / "jobs"

    @property
    def run_state_path(self) -> Path:
        return self.run_dir / "run_state.json"

    @property
    def resolved_config_path(self) -> Path:
        return self.run_dir / "resolved_config.yaml"

    def stage_dir(self, stage_id: str) -> Path:
        return self.run_dir / stage_id

    def manifest_path(self, stage_id: str) -> Path:
        return self.stage_dir(stage_id) / "manifest.json"

    def cases_dir(self, stage_id: str) -> Path:
        return self.stage_dir(stage_id) / "cases"

    # -- run_state -----------------------------------------------------------------
    def read_state(self) -> list[StageState]:
        if not self.run_state_path.exists():
            return []
        return _STATES.validate_json(self.run_state_path.read_bytes())

    def write_state(self, states: list[StageState]) -> None:
        write_atomic(self.run_state_path, _STATES.dump_json(states, indent=2).decode("utf-8"))

    def state(self, stage_id: str) -> StageState | None:
        return next((s for s in self.read_state() if s.stage_id == stage_id), None)

    def ensure(self, stage_id: str, pipeline_id: str) -> StageState:
        """The stage's entry, appended as pending when new; another pipeline's id is an error."""
        states = self.read_state()
        existing = next((s for s in states if s.stage_id == stage_id), None)
        if existing is not None:
            if existing.pipeline_id != pipeline_id:
                raise ValueError(
                    f"stage id {stage_id!r} already belongs to pipeline"
                    f" {existing.pipeline_id!r} in {self.run_dir}"
                )
            return existing
        entry = StageState(stage_id=stage_id, pipeline_id=pipeline_id)
        self.write_state([*states, entry])
        return entry

    def begin(self, stage_id: str, pipeline_id: str) -> StageState:
        """Mark the stage running and move it to the end of the execution order.

        When the stage was done, every done stage executed after it read its previous
        output and becomes stale (a stage that is not done is in nobody's view).
        """
        self.ensure(stage_id, pipeline_id)
        states = self.read_state()
        index = _index(states, stage_id)
        tail = states[index + 1 :]
        if states[index].status == "done":
            tail = [
                s.model_copy(update={"status": "stale"}) if s.status == "done" else s for s in tail
            ]
        running = StageState(
            stage_id=stage_id, pipeline_id=pipeline_id, status="running", started=now()
        )
        self.write_state([*states[:index], *tail, running])
        return running

    def update(self, stage_id: str, **changes: Any) -> StageState:
        states = self.read_state()
        index = _index(states, stage_id)
        states[index] = StageState.model_validate({**states[index].model_dump(), **changes})
        self.write_state(states)
        return states[index]

    def mark_stale(self, stage_id: str) -> None:
        """``--from``: the stage and every stage executed after it become stale."""
        states = self.read_state()
        index = next(
            (i for i, s in enumerate(states) if s.stage_id == stage_id and s.status != "pending"),
            len(states),
        )
        tail = [
            s if s.status == "pending" else s.model_copy(update={"status": "stale"})
            for s in states[index:]
        ]
        self.write_state(states[:index] + tail)

    # -- views ---------------------------------------------------------------------
    def _inputs(self, stage_id: str | None) -> list[str]:
        """Done stages executed before ``stage_id`` (all of them when it has not run yet)."""
        done: list[str] = []
        for state in self.read_state():
            if state.stage_id == stage_id and state.status != "pending":
                break
            if state.status == "done":
                done.append(state.stage_id)
        return done

    def view(self, stage_id: str | None = None) -> Manifest:
        manifests = [load_manifest(self.manifest_path(s)) for s in self._inputs(stage_id)]
        return Manifest.union(manifests, run_id=self.run_id, stage_id=stage_id or "view")

    def input_sha(self, stage_id: str) -> str:
        """sha256 over the manifests that make up ``view(stage_id)``, in order."""
        parts = [[s, sha256_file(self.manifest_path(s))] for s in self._inputs(stage_id)]
        return hashlib.sha256(json.dumps(parts).encode("utf-8")).hexdigest()

    def write_resolved_config(self, resolved: ResolvedConfig) -> Path:
        """Record the resolved config under its pipeline_id, keeping other pipelines' records."""
        path = self.resolved_config_path
        records = yaml.safe_load(path.read_text(encoding="utf-8")) if path.exists() else None
        records = dict(records or {})
        records[resolved.pipeline.pipeline_id] = resolved.model_dump(mode="json")
        write_atomic(path, yaml.safe_dump(records, sort_keys=False, allow_unicode=True))
        return path

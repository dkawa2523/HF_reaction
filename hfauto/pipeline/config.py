"""Pipeline, site and method files (design §7.4, §7.5).

The four layers (site / method / system / pipeline) do not overlap, so nothing is merged
and there is no precedence order. Stage-specific keys are validated by the stage's own
``StageConfig`` when the stage runs.
"""

from __future__ import annotations

import dataclasses
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from hfauto.chemistry.gates import Policy
from hfauto.core.method import EngineSite, MethodSpec
from hfauto.core.system import Conditions, SystemConfig, load_system
from hfauto.execution.process import Command, resolve_executable, run_command

_POLICY_FIELDS = frozenset(f.name for f in dataclasses.fields(Policy))
_REPO_ROOT = Path(__file__).resolve().parents[2]
_GIT_TIMEOUT_S = 30.0


class SiteConfig(BaseModel):
    """configs/sites/*.yaml: the only place for absolute paths; keyed by registry name."""

    model_config = ConfigDict(frozen=True, extra="forbid")
    site: str
    scratch_root: str
    cores: int = Field(ge=1)
    memory_mb: int = Field(ge=1)
    engines: dict[str, EngineSite] = {}


class StageEntry(BaseModel):
    """One ``stages:`` item; keys other than id / stage belong to the stage's config."""

    model_config = ConfigDict(frozen=True, extra="allow")
    id: str
    stage: str

    def settings(self) -> dict[str, Any]:
        return dict(self.model_extra or {})


class PipelineConfig(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")
    pipeline_id: str
    conditions: Conditions = Conditions()
    gates: dict[str, float] = {}  # Policy overrides
    stages: list[StageEntry]

    @field_validator("gates")
    @classmethod
    def _known_gates(cls, value: dict[str, float]) -> dict[str, float]:
        unknown = sorted(set(value) - _POLICY_FIELDS)
        if unknown:
            raise ValueError(f"unknown gates {unknown}; allowed: {sorted(_POLICY_FIELDS)}")
        return value

    @model_validator(mode="after")
    def _check_stages(self) -> PipelineConfig:
        counts = Counter(entry.id for entry in self.stages)
        duplicated = sorted(k for k, n in counts.items() if n > 1)
        if duplicated:
            raise ValueError(f"duplicate stage ids: {duplicated}")
        # DFT minima and reaction paths must share one stationary-point level of theory
        stationary = {e.id: e.settings().get("method") for e in self.stages
                      if e.stage == "reaction-paths"
                      or (e.stage == "minima" and e.settings().get("level") == "dft")}
        if len(set(stationary.values())) > 1:
            raise ValueError(f"stationary-point stages use different methods: {stationary}")
        return self

    def entry(self, stage_id: str) -> StageEntry:
        for entry in self.stages:
            if entry.id == stage_id:
                return entry
        raise KeyError(f"stage {stage_id!r} is not in pipeline {self.pipeline_id!r}")


class ResolvedConfig(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")
    pipeline: PipelineConfig
    system: SystemConfig
    site: SiteConfig
    methods: dict[str, MethodSpec]
    code_version: str  # "<git sha>", "<git sha>-dirty" or "unknown"

    def policy(self) -> Policy:
        return Policy(**self.pipeline.gates)


def method_ids(value: object) -> set[str]:
    """Method ids below ``value`` under ``method`` / ``energy_method`` / ``methods`` keys."""
    if isinstance(value, list):
        return set().union(*(method_ids(item) for item in value))
    if not isinstance(value, dict):
        return set()
    found: set[str] = set()
    for key, item in value.items():
        if key in ("method", "energy_method") and isinstance(item, str):
            found.add(item)
        elif key == "methods" and isinstance(item, list):
            found.update(str(i) for i in item)
        else:
            found |= method_ids(item)
    return found


def _read_yaml(path: Path) -> Any:
    return yaml.safe_load(Path(path).read_text(encoding="utf-8"))


def load_site(path: Path) -> SiteConfig:
    return SiteConfig.model_validate(_read_yaml(path))


def load_method(methods_dir: Path, method_id: str) -> MethodSpec:
    path = Path(methods_dir) / f"{method_id}.yaml"
    method = MethodSpec.model_validate(_read_yaml(path))
    if method.id != method_id:
        raise ValueError(f"{path}: id {method.id!r} does not match the file name")
    return method


def load(pipeline_path: Path, system_path: Path, site_path: Path) -> ResolvedConfig:
    """Read the four layers; methods come from ``<configs>/methods/<id>.yaml``."""
    pipeline_path = Path(pipeline_path)
    pipeline = PipelineConfig.model_validate(_read_yaml(pipeline_path))
    methods_dir = pipeline_path.resolve().parent.parent / "methods"
    ids = sorted(method_ids([entry.settings() for entry in pipeline.stages]))
    return ResolvedConfig(
        pipeline=pipeline,
        system=load_system(Path(system_path)),
        site=load_site(Path(site_path)),
        methods={method_id: load_method(methods_dir, method_id) for method_id in ids},
        code_version=code_version(),
    )


def _git(tmp: Path, git: str, *args: str) -> str | None:
    command = Command(argv=(git, "-C", str(_REPO_ROOT), *args), cwd=tmp)
    result = run_command(command, timeout_s=_GIT_TIMEOUT_S)
    if result.returncode != 0:
        return None
    return result.stdout.read_text(encoding="utf-8", errors="replace").strip()


def code_version() -> str:
    """The git sha of the hfauto checkout plus "-dirty" for tracked changes; else "unknown"."""
    git = resolve_executable("git")
    if git is None:
        return "unknown"
    with tempfile.TemporaryDirectory(prefix="hfauto_git_") as tmp:
        sha = _git(Path(tmp), git, "rev-parse", "HEAD")
        status = _git(Path(tmp), git, "status", "--porcelain", "--untracked-files=no")
    if not sha or status is None:
        return "unknown"
    return f"{sha}-dirty" if status else sha

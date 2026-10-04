"""Stage types (design §7.4): the only runtime interface a stage may depend on.

The pipeline implements ``StageRuntime``; stages never import the pipeline, the engine
registry or the execution layer.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar, Protocol, TypeVar

from pydantic import BaseModel, ConfigDict

from hfauto.backends.protocols import Capability, Engine
from hfauto.chemistry.gates import Policy
from hfauto.chemistry.xyz import XYZ
from hfauto.core.evidence import FileRef, Geometry
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec
from hfauto.core.records import ArtifactType
from hfauto.core.system import SystemConfig

ItemT = TypeVar("ItemT")
ResultT = TypeVar("ResultT")


class StageConfig(BaseModel):
    """Base of every stage configuration; unknown keys in the pipeline YAML are errors."""

    model_config = ConfigDict(extra="forbid")


@dataclass(frozen=True)
class StageSpec:
    name: str
    config: type[StageConfig]
    consumes: tuple[ArtifactType, ...]  # each type needs at least one successful input artifact
    produces: tuple[ArtifactType, ...]


class StageRuntime(Protocol):
    run_id: str
    stage_id: str
    stage_dir: Path
    system: SystemConfig
    policy: Policy

    def engine(self, capability: Capability, name: str) -> Engine: ...

    def method(self, method_id: str) -> MethodSpec: ...

    def load_xyz(self, geometry: Geometry) -> XYZ: ...

    def file_ref(self, path: Path) -> FileRef: ...  # path relative to the run directory + sha

    def resolve(self, ref: FileRef) -> Path: ...  # the file a FileRef of the run names

    def thread_map(  # in input order; each item's jobs get cores // (items in its wave) ranks
        self, fn: Callable[[ItemT], ResultT], items: Sequence[ItemT]
    ) -> list[ResultT]: ...


class Stage(Protocol):
    spec: ClassVar[StageSpec]

    def run(self, inputs: Manifest, config: StageConfig, rt: StageRuntime) -> list[Artifact]: ...

"""Stage manifests (design §5.3): one file per stage holding only that stage's artifacts."""

from __future__ import annotations

import os
import tempfile
from collections.abc import Sequence
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal, TypeVar

from pydantic import BaseModel, ConfigDict, model_validator

from hfauto.core.evidence import Evidence, Failure
from hfauto.core.records import ArtifactType, Payload

T = TypeVar("T", bound=BaseModel)


class Artifact(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")
    artifact_id: str
    type: ArtifactType
    parents: tuple[str, ...] = ()
    status: Literal["success", "failed"] = "success"
    payload: Payload | None = None
    failure: Failure | None = None

    @model_validator(mode="after")
    def _check_status(self) -> Artifact:
        if self.status == "success":
            if self.payload is None:
                raise ValueError(f"{self.artifact_id}: successful artifact without payload")
            if self.payload.kind != self.type.value:
                raise ValueError(
                    f"{self.artifact_id}: payload kind {self.payload.kind!r}"
                    f" != type {self.type.value!r}"
                )
        elif self.failure is None:
            raise ValueError(f"{self.artifact_id}: failed artifact without failure")
        return self


class Manifest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: Literal["hfauto.manifest.v2"] = "hfauto.manifest.v2"
    run_id: str
    stage_id: str
    created_at: str
    artifacts: list[Artifact] = []

    def of(self, type: ArtifactType, *, ok_only: bool = True) -> list[Artifact]:
        return [
            a for a in self.artifacts
            if a.type == type and (not ok_only or a.status == "success")
        ]

    def get(self, artifact_id: str) -> Artifact:
        """Return the artifact with this id; the last one wins when ids repeat."""
        for artifact in reversed(self.artifacts):
            if artifact.artifact_id == artifact_id:
                return artifact
        raise KeyError(artifact_id)

    def records(self, type: ArtifactType, model: type[T]) -> list[T]:
        out: list[T] = []
        for artifact in self.of(type):
            if not isinstance(artifact.payload, model):
                raise TypeError(f"{artifact.artifact_id}: payload is not {model.__name__}")
            out.append(artifact.payload)
        return out

    def evidence(self, calc_id: str) -> Evidence:
        artifact = self.get(calc_id)
        if artifact.status != "success" or not isinstance(artifact.payload, Evidence):
            raise ValueError(f"{calc_id}: not a successful calculation artifact")
        return artifact.payload

    @classmethod
    def union(cls, manifests: Sequence[Manifest], *, run_id: str, stage_id: str) -> Manifest:
        """Merge artifacts keeping first-seen order; a later artifact with the same id wins."""
        merged: dict[str, Artifact] = {}
        for manifest in manifests:
            for artifact in manifest.artifacts:
                merged[artifact.artifact_id] = artifact
        return cls(
            run_id=run_id,
            stage_id=stage_id,
            created_at=datetime.now(UTC).isoformat(timespec="seconds"),
            artifacts=list(merged.values()),
        )


def load_manifest(path: Path) -> Manifest:
    return Manifest.model_validate_json(Path(path).read_text(encoding="utf-8"))


def save_manifest(manifest: Manifest, path: Path) -> Path:
    """Write atomically: a temporary file in the same directory, then os.replace."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(manifest.model_dump_json(indent=2))
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise
    return path

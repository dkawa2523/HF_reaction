from __future__ import annotations

from collections.abc import Iterable
from datetime import datetime, timezone
from typing import Any
from uuid import uuid4

from pydantic import BaseModel, Field

from hfauto.core.schemas.artifact import Artifact


class Manifest(BaseModel):
    schema_version: str = "hfauto.manifest.v1"
    manifest_id: str = Field(default_factory=lambda: f"man_{uuid4().hex[:12]}")
    run_id: str
    stage: str
    created_at: str = Field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    parents: list[str] = Field(default_factory=list)
    metadata: dict[str, Any] = Field(default_factory=dict)
    artifacts: list[Artifact] = Field(default_factory=list)

    @classmethod
    def new(
        cls,
        run_id: str,
        stage: str,
        parents: list[str] | None = None,
        metadata: dict[str, Any] | None = None,
    ) -> Manifest:
        return cls(run_id=run_id, stage=stage, parents=parents or [], metadata=metadata or {})

    @classmethod
    def merge(
        cls,
        manifests: Iterable[Manifest],
        *,
        run_id: str,
        stage: str,
        metadata: dict[str, Any] | None = None,
    ) -> Manifest:
        """Join independent evidence branches into one downstream manifest.

        Later manifests supply the latest revision of a repeated artifact ID;
        unrelated calculations from every branch remain available to thermo
        and scientific gates.
        """

        sources = list(manifests)
        if not sources:
            raise ValueError("at least one source manifest is required")
        out = cls.new(
            run_id=run_id,
            stage=stage,
            parents=[source.manifest_id for source in sources],
            metadata=metadata,
        )
        for source in sources:
            out.extend(source.latest_artifacts())
        return out

    def add_artifact(self, artifact: Artifact) -> None:
        self.artifacts.append(artifact)

    def extend(self, artifacts: Iterable[Artifact]) -> None:
        self.artifacts.extend(list(artifacts))

    def iter_artifacts(self, artifact_type: str | None = None) -> Iterable[Artifact]:
        for artifact in self.artifacts:
            if artifact_type is None or artifact.artifact_type == artifact_type:
                yield artifact

    def find(self, artifact_id: str) -> Artifact | None:
        """Return the latest artifact with this id.

        Stages may add a revised Artifact with the same id, for example a
        species whose xyz_path was updated by preopt. Returning the latest
        revision keeps CLI inspection and downstream stages aligned.
        """
        for artifact in reversed(self.artifacts):
            if artifact.artifact_id == artifact_id:
                return artifact
        return None

    def latest_artifacts(self, artifact_type: str | None = None) -> list[Artifact]:
        """Return latest revision of each artifact id, preserving first-seen order."""
        order: list[str] = []
        latest: dict[str, Artifact] = {}
        for artifact in self.artifacts:
            if artifact_type is not None and artifact.artifact_type != artifact_type:
                continue
            if artifact.artifact_id not in latest:
                order.append(artifact.artifact_id)
            latest[artifact.artifact_id] = artifact
        return [latest[k] for k in order]

    def carry_forward(self, stage: str) -> Manifest:
        """Create a new manifest carrying current artifacts and run metadata forward."""
        out = Manifest.new(run_id=self.run_id, stage=stage, parents=[self.manifest_id], metadata=dict(self.metadata or {}))
        out.artifacts = list(self.artifacts)
        return out

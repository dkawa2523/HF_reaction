from __future__ import annotations

from datetime import datetime, timezone
from typing import Any

from pydantic import BaseModel, Field


class ArtifactStatus(BaseModel):
    """Machine-readable status. Failures are first-class data, not only exceptions."""

    status: str = "success"  # success | failed | skipped | partial
    category: str | None = None
    reason: str | None = None
    recoverable: bool = True
    recommended_fallback: str | None = None


class Artifact(BaseModel):
    artifact_id: str
    artifact_type: str
    parents: list[str] = Field(default_factory=list)
    paths: dict[str, str] = Field(default_factory=dict)
    data: dict[str, Any] = Field(default_factory=dict)
    method: dict[str, Any] | None = None
    provenance: dict[str, Any] = Field(default_factory=dict)
    qc: dict[str, Any] = Field(default_factory=dict)
    status: ArtifactStatus = Field(default_factory=ArtifactStatus)

    @staticmethod
    def now_iso() -> str:
        return datetime.now(timezone.utc).isoformat()

    @classmethod
    def failure(
        cls,
        artifact_id: str,
        artifact_type: str,
        reason: str,
        category: str = "unknown",
        parents: list[str] | None = None,
        recoverable: bool = True,
        recommended_fallback: str | None = None,
        **data: Any,
    ) -> Artifact:
        payload = data.pop("data") if set(data.keys()) == {"data"} and isinstance(data.get("data"), dict) else data
        return cls(
            artifact_id=artifact_id,
            artifact_type=artifact_type,
            parents=parents or [],
            data=payload,
            status=ArtifactStatus(
                status="failed",
                category=category,
                reason=reason,
                recoverable=recoverable,
                recommended_fallback=recommended_fallback,
            ),
            provenance={"created_at": cls.now_iso()},
        )

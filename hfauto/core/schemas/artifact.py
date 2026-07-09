from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field


class ArtifactStatus(BaseModel):
    """Machine-readable status. Failures are first-class data, not only exceptions."""

    status: str = "success"  # success | failed | skipped | partial
    category: Optional[str] = None
    reason: Optional[str] = None
    recoverable: bool = True
    recommended_fallback: Optional[str] = None


class Artifact(BaseModel):
    artifact_id: str
    artifact_type: str
    parents: List[str] = Field(default_factory=list)
    paths: Dict[str, str] = Field(default_factory=dict)
    data: Dict[str, Any] = Field(default_factory=dict)
    method: Optional[Dict[str, Any]] = None
    provenance: Dict[str, Any] = Field(default_factory=dict)
    qc: Dict[str, Any] = Field(default_factory=dict)
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
        parents: Optional[list[str]] = None,
        recoverable: bool = True,
        recommended_fallback: Optional[str] = None,
        **data: Any,
    ) -> "Artifact":
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

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field


class SiteRecord(BaseModel):
    """Candidate HF interaction/protonation site on a molecule."""

    site_id: str
    mol_id: str
    atom_index: int
    site_type: str
    site_smarts: str
    priority: int = 1
    site_confidence: float = 0.0
    excluded: bool = False
    exclude_reason: str | None = None
    local_environment: dict[str, Any] = Field(default_factory=dict)
    basicity_proxy: dict[str, Any] = Field(default_factory=dict)

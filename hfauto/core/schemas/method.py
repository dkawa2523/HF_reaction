"""Normalized electronic-structure method records."""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, ConfigDict, Field


class ElectronicStructureMethodRecord(BaseModel):
    """The scientifically relevant method lineage shared by all QM backends."""

    model_config = ConfigDict(extra="allow")

    engine: str
    task: str
    backend: str | None = None
    method_id: str | None = None
    functional: str | None = None
    basis: str | None = None
    dispersion: str | int | None = None
    solvation_model: str | None = None
    dielectric: float | None = None
    required_program_version: str | None = None
    charge: int | None = None
    multiplicity: int | None = None
    settings: dict[str, Any] = Field(default_factory=dict)

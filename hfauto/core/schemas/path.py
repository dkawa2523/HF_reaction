"""Normalized records for reaction-path calculations."""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, Field, field_validator

PathClassification = Literal[
    "missing_profile",
    "resolved_internal_maximum",
    "monotonic_no_internal_maximum",
    "unresolved_internal_profile",
]
PathDiagnosis = Literal[
    "resolved_saddle_candidate",
    "multiple_step_candidate",
    "same_basin",
    "endpoint_path_inconsistency",
    "path_stagnation",
    "path_unresolved",
]


class PathImageRecord(BaseModel):
    """Energy evidence for one ordered image along a reaction path."""

    index: int
    energy_hartree: float | None = None
    comment: str | None = None

    @field_validator("index")
    @classmethod
    def validate_index(cls, value: int) -> int:
        if int(value) < 0:
            raise ValueError("path image index must be non-negative")
        return int(value)


class IntermediateWellRecord(BaseModel):
    """One path-local minimum bracketed by two resolved maxima."""

    image_index: int
    left_maximum_index: int
    right_maximum_index: int
    left_drop_kcal_mol: float
    right_drop_kcal_mol: float
    resolved: bool


class ReactionPathRecord(BaseModel):
    """Backend-neutral NEB/string result used by TS selection and reporting."""

    reaction_id: str
    engine: str
    converged: bool
    images: list[PathImageRecord] = Field(default_factory=list)
    classification: PathClassification
    highest_internal_image_index: int | None = None
    internal_maximum_indices: list[int] = Field(default_factory=list)
    intermediate_wells: list[IntermediateWellRecord] = Field(default_factory=list)
    ts_guess_image_index: int | None = None
    barrier_from_reactant_kcal_mol: float | None = None
    maximum_internal_rise_from_reactant_kcal_mol: float | None = None
    reaction_energy_kcal_mol: float | None = None
    barrier_threshold_kcal_mol: float
    energy_source: str
    unresolved_reasons: list[str] = Field(default_factory=list)


class PathOptimizationRecord(BaseModel):
    """Backend-neutral convergence evidence for one path attempt."""

    iterations: int | None = None
    max_gradient: float | None = None
    rms_gradient: float | None = None
    max_displacement: float | None = None
    rms_displacement: float | None = None
    max_gradient_history: list[float] = Field(default_factory=list)
    converged: bool = False
    stagnant: bool = False


class PathAttemptRecord(BaseModel):
    """One auditable path-search attempt and its method-neutral diagnosis."""

    attempt_id: str
    reaction_id: str
    strategy: str
    engine: str
    endpoint_basin_status: str
    path: ReactionPathRecord
    optimization: PathOptimizationRecord = Field(default_factory=PathOptimizationRecord)
    diagnosis: PathDiagnosis
    next_action: str
    reasons: list[str] = Field(default_factory=list)
    evidence: dict[str, Any] = Field(default_factory=dict)

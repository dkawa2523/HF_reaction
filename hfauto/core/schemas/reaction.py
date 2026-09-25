from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field

from hfauto.core.schemas.chemistry import (
    BondChangeRecord,
    ReactionCoordinateRecord,
    ReactionHypothesisRecord,
)


class ReactionRecord(BaseModel):
    """Reaction-level contract for an atom-mapped pair of chemical states."""

    reaction_id: str
    reaction_type: str
    mol_id: str | None = None
    site_id: str | None = None
    # Compatibility field for existing HF runs.
    hf_n: int = 0
    conformer_id: str | None = None
    candidate_species_id: str | None = None
    hf_cluster_species_id: str | None = None
    reactant_species_id: str
    product_species_id: str
    ts_species_id: str | None = None
    mechanism_family: str | None = None
    bond_changes: list[BondChangeRecord | dict[str, Any]] = Field(default_factory=list)
    reaction_coordinate: ReactionCoordinateRecord | dict[str, Any] = Field(
        default_factory=dict
    )
    stoichiometry: dict[str, list[dict[str, Any]]] = Field(default_factory=dict)
    hypothesis: ReactionHypothesisRecord | None = None
    environment: dict[str, Any] = Field(default_factory=dict)
    qc: dict[str, Any] = Field(default_factory=dict)

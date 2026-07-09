from __future__ import annotations

from typing import Any, Optional

from pydantic import BaseModel, Field


class ReactionRecord(BaseModel):
    """Reaction-level contract for a single HF process around one site and conformer."""

    reaction_id: str
    reaction_type: str
    mol_id: str
    site_id: str
    hf_n: int
    conformer_id: Optional[str] = None
    candidate_species_id: Optional[str] = None
    hf_cluster_species_id: Optional[str] = None
    reactant_species_id: str
    product_species_id: str
    ts_species_id: Optional[str] = None
    reaction_coordinate: dict[str, Any] = Field(default_factory=dict)
    qc: dict[str, Any] = Field(default_factory=dict)

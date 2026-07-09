from __future__ import annotations

from typing import Any, Literal, Optional

from pydantic import BaseModel, Field


SpeciesState = Literal[
    "bare_candidate",
    "hf_cluster",
    "reactant_complex",
    "ion_pair",
    "transition_state",
    "probe_reactant",
    "probe_product",
]


class SpeciesRecord(BaseModel):
    """A concrete molecular system with one XYZ geometry and a charge/spin state."""

    species_id: str
    mol_id: Optional[str] = None
    site_id: Optional[str] = None
    conformer_id: Optional[str] = None
    state: SpeciesState
    hf_n: int = 0
    charge: int = 0
    multiplicity: int = 1
    xyz_path: str
    atom_order_key: Optional[str] = None
    components: list[dict[str, Any]] = Field(default_factory=list)

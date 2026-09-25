from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field

from hfauto.core.schemas.chemistry import ChemicalStateRecord, ComponentRecord


class SpeciesRecord(BaseModel):
    """A concrete molecular system with one XYZ geometry and a charge/spin state."""

    species_id: str
    mol_id: str | None = None
    site_id: str | None = None
    conformer_id: str | None = None
    state: str
    # Compatibility field for existing HF runs. New chemistry is described by
    # ``chemical_state.components`` and must not branch on this value.
    hf_n: int = 0
    charge: int = 0
    multiplicity: int = 1
    xyz_path: str
    atom_order_key: str | None = None
    components: list[ComponentRecord | dict[str, Any]] = Field(default_factory=list)
    chemical_state: ChemicalStateRecord | None = None
    environment: dict[str, Any] = Field(default_factory=dict)

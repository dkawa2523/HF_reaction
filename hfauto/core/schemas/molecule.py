from __future__ import annotations

from typing import Any, Optional

from pydantic import BaseModel, Field


class MoleculeRecord(BaseModel):
    """Canonical molecule-level data extracted from the input SDF and public DBs."""

    mol_id: str
    source_sdf_index: int
    name: str
    canonical_smiles: Optional[str] = None
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    formal_charge: Optional[int] = 0
    multiplicity: int = 1
    num_atoms: Optional[int] = None
    has_3d: bool = False
    sdf_props: dict[str, Any] = Field(default_factory=dict)
    identity: dict[str, Any] = Field(default_factory=dict)
    public_data: dict[str, Any] = Field(default_factory=dict)
    ingest_qc: dict[str, Any] = Field(default_factory=dict)

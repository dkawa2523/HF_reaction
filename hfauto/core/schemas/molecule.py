from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field


class MoleculeRecord(BaseModel):
    """Canonical molecule-level data extracted from the input SDF and public DBs."""

    mol_id: str
    source_sdf_index: int
    name: str
    canonical_smiles: str | None = None
    isomeric_smiles: str | None = None
    inchi: str | None = None
    inchikey: str | None = None
    formula: str | None = None
    exact_mw: float | None = None
    formal_charge: int | None = 0
    multiplicity: int = 1
    num_atoms: int | None = None
    num_fragments: int = 1
    has_3d: bool = False
    sdf_props: dict[str, Any] = Field(default_factory=dict)
    identity: dict[str, Any] = Field(default_factory=dict)
    public_data: dict[str, Any] = Field(default_factory=dict)
    ingest_qc: dict[str, Any] = Field(default_factory=dict)
    extras: dict[str, Any] = Field(default_factory=dict)

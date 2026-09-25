from __future__ import annotations

from typing import Any

from pydantic import BaseModel, ConfigDict, Field


class ThermoRecord(BaseModel):
    """Method-neutral reaction thermochemistry with legacy fields as optional metadata."""

    model_config = ConfigDict(extra="allow")

    reaction_id: str
    T_K: float
    quality_tier: str
    stoichiometry: dict[str, list[dict[str, Any]]] = Field(default_factory=dict)
    delta_G_reaction_standard_kcal_mol: float | None = None
    delta_G_reaction_process_kcal_mol: float | None = None
    delta_G_activation_standard_kcal_mol: float | None = None
    delta_G_activation_process_kcal_mol: float | None = None
    delta_E_activation_kcal_mol: float | None = None
    delta_E_activation_unresolved_kcal_mol: float | None = None
    activation_energy_resolved: bool | None = None
    delta_E0_activation_kcal_mol: float | None = None
    delta_H_activation_standard_kcal_mol: float | None = None
    activation_connectivity_validated: bool = False
    # Compatibility metadata for existing HF screening runs.
    mol_id: str | None = None
    site_id: str | None = None
    hf_n: int = 0
    G_candidate_hartree: float | None = None
    G_hf_cluster_hartree: float | None = None
    G_reactant_complex_hartree: float | None = None
    G_product_ionpair_hartree: float | None = None
    G_TS_hartree: float | None = None
    delta_G_assoc_kcal_mol: float | None = None
    delta_G_ionpair_kcal_mol: float | None = None
    delta_G_act_kcal_mol: float | None = None

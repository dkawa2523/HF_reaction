from __future__ import annotations

from typing import Optional

from pydantic import BaseModel


class ThermoRecord(BaseModel):
    reaction_id: str
    mol_id: str
    site_id: str
    hf_n: int
    T_K: float
    quality_tier: str
    G_candidate_hartree: Optional[float] = None
    G_hf_cluster_hartree: Optional[float] = None
    G_reactant_complex_hartree: Optional[float] = None
    G_product_ionpair_hartree: Optional[float] = None
    G_TS_hartree: Optional[float] = None
    delta_G_assoc_kcal_mol: Optional[float] = None
    delta_G_ionpair_kcal_mol: Optional[float] = None
    delta_G_act_kcal_mol: Optional[float] = None

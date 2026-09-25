from __future__ import annotations

from pydantic import BaseModel


class DescriptorRecord(BaseModel):
    species_id: str
    mol_id: str
    site_id: str
    hf_n: int
    r_HF_A: float | None = None
    delta_r_HF_A: float | None = None
    nu_HF_cm1: float | None = None
    delta_nu_HF_cm1: float | None = None
    B_H_distance_A: float | None = None
    B_H_F_angle_deg: float | None = None
    min_nonbonded_distance_A: float | None = None

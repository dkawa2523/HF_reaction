from __future__ import annotations

from typing import Optional

from pydantic import BaseModel


class DescriptorRecord(BaseModel):
    species_id: str
    mol_id: str
    site_id: str
    hf_n: int
    r_HF_A: Optional[float] = None
    delta_r_HF_A: Optional[float] = None
    nu_HF_cm1: Optional[float] = None
    delta_nu_HF_cm1: Optional[float] = None
    B_H_distance_A: Optional[float] = None
    B_H_F_angle_deg: Optional[float] = None
    min_nonbonded_distance_A: Optional[float] = None

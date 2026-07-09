from __future__ import annotations

from typing import Any, Optional

from pydantic import BaseModel, Field


class CalculationRecord(BaseModel):
    """Normalized quantum-calculation result payload stored inside calculation artifacts."""

    calc_id: str
    species_id: str
    task: str
    engine: str
    method_id: str
    status: str = "success"
    electronic_energy_hartree: Optional[float] = None
    zpe_hartree: Optional[float] = None
    enthalpy_298K_hartree: Optional[float] = None
    gibbs_298K_hartree: Optional[float] = None
    n_imag: Optional[int] = None
    imag_freq_cm1: Optional[float] = None
    hf_stretch_cm1: Optional[float] = None
    qc: dict[str, Any] = Field(default_factory=dict)

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field


class CalculationRecord(BaseModel):
    """Normalized quantum-calculation result payload stored inside calculation artifacts."""

    calc_id: str
    species_id: str
    task: str
    engine: str
    method_id: str
    status: str = "success"
    electronic_energy_hartree: float | None = None
    zpe_hartree: float | None = None
    enthalpy_298K_hartree: float | None = None
    gibbs_298K_hartree: float | None = None
    n_imag: int | None = None
    imag_freq_cm1: float | None = None
    hf_stretch_cm1: float | None = None
    qc: dict[str, Any] = Field(default_factory=dict)

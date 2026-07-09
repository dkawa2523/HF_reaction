from __future__ import annotations

import math

from hfauto.core.constants import HARTREE_TO_KCAL_MOL, HARTREE_TO_KJ_MOL, R_KCAL_MOL_K


def hartree_to_kcal_mol(value: float) -> float:
    return float(value) * HARTREE_TO_KCAL_MOL


def hartree_to_kj_mol(value: float) -> float:
    return float(value) * HARTREE_TO_KJ_MOL


def kcal_mol_to_hartree(value: float) -> float:
    return float(value) / HARTREE_TO_KCAL_MOL


def pressure_correction_hartree(T_K: float, p_bar: float | None, p_standard_bar: float = 1.0) -> float:
    """Ideal-gas chemical-potential correction RT ln(p/p°), in Hartree."""
    if p_bar is None:
        return 0.0
    p = max(float(p_bar), 1e-300)
    p0 = max(float(p_standard_bar), 1e-300)
    return kcal_mol_to_hartree(R_KCAL_MOL_K * float(T_K) * math.log(p / p0))

from __future__ import annotations

"""Small, auditable thermochemistry/kinetics helpers.

The functions here are intentionally conservative.  They do not try to replace
full thermochemistry packages such as GoodVibes or Arkane; instead they provide
an internal fallback model with explicit provenance so the workflow remains
reviewable when external tools are unavailable.
"""

from collections.abc import Iterable
from dataclasses import asdict, dataclass
from math import exp, log
from typing import Any

from hfauto.core.constants import KB_OVER_H_PER_K_S, R_KCAL_MOL_K
from hfauto.core.units import kcal_mol_to_hartree

CM1_TO_K = 1.438776877  # h*c/kB in K cm


@dataclass
class SpeciesThermoResult:
    species_id: str | None
    T_K: float
    p_bar: float | None
    standard_pressure_bar: float
    electronic_energy_hartree: float | None
    thermal_correction_gibbs_hartree: float
    G_standard_hartree: float | None
    pressure_correction_hartree: float
    G_final_hartree: float | None
    thermal_source_artifact_id: str | None
    electronic_source_artifact_id: str | None
    thermo_backend: str
    thermo_model: str
    quasi_rrho_applied: bool
    quasi_rrho_delta_hartree: float
    low_frequency_count: int
    frequency_scale_factor: float
    frequencies_used_count: int
    scientific_thermo: bool

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def positive_frequencies(freqs: Iterable[Any] | None, scale: float = 1.0) -> list[float]:
    out: list[float] = []
    for f in freqs or []:
        try:
            x = float(f) * float(scale)
        except Exception:
            continue
        if x > 0.0:
            out.append(x)
    return out


def low_frequency_summary(freqs: Iterable[Any] | None, cutoff_cm1: float = 100.0, scale: float = 1.0) -> tuple[int, float]:
    """Return count and a simple qRRHO free-energy penalty in Hartree.

    Low-frequency modes in weak HF complexes inflate RRHO entropy.  The fallback
    model adds a small positive ΔG penalty for modes below ``cutoff_cm1``.  This
    is a proxy model, not a GoodVibes replacement; the output explicitly records
    ``thermo_model=internal_quasirrho_proxy``.
    """
    pos = positive_frequencies(freqs, scale=scale)
    weights = [max(0.0, (float(cutoff_cm1) - f) / float(cutoff_cm1)) for f in pos if f < float(cutoff_cm1)]
    # 0.06 kcal/mol per fully floppy low mode at 298 K; scaled by weight.
    penalty_kcal = 0.06 * sum(weights)
    return len(weights), kcal_mol_to_hartree(penalty_kcal)


def correction_from_calc(calc_data: dict[str, Any] | None) -> float:
    calc_data = calc_data or {}
    if calc_data.get("thermal_correction_gibbs_hartree") is not None:
        return float(calc_data["thermal_correction_gibbs_hartree"])
    if calc_data.get("gibbs_298K_hartree") is not None and calc_data.get("electronic_energy_hartree") is not None:
        return float(calc_data["gibbs_298K_hartree"]) - float(calc_data["electronic_energy_hartree"])
    return 0.0


def temperature_adjustment_proxy(T_K: float, base_T_K: float = 298.15, low_frequency_count: int = 0) -> float:
    """Very small, explicit proxy for T dependence when only 298 K data exist.

    Full T-dependent thermochemistry requires partition functions; ORCA/GoodVibes
    should be used in production.  This fallback intentionally applies only a
    transparent low-frequency entropy proxy so multi-temperature tables remain
    differentiable without pretending to be high-accuracy thermochemistry.
    """
    dT = float(T_K) - float(base_T_K)
    if abs(dT) < 1e-12 or low_frequency_count <= 0:
        return 0.0
    # Penalize high-T entropy from floppy modes modestly: +0.00003 kcal/mol/K/mode.
    return kcal_mol_to_hartree(0.00003 * dT * float(low_frequency_count))


def equilibrium_constant_from_delta_g(delta_g_kcal_mol: float | None, T_K: float) -> float | None:
    if delta_g_kcal_mol is None:
        return None
    exponent = -float(delta_g_kcal_mol) / (R_KCAL_MOL_K * float(T_K))
    # Avoid overflow in pathological dummy/fallback cases.
    exponent = max(-700.0, min(700.0, exponent))
    return exp(exponent)


def eyring_rate(delta_g_act_kcal_mol: float, T_K: float) -> float:
    exponent = -float(delta_g_act_kcal_mol) / (R_KCAL_MOL_K * float(T_K))
    exponent = max(-700.0, min(700.0, exponent))
    return KB_OVER_H_PER_K_S * float(T_K) * exp(exponent)


def wigner_tunneling_factor(imag_freq_cm1: float | None, T_K: float) -> float:
    """Wigner tunneling correction for one imaginary frequency.

    kappa = 1 + (1/24) * (h c |ν‡| / kB T)^2
    """
    if imag_freq_cm1 is None:
        return 1.0
    nu = abs(float(imag_freq_cm1))
    if nu <= 0.0 or T_K <= 0.0:
        return 1.0
    x = CM1_TO_K * nu / float(T_K)
    return 1.0 + (x * x) / 24.0


def arrhenius_fit_from_rates(records: list[dict[str, Any]]) -> dict[str, float | None]:
    """Fit ln(k) = ln(A) - Ea/R * 1/T for records sharing one reaction.

    Returns a minimal fit dictionary.  With fewer than two positive rates, the
    fit is undefined and values are ``None``.
    """
    pts: list[tuple[float, float]] = []
    for r in records:
        try:
            T = float(r["T_K"])
            k = float(r.get("k_corrected_s-1") or r.get("k_TST_s-1"))
        except Exception:
            continue
        if T > 0 and k > 0:
            pts.append((1.0 / T, log(k)))
    if len(pts) < 2:
        return {"arrhenius_A_s-1": None, "arrhenius_Ea_kcal_mol": None, "arrhenius_fit_points": len(pts)}
    # simple least squares without numpy dependency here
    n = float(len(pts))
    sx = sum(p[0] for p in pts)
    sy = sum(p[1] for p in pts)
    sxx = sum(p[0] * p[0] for p in pts)
    sxy = sum(p[0] * p[1] for p in pts)
    denom = n * sxx - sx * sx
    if abs(denom) < 1e-30:
        return {"arrhenius_A_s-1": None, "arrhenius_Ea_kcal_mol": None, "arrhenius_fit_points": len(pts)}
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n
    Ea = -slope * R_KCAL_MOL_K
    A = exp(max(-700.0, min(700.0, intercept)))
    return {"arrhenius_A_s-1": A, "arrhenius_Ea_kcal_mol": Ea, "arrhenius_fit_points": len(pts)}

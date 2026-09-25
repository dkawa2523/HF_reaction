from __future__ import annotations

"""Auditable thermochemistry and kinetics helpers for hfauto.

The functions in this module provide deterministic fallback behaviour for CI,
code review and lightweight screening.  They are intentionally transparent and
metadata-rich rather than a replacement for GoodVibes/Arkane/Cantera production
calculations.
"""

import math
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any

from hfauto.core.constants import KB_OVER_H_PER_K_S, R_KCAL_MOL_K
from hfauto.core.units import kcal_mol_to_hartree

HC_OVER_KB_CM_K = 1.438776877


def finite_float(value: Any, default: float | None = None) -> float | None:
    try:
        x = float(value)
    except Exception:
        return default
    if math.isnan(x) or math.isinf(x):
        return default
    return x


def positive_frequencies(frequencies: Iterable[Any] | None, scale_factor: float = 1.0) -> list[float]:
    out: list[float] = []
    for f in frequencies or []:
        x = finite_float(f)
        if x is not None and x > 0.0:
            out.append(x * float(scale_factor))
    return out


def low_frequency_count(frequencies_cm1: Iterable[Any] | None, cutoff_cm1: float = 100.0) -> int:
    cutoff = max(float(cutoff_cm1), 1.0)
    return sum(1 for f in positive_frequencies(frequencies_cm1) if f < cutoff)


def low_frequency_summary(frequencies_cm1: Iterable[Any] | None, cutoff_cm1: float = 100.0) -> dict[str, Any]:
    cutoff = max(float(cutoff_cm1), 1.0)
    freqs = positive_frequencies(frequencies_cm1)
    low = [f for f in freqs if f < cutoff]
    weights = [max(0.0, (cutoff - f) / cutoff) for f in low]
    return {
        "frequency_count": len(freqs),
        "low_frequency_cutoff_cm1": cutoff,
        "low_frequency_count": len(low),
        "low_frequency_weighted_count": float(sum(weights)),
        "lowest_positive_frequency_cm1": min(freqs) if freqs else None,
    }


@dataclass(frozen=True)
class QuasiRRHOResult:
    frequency_scale_factor: float
    low_freq_cutoff_cm1: float
    low_frequency_count: int
    positive_frequency_count: int
    quasi_rrho_delta_G_kcal_mol: float
    quasi_rrho_delta_G_hartree: float
    model: str


def quasi_rrho_low_frequency_correction(
    freqs_cm1: Iterable[Any] | None,
    T_K: float = 298.15,
    scale_factor: float = 1.0,
    cutoff_cm1: float = 100.0,
    alpha: float = 0.50,
    enabled: bool = True,
) -> QuasiRRHOResult:
    """Conservative low-frequency correction to G.

    For each positive frequency below the cutoff, add
    ``alpha * R*T*ln(cutoff/nu)`` to G.  This moves floppy weak complexes in the
    same direction as quasi-RRHO entropy treatments, while clearly labelling the
    model as an internal fallback rather than full GoodVibes.
    """
    pos = positive_frequencies(freqs_cm1, scale_factor=scale_factor)
    if not enabled:
        return QuasiRRHOResult(float(scale_factor), float(cutoff_cm1), 0, len(pos), 0.0, 0.0, "disabled")
    cutoff = max(float(cutoff_cm1), 1e-9)
    low = [max(f, 1e-9) for f in pos if f < cutoff]
    delta = sum(float(alpha) * R_KCAL_MOL_K * float(T_K) * math.log(cutoff / f) for f in low)
    return QuasiRRHOResult(
        frequency_scale_factor=float(scale_factor),
        low_freq_cutoff_cm1=cutoff,
        low_frequency_count=len(low),
        positive_frequency_count=len(pos),
        quasi_rrho_delta_G_kcal_mol=float(delta),
        quasi_rrho_delta_G_hartree=kcal_mol_to_hartree(delta),
        model="hfauto_clamped_quasi_rrho",
    )


def quasi_rrho_proxy_correction_hartree(
    frequencies_cm1: Iterable[Any] | None,
    T_K: float,
    cutoff_cm1: float = 100.0,
    penalty_kcal_mol_per_weighted_mode_298K: float = 0.10,
) -> tuple[float, dict[str, Any]]:
    summary = low_frequency_summary(frequencies_cm1, cutoff_cm1=cutoff_cm1)
    weighted = float(summary["low_frequency_weighted_count"])
    penalty = float(penalty_kcal_mol_per_weighted_mode_298K) * weighted * (float(T_K) / 298.15)
    return kcal_mol_to_hartree(penalty), {
        **summary,
        "quasi_rrho_model": "hfauto_low_frequency_proxy",
        "quasi_rrho_warning": "Screening proxy only; replace with GoodVibes/ORCA quasi-RRHO for production thermochemistry.",
        "quasi_rrho_correction_kcal_mol": penalty,
        "quasi_rrho_correction_hartree": kcal_mol_to_hartree(penalty),
    }


def zpe_scale_correction_hartree(zpe_hartree: float | None, frequency_scale_factor: float = 1.0) -> float:
    zpe = finite_float(zpe_hartree, 0.0) or 0.0
    return zpe * ((finite_float(frequency_scale_factor, 1.0) or 1.0) - 1.0)


def vibrational_temperature_delta_hartree(
    frequencies_cm1: Iterable[Any] | None,
    T_K: float,
    reference_T_K: float = 298.15,
    scale_factor: float = 1.0,
    quasi_rrho: bool = False,
    cutoff_cm1: float = 100.0,
) -> float:
    """Small transparent temperature proxy for optional offline screening.

    Full temperature-dependent thermochemistry requires complete partition
    functions.  The fallback only adjusts low-frequency contributions mildly, so
    multi-temperature tables are not identical while remaining conservative.
    """
    if abs(float(T_K) - float(reference_T_K)) < 1e-12:
        return 0.0
    count = low_frequency_count(positive_frequencies(frequencies_cm1, scale_factor), cutoff_cm1) if quasi_rrho else 0
    return kcal_mol_to_hartree(0.00003 * (float(T_K) - float(reference_T_K)) * float(count))


def thermal_correction_from_calc(calc_data: dict[str, Any], T_K: float, settings: dict[str, Any] | None = None) -> dict[str, Any]:
    settings = settings or {}
    e = calc_data.get("electronic_energy_hartree")
    g = calc_data.get("gibbs_298K_hartree")
    if e is not None and g is not None:
        base_corr = float(g) - float(e)
        source = "gibbs_minus_electronic"
    elif calc_data.get("thermal_correction_gibbs_hartree") is not None:
        base_corr = float(calc_data["thermal_correction_gibbs_hartree"])
        source = "thermal_correction_gibbs_hartree"
    else:
        base_corr = 0.0
        source = "missing_gibbs_default_zero"
    scale = float(settings.get("frequency_scale_factor", 1.0))
    cutoff = float(settings.get("low_freq_cutoff_cm1", settings.get("quasi_rrho_cutoff_cm1", 100.0)))
    alpha = float(settings.get("quasi_rrho_alpha", 0.50))
    apply = bool(settings.get("apply_quasi_rrho", settings.get("quasi_rrho", True)))
    q = quasi_rrho_low_frequency_correction(calc_data.get("frequencies_cm1", []), T_K=T_K, scale_factor=scale, cutoff_cm1=cutoff, alpha=alpha, enabled=apply)
    return {
        "T_K": float(T_K),
        "thermal_correction_gibbs_hartree": base_corr + q.quasi_rrho_delta_G_hartree,
        "base_thermal_correction_gibbs_hartree": base_corr,
        "thermal_correction_source": source,
        "quasi_rrho_applied": apply,
        "quasi_rrho_model": q.model,
        "frequency_scale_factor": q.frequency_scale_factor,
        "low_freq_cutoff_cm1": q.low_freq_cutoff_cm1,
        "low_frequency_count": q.low_frequency_count,
        "positive_frequency_count": q.positive_frequency_count,
        "quasi_rrho_delta_G_kcal_mol": q.quasi_rrho_delta_G_kcal_mol,
        "quasi_rrho_delta_G_hartree": q.quasi_rrho_delta_G_hartree,
    }


def eyring_prefactor_s(T_K: float) -> float:
    return KB_OVER_H_PER_K_S * float(T_K)


def tst_rate_s(T_K: float, delta_G_act_kcal_mol: float, kappa: float = 1.0) -> float:
    exponent = -float(delta_G_act_kcal_mol) / (R_KCAL_MOL_K * float(T_K))
    exponent = max(min(exponent, 700.0), -700.0)
    return float(kappa) * eyring_prefactor_s(T_K) * math.exp(exponent)


def tst_rate(delta_G_act_kcal_mol: float, T_K: float, transmission_coefficient: float = 1.0) -> float:
    return tst_rate_s(T_K, delta_G_act_kcal_mol, transmission_coefficient)


def eyring_rate_s_inv(delta_G_act_kcal_mol: float, T_K: float, kappa: float = 1.0) -> float:
    return tst_rate_s(T_K, delta_G_act_kcal_mol, kappa)


def wigner_tunneling_factor(arg1: float | None, arg2: float | None) -> float:
    """Wigner factor accepting either (T, imag_freq) or (imag_freq, T).

    Older phase code used both conventions.  Values above 2000 are usually
    frequencies; values between 50 and 1500 may be either, so keyword-free calls
    from current stages should pass (imag_freq, T) where T is positive.
    """
    if arg1 is None or arg2 is None:
        return 1.0
    a = float(arg1)
    b = float(arg2)
    # Heuristic: temperature is typically 200-1000 K; imaginary frequency often
    # negative or hundreds-thousands cm^-1.  Negative value is definitely freq.
    if a < 0.0:
        imag, T = a, b
    elif b < 0.0:
        imag, T = b, a
    elif a > 1500.0 and b < 1500.0:
        imag, T = a, b
    elif b > 1500.0 and a < 1500.0:
        imag, T = b, a
    else:
        # Default to current SimpleKinetics convention: (imag, T).
        imag, T = a, b
    nu = abs(float(imag))
    if nu <= 0.0 or T <= 0.0:
        return 1.0
    x = HC_OVER_KB_CM_K * nu / T
    return 1.0 + (x * x) / 24.0


def equilibrium_constant_from_delta_g(*args) -> float | None:
    """Return exp(-ΔG/RT), accepting (deltaG,T) or (T,deltaG)."""
    if len(args) != 2:
        raise TypeError("equilibrium_constant_from_delta_g expects two arguments")
    a, b = args
    # Temperatures are positive and usually > 100.  ΔG can be negative/small.
    if a is not None and float(a) > 100.0 and (b is None or abs(float(b)) < 100.0):
        T, dg = float(a), b
    else:
        dg, T = a, float(b)
    if dg is None:
        return None
    exponent = -float(dg) / (R_KCAL_MOL_K * float(T))
    exponent = max(min(exponent, 700.0), -700.0)
    return math.exp(exponent)


def arrhenius_fit(records: list[dict[str, Any]], k_field: str = "k_corrected_s-1") -> dict[str, Any]:
    xs: list[float] = []
    ys: list[float] = []
    for r in records:
        T = finite_float(r.get("T_K"))
        k = finite_float(r.get(k_field))
        if T and k and k > 0.0:
            xs.append(1.0 / T)
            ys.append(math.log(k))
    if len(xs) < 2:
        return {"fit_status": "insufficient_points", "A_s-1": None, "Ea_kcal_mol": None, "arrhenius_b": 0.0, "r2": None, "n_points": len(xs)}
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    denom = sum((x - mx) ** 2 for x in xs)
    if denom < 1e-30:
        return {"fit_status": "singular", "A_s-1": None, "Ea_kcal_mol": None, "arrhenius_b": 0.0, "r2": None, "n_points": n}
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / denom
    intercept = my - slope * mx
    Ea = -slope * R_KCAL_MOL_K
    A = math.exp(max(min(intercept, 700.0), -700.0))
    ss_tot = sum((y - my) ** 2 for y in ys)
    ss_res = sum((y - (intercept + slope * x)) ** 2 for x, y in zip(xs, ys))
    r2 = None if ss_tot < 1e-30 else 1.0 - ss_res / ss_tot
    return {"fit_status": "success", "A_s-1": A, "Ea_kcal_mol": Ea, "arrhenius_b": 0.0, "r2": r2, "n_points": n}

def equilibrium_constant(delta_G_kcal_mol: float | None, T_K: float) -> float | None:
    return equilibrium_constant_from_delta_g(delta_G_kcal_mol, T_K)

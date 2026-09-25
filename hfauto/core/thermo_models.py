from __future__ import annotations

"""Small, auditable thermochemistry helpers.

These helpers do not try to replace full-featured thermochemistry packages such
as GoodVibes or Arkane.  They provide deterministic in-package behavior so the
HF reaction pipeline can be run, reviewed, and regression-tested without those
optional tools.  External backends can replace these calculations while keeping
record shapes stable.
"""

import math
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from typing import Any

from hfauto.core.constants import HARTREE_TO_KCAL_MOL, KB_OVER_H_PER_K_S, R_KCAL_MOL_K
from hfauto.core.units import kcal_mol_to_hartree, pressure_correction_hartree

CM_TO_K = 1.438776877  # h*c/kB in K cm


@dataclass(frozen=True)
class LowFrequencyCorrection:
    """Transparent low-frequency quasi-RRHO style correction summary.

    The correction is intentionally conservative and explicitly marked as an
    internal model.  It increases G for floppy low-frequency modes, mimicking the
    direction of GoodVibes/Grimme/Truhlar corrections without claiming numerical
    equivalence.  Production workflows should prefer GoodVibes/ORCA quasi-RRHO.
    """

    model: str
    cutoff_cm1: float
    low_frequency_count: int
    correction_kcal_mol: float
    correction_hartree: float
    applied: bool

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def positive_frequencies(freqs: Iterable[float] | None) -> list[float]:
    return [float(f) for f in (freqs or []) if float(f) > 0.0]


def low_frequency_correction(
    freqs_cm1: Iterable[float] | None,
    *,
    cutoff_cm1: float = 100.0,
    max_correction_per_mode_kcal_mol: float = 0.20,
    enabled: bool = True,
) -> LowFrequencyCorrection:
    """Return a simple low-frequency G correction.

    For each positive frequency below the cutoff, add a fractional penalty:
    ``max_per_mode * (1 - nu/cutoff)^2``.  The output records the model name and
    magnitude so downstream rankings can distinguish this internal estimate from
    external GoodVibes results.
    """

    freqs = positive_frequencies(freqs_cm1)
    if not enabled:
        return LowFrequencyCorrection("none", float(cutoff_cm1), 0, 0.0, 0.0, False)
    low = [f for f in freqs if f < float(cutoff_cm1)]
    corr = 0.0
    for f in low:
        corr += float(max_correction_per_mode_kcal_mol) * (1.0 - f / float(cutoff_cm1)) ** 2
    return LowFrequencyCorrection(
        model="internal_low_frequency_penalty_v1",
        cutoff_cm1=float(cutoff_cm1),
        low_frequency_count=len(low),
        correction_kcal_mol=float(corr),
        correction_hartree=float(kcal_mol_to_hartree(corr)),
        applied=bool(enabled and low),
    )


def scale_thermal_correction(
    thermal_correction_hartree: float,
    *,
    T_K: float,
    source_T_K: float = 298.15,
    enabled: bool = False,
) -> tuple[float, str]:
    """Optionally apply a simple T scaling to a thermal Gibbs correction.

    Full thermal re-evaluation needs rotational constants, frequencies and
    symmetry.  Since most QM backends only expose a 298 K Gibbs correction in
    this project skeleton, the default is to reuse the parsed correction and
    record that fact.  A linear scaling option is provided for exploratory runs,
    but is marked explicitly.
    """

    if not enabled or abs(float(T_K) - float(source_T_K)) < 1e-9:
        return float(thermal_correction_hartree), "source_temperature_correction_reused"
    return float(thermal_correction_hartree) * float(T_K) / float(source_T_K), "linear_thermal_correction_scaling_experimental"


def wigner_tunneling_factor(imag_freq_cm1: float | None, T_K: float) -> float:
    """Wigner tunneling correction for a TS imaginary frequency."""

    if imag_freq_cm1 is None:
        return 1.0
    theta = CM_TO_K * abs(float(imag_freq_cm1)) / max(float(T_K), 1e-12)
    return float(1.0 + (theta * theta) / 24.0)


def eyring_rate_s(delta_g_act_kcal_mol: float, T_K: float, transmission_coeff: float = 1.0) -> float:
    """Eyring transition-state-theory rate constant in s^-1."""

    return float(transmission_coeff) * KB_OVER_H_PER_K_S * float(T_K) * math.exp(-float(delta_g_act_kcal_mol) / (R_KCAL_MOL_K * float(T_K)))


def equilibrium_constant(delta_g_kcal_mol: float | None, T_K: float) -> float | None:
    if delta_g_kcal_mol is None:
        return None
    try:
        return float(math.exp(-float(delta_g_kcal_mol) / (R_KCAL_MOL_K * float(T_K))))
    except OverflowError:
        return math.inf if float(delta_g_kcal_mol) < 0 else 0.0


def process_adjusted_reaction_delta_g(
    delta_g_standard_kcal_mol: float | None,
    *,
    T_K: float,
    reactant_pressures_bar: list[float | None],
    product_pressures_bar: list[float | None] | None = None,
    p_standard_bar: float = 1.0,
) -> float | None:
    """Return ΔG under ideal-gas pressure assumptions.

    ΔG = ΔG° + RT ln(Q).  For screening we often know reactant partial
    pressures but not the product complex pressure.  Missing pressures are
    omitted and the output status in the caller should record the assumption.
    """

    if delta_g_standard_kcal_mol is None:
        return None
    dg_h = kcal_mol_to_hartree(float(delta_g_standard_kcal_mol))
    for p in product_pressures_bar or []:
        if p is not None:
            dg_h += pressure_correction_hartree(T_K, p, p_standard_bar)
    for p in reactant_pressures_bar:
        if p is not None:
            dg_h -= pressure_correction_hartree(T_K, p, p_standard_bar)
    return float(dg_h * HARTREE_TO_KCAL_MOL)

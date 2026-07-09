from __future__ import annotations

from typing import Any

from hfauto.core.thermo import low_frequency_count, vibrational_temperature_delta_hartree


class SimpleThermoEngine:
    """In-package thermal correction backend.

    This backend uses QM-code Gibbs corrections when present and adds an optional
    frequency-only temperature adjustment.  It records enough metadata to make
    the approximation explicit.  GoodVibes/Arkane adapters can return the same
    keys with higher-fidelity values.
    """

    name = "simple"

    def __init__(self, **kwargs: Any):
        self.config = kwargs

    def correction_hartree(self, thermal_calc, T_K: float, config: dict[str, Any]) -> tuple[float, dict[str, Any]]:
        data = thermal_calc.data if thermal_calc is not None else {}
        electronic = data.get("electronic_energy_hartree")
        if data.get("thermal_correction_gibbs_hartree") is not None:
            base_corr = float(data["thermal_correction_gibbs_hartree"])
            source = "thermal_correction_gibbs_hartree"
        elif data.get("gibbs_298K_hartree") is not None and electronic is not None:
            base_corr = float(data["gibbs_298K_hartree"]) - float(electronic)
            source = "gibbs_minus_electronic"
        else:
            base_corr = 0.0
            source = "missing_thermal_correction_zero_used"
        freqs = data.get("frequencies_cm1") or []
        quasi = bool(config.get("quasi_rrho", False))
        cutoff = float(config.get("quasi_rrho_cutoff_cm1", 100.0))
        scale = float(config.get("frequency_scale_factor", 1.0))
        ref_T = float(config.get("reference_temperature_K", 298.15))
        apply_temp = bool(config.get("temperature_adjustment", True))
        delta = 0.0
        if apply_temp and freqs:
            delta = vibrational_temperature_delta_hartree(
                freqs,
                T_K,
                ref_T,
                scale_factor=scale,
                quasi_rrho=quasi,
                cutoff_cm1=cutoff,
            )
        meta = {
            "thermal_backend": self.name,
            "thermal_correction_source": source,
            "reference_temperature_K": ref_T,
            "temperature_adjustment_model": "frequency_only_delta_Gvib" if apply_temp and freqs else "none",
            "frequency_scale_factor": scale,
            "quasi_rrho_applied": quasi,
            "quasi_rrho_cutoff_cm1": cutoff,
            "low_frequency_count": low_frequency_count(freqs, cutoff),
            "thermal_temperature_delta_hartree": delta,
        }
        return base_corr + delta, meta

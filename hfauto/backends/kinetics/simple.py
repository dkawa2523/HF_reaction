from __future__ import annotations

from typing import Any

from hfauto.core.thermo import eyring_rate_s_inv, wigner_tunneling_factor


class SimpleKineticsEngine:
    name = "simple"

    def __init__(self, **kwargs: Any):
        self.config = kwargs

    def rate_record(self, thermo_record: dict[str, Any], config: dict[str, Any]) -> dict[str, Any] | None:
        dg = thermo_record.get("delta_G_act_kcal_mol")
        T = thermo_record.get("T_K")
        if dg is None or T is None:
            return None
        k_tst = eyring_rate_s_inv(float(dg), float(T))
        tunneling = str(config.get("tunneling_model", "none")).lower()
        factor = 1.0
        if tunneling == "wigner":
            factor = wigner_tunneling_factor(thermo_record.get("imag_freq_cm1"), float(T))
        elif tunneling in {"none", "off", "false"}:
            factor = 1.0
        else:
            # Unknown correction models are kept explicit and neutral.
            factor = 1.0
        k_corr = k_tst * factor
        return {
            "reaction_id": thermo_record.get("reaction_id"),
            "mol_id": thermo_record.get("mol_id"),
            "site_id": thermo_record.get("site_id"),
            "hf_n": thermo_record.get("hf_n"),
            "T_K": float(T),
            "delta_G_act_kcal_mol": float(dg),
            "imag_freq_cm1": thermo_record.get("imag_freq_cm1"),
            "k_TST_s-1": k_tst,
            "tunneling_model": tunneling,
            "tunneling_factor": factor,
            "k_corrected_s-1": k_corr,
            "quality_tier": thermo_record.get("quality_tier"),
            "confidence_score": thermo_record.get("confidence_score"),
            "kinetics_quality": "TST_from_validated_thermo" if thermo_record.get("quality_tier") in {"Q4", "Q5"} else "TST_from_unvalidated_or_proxy_TS",
        }

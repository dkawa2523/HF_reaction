from __future__ import annotations

from typing import Any

from hfauto.core.thermo_models import equilibrium_constant, eyring_rate_s, wigner_tunneling_factor


class TSTKineticsEngine:
    name = "tst"

    def __init__(self, **config: Any):
        self.config = config

    def rate_record(self, thermo: dict[str, Any]) -> dict[str, Any] | None:
        dg = thermo.get("delta_G_act_kcal_mol")
        T = thermo.get("T_K")
        if dg is None or T is None:
            return None
        model = str(self.config.get("tunneling_model", "none")).lower()
        kappa = 1.0
        if model == "wigner":
            kappa = wigner_tunneling_factor(thermo.get("imag_freq_cm1"), float(T))
        elif model in {"none", "false", "0"}:
            kappa = 1.0
        else:
            # Placeholder for future Eckart/instanton interfaces; keep explicit.
            kappa = 1.0
            model = f"unsupported_{model}_treated_as_none"
        k_tst = eyring_rate_s(float(dg), float(T), 1.0)
        k_corr = eyring_rate_s(float(dg), float(T), kappa)
        record = {
            "reaction_id": thermo.get("reaction_id"),
            "candidate_id": thermo.get("candidate_id"),
            "comparison_group": thermo.get("comparison_group"),
            "mechanism_family": thermo.get("mechanism_family"),
            "T_K": float(T),
            "delta_G_act_kcal_mol": float(dg),
            "delta_G_reaction_kcal_mol": thermo.get("delta_G_reaction_kcal_mol"),
            "imag_freq_cm1": thermo.get("imag_freq_cm1"),
            "k_TST_s-1": float(k_tst),
            "transmission_coefficient": float(kappa),
            "tunneling_model": model,
            "k_corrected_s-1": float(k_corr),
            "K_assoc_standard": equilibrium_constant(thermo.get("delta_G_assoc_standard_kcal_mol"), float(T)),
            "K_assoc_process_adjusted": equilibrium_constant(thermo.get("delta_G_assoc_process_kcal_mol"), float(T)),
            "K_ionpair": equilibrium_constant(thermo.get("delta_G_ionpair_kcal_mol"), float(T)),
            "K_reaction_standard": thermo.get("K_reaction_standard"),
            "kinetics_backend": self.name,
            "kinetics_quality": "TST_from_validated_thermo" if thermo.get("quality_tier") in {"Q4", "Q5"} else "TST_from_unvalidated_or_proxy_TS",
        }
        for key in ("mol_id", "site_id", "site_type", "hf_n"):
            if thermo.get(key) is not None:
                record[key] = thermo[key]
        return record

from __future__ import annotations

from typing import Any

import pandas as pd

from hfauto.backends.registry import get_thermo_engine
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.qc import cap_confidence, confidence_score, production_rank_eligible, quality_tier, scientific_rank_eligible
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.core.thermo_models import equilibrium_constant
from hfauto.core.units import HARTREE_TO_KCAL_MOL
from hfauto.stages.base import Stage, StageContext


class ThermoStage(Stage):
    """Build species and reaction thermochemistry records.

    Phase 6 separates thermochemistry into two levels:
    - species_thermo: G° and pressure-adjusted G for each species/T;
    - thermo: reaction-level ΔG values, K values, quality and source IDs.
    """

    name = "thermo"

    @staticmethod
    def _is_scientific_calc(calc: Artifact) -> bool:
        if calc.status.status != "success":
            return False
        if calc.qc.get("engine_is_dummy") or calc.qc.get("fallback_dummy"):
            return False
        engine = (calc.method or {}).get("engine") or calc.data.get("engine")
        return engine not in {"dummy", "xtb"}

    @staticmethod
    def _calc_task(calc: Artifact) -> str:
        return str(calc.data.get("task") or (calc.method or {}).get("task") or "")

    @staticmethod
    def _has_real_irc(irc_artifact: Artifact | None) -> bool:
        if irc_artifact is None:
            return False
        if irc_artifact.qc.get("irc_backend_is_dummy") or irc_artifact.qc.get("fallback_dummy"):
            return False
        if irc_artifact.data.get("irc_backend_is_dummy") or irc_artifact.data.get("fallback_dummy"):
            return False
        real = bool(
            irc_artifact.qc.get("real_irc_executed", False)
            or irc_artifact.qc.get("real_orca_executed", False)
            or irc_artifact.data.get("real_orca_executed", False)
        )
        return bool(irc_artifact.qc.get("irc_validated", False) and real)

    def run(self, manifest: Manifest | None, config: dict[str, Any], context: StageContext) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        engine_name = config.get("engine", "goodvibes")
        engine_settings = {**(config.get("settings", {}) or {}), **{k: v for k, v in config.items() if k not in {"name", "enabled", "engine", "settings"}}}
        thermo_engine = get_thermo_engine(engine_name, **engine_settings)

        calcs_by_species: dict[str, list[Artifact]] = {}
        for calc in manifest.iter_artifacts("calculation"):
            if calc.status.status != "success" or calc.data.get("electronic_energy_hartree") is None:
                continue
            sid = calc.data.get("species_id")
            if sid:
                calcs_by_species.setdefault(str(sid), []).append(calc)

        def calcs_for(species_id: str | None) -> list[Artifact]:
            return calcs_by_species.get(str(species_id or ""), [])

        thermal_tasks = {"opt_freq", "neb_ts", "optts", "ts_search", "neb_ts_optts_freq", "neb_ts_freq"}
        ts_tasks = {"opt_freq", "neb_ts", "optts", "ts_search", "neb_ts_optts_freq", "neb_ts_freq"}

        def select_energy_calc(species_id: str | None) -> Artifact | None:
            calcs = calcs_for(species_id)
            if not calcs:
                return None
            scientific_sp = [c for c in calcs if self._calc_task(c) == "single_point" and self._is_scientific_calc(c)]
            if scientific_sp:
                return scientific_sp[-1]
            scientific_opt = [c for c in calcs if (self._calc_task(c) in thermal_tasks or (c.method or {}).get("stage") in {"dft-minima", "ts-search"}) and self._is_scientific_calc(c)]
            if scientific_opt:
                return scientific_opt[-1]
            sp = [c for c in calcs if self._calc_task(c) == "single_point"]
            if sp:
                return sp[-1]
            return calcs[-1]

        def select_thermal_calc(species_id: str | None) -> Artifact | None:
            calcs = calcs_for(species_id)
            thermal = [c for c in calcs if self._calc_task(c) in thermal_tasks or (c.method or {}).get("stage") in {"dft-minima", "ts-search"}]
            return thermal[-1] if thermal else select_energy_calc(species_id)

        def combined_calc(species_id: str | None) -> Artifact | None:
            ecalc = select_energy_calc(species_id)
            tcalc = select_thermal_calc(species_id)
            if ecalc is None:
                return None
            data = dict(ecalc.data)
            if tcalc is not None:
                for key in [
                    "thermal_correction_gibbs_hartree",
                    "thermal_correction_enthalpy_hartree",
                    "gibbs_298K_hartree",
                    "enthalpy_298K_hartree",
                    "zpe_hartree",
                    "frequencies_cm1",
                    "n_imag",
                    "lowest_freq_cm1",
                    "imag_freq_cm1",
                    "hf_stretch_cm1",
                    "dipole_D",
                ]:
                    if key in tcalc.data and tcalc.data.get(key) is not None:
                        data[key] = tcalc.data[key]
            if "thermal_correction_gibbs_hartree" not in data and data.get("gibbs_298K_hartree") is not None:
                data["thermal_correction_gibbs_hartree"] = float(data["gibbs_298K_hartree"]) - float(data["electronic_energy_hartree"])
            data.update(
                {
                    "species_id": species_id,
                    "energy_source_calc_id": ecalc.artifact_id,
                    "thermal_source_calc_id": tcalc.artifact_id if tcalc else None,
                    "energy_source_is_scientific": self._is_scientific_calc(ecalc),
                    "thermal_source_is_scientific": self._is_scientific_calc(tcalc) if tcalc else False,
                }
            )
            qc = {**(ecalc.qc or {})}
            if tcalc is not None:
                qc.update({f"thermal_{k}": v for k, v in (tcalc.qc or {}).items() if k in {"fallback_dummy", "engine_is_dummy", "real_orca_executed"}})
            return Artifact(
                artifact_id=f"combined_thermo_source_{species_id}",
                artifact_type="calculation",
                parents=[x for x in [ecalc.artifact_id, tcalc.artifact_id if tcalc else None] if x],
                paths={**(ecalc.paths or {})},
                data=data,
                method={"engine": (ecalc.method or {}).get("engine"), "method_id": (ecalc.method or {}).get("method_id"), "task": "combined_sp_plus_freq"},
                qc=qc,
            )

        species_by_id: dict[str, Artifact] = {}
        for art_type in ["species", "species_preopt", "species_optimized"]:
            for sp in manifest.latest_artifacts(art_type):
                sid = sp.data.get("species_id") or sp.data.get("source_species_id") or sp.artifact_id
                species_by_id[str(sid)] = sp

        irc_by_reaction = {a.data.get("reaction_id"): a for a in manifest.latest_artifacts("irc") if a.status.status == "success"}
        reactions = list(manifest.latest_artifacts("reaction_validated")) or list(manifest.latest_artifacts("reaction"))

        temps = [float(t) for t in context.global_config.get("temperature_K", [298.15])]
        pstd = float(context.global_config.get("standard_pressure_bar", 1.0))
        partial = context.global_config.get("partial_pressures", {}) or {}
        p_hf = partial.get("HF")
        p_candidate = partial.get("candidate")
        p_complex = partial.get("complex")
        p_ts = partial.get("TS")

        # Build species thermo records only for species that appear in reactions.
        reaction_species: set[str] = set()
        for rxn in reactions:
            for key in ["reactant_species_id", "product_species_id", "ts_species_id", "isolated_candidate_species_id", "candidate_species_id", "hf_cluster_species_id"]:
                if rxn.data.get(key):
                    reaction_species.add(str(rxn.data[key]))

        species_thermo_index: dict[tuple[str, float, str], dict[str, Any]] = {}
        species_thermo_records: list[dict[str, Any]] = []

        def p_for_species(species_id: str, role: str) -> float | None:
            if role == "candidate":
                return p_candidate
            if role == "hf_cluster":
                return partial.get(f"HF_cluster_{species_by_id.get(species_id, Artifact(artifact_id='x', artifact_type='x')).data.get('hf_n')}") or partial.get("HF_cluster") or p_hf
            if role in {"reactant_complex", "ion_pair"}:
                return p_complex
            if role == "transition_state":
                return p_ts
            return None

        def make_species_thermo(species_id: str, T: float, role: str) -> dict[str, Any] | None:
            key = (species_id, float(T), role)
            if key in species_thermo_index:
                return species_thermo_index[key]
            calc = combined_calc(species_id)
            if calc is None:
                return None
            rec = thermo_engine.species_thermo(
                species_id=species_id,
                source_calc=calc,
                T_K=T,
                p_bar=p_for_species(species_id, role),
                p_standard_bar=pstd,
            )
            rec["role"] = role
            rec["state"] = (species_by_id.get(species_id) or Artifact(artifact_id="missing", artifact_type="missing")).data.get("state")
            rec["is_scientific_energy"] = bool(calc.data.get("energy_source_is_scientific"))
            rec["is_scientific_thermal"] = bool(calc.data.get("thermal_source_is_scientific"))
            species_thermo_index[key] = rec
            species_thermo_records.append(rec)
            out.add_artifact(
                Artifact(
                    artifact_id=f"species_thermo_{species_id}_{role}_{int(T)}K",
                    artifact_type="species_thermo",
                    parents=calc.parents,
                    data=rec,
                    qc={
                        "thermo_status": "success",
                        "quasi_rrho_applied": rec.get("quasi_rrho_applied"),
                        "low_frequency_count": rec.get("low_frequency_count"),
                        "energy_source_is_scientific": rec.get("is_scientific_energy"),
                        "thermal_source_is_scientific": rec.get("is_scientific_thermal"),
                        "production_thermo_ready": rec.get("production_thermo_ready"),
                        "external_goodvibes_executed": rec.get("external_goodvibes_executed"),
                        "external_goodvibes_status": rec.get("external_goodvibes_status"),
                    },
                )
            )
            return rec

        thermo_records: list[dict[str, Any]] = []
        for rxn in reactions:
            rc_id = rxn.data.get("reactant_species_id")
            ip_id = rxn.data.get("product_species_id")
            ts_id = rxn.data.get("ts_species_id")
            b_id = rxn.data.get("isolated_candidate_species_id") or rxn.data.get("candidate_species_id")
            hfn_id = rxn.data.get("hf_cluster_species_id")
            reaction_id = rxn.data.get("reaction_id", rxn.artifact_id)
            irc_artifact = irc_by_reaction.get(reaction_id)
            has_irc = self._has_real_irc(irc_artifact)

            has_dft = all(
                any(self._is_scientific_calc(c) for c in calcs_for(sid))
                for sid in [rc_id, ip_id]
                if sid
            )
            has_high_level_sp = all(
                any(self._calc_task(c) == "single_point" and self._is_scientific_calc(c) for c in calcs_for(sid))
                for sid in [rc_id, ip_id, ts_id]
                if sid
            )
            has_ts = bool(
                ts_id
                and any(
                    (self._calc_task(c) in ts_tasks or (c.method or {}).get("stage") == "ts-search")
                    and self._is_scientific_calc(c)
                    and bool(c.qc.get("ts_validated_by_frequency", False))
                    for c in calcs_for(ts_id)
                )
            )
            q_tier = quality_tier(has_dft=has_dft, has_ts=has_ts, has_irc=bool(has_ts and has_irc), has_high_level_sp=has_high_level_sp)
            raw_conf = confidence_score(1.0 if q_tier in {"Q4", "Q5"} else 0.7, 0.9, 1.0 if has_irc else 0.65, 0.5, 0.2)
            # If the main free-energy values are still development/dummy values,
            # keep confidence low even when the workflow wiring succeeded.
            main_values_are_dummy = not has_dft
            main_values_are_fallback = bool(main_values_are_dummy)
            conf = cap_confidence(raw_conf, q_tier, has_dummy_or_fallback=main_values_are_dummy or main_values_are_fallback)

            for T in temps:
                b = make_species_thermo(str(b_id), T, "candidate") if b_id else None
                hfn = make_species_thermo(str(hfn_id), T, "hf_cluster") if hfn_id else None
                rc = make_species_thermo(str(rc_id), T, "reactant_complex") if rc_id else None
                ip = make_species_thermo(str(ip_id), T, "ion_pair") if ip_id else None
                ts = make_species_thermo(str(ts_id), T, "transition_state") if ts_id else None
                if rc is None or ip is None:
                    out.add_artifact(Artifact.failure(f"thermo_failed_{reaction_id}_{int(T)}K", "thermo", "missing_rc_or_ip_energy", category="missing_input", parents=[rxn.artifact_id]))
                    continue
                dg_assoc_std = None if b is None or hfn is None else (rc["G_standard_hartree"] - b["G_standard_hartree"] - hfn["G_standard_hartree"]) * HARTREE_TO_KCAL_MOL
                dg_assoc_proc = None if b is None or hfn is None else (rc["G_standard_hartree"] - b["G_process_hartree"] - hfn["G_process_hartree"]) * HARTREE_TO_KCAL_MOL
                dg_ion = (ip["G_standard_hartree"] - rc["G_standard_hartree"]) * HARTREE_TO_KCAL_MOL
                dg_act = None if ts is None else (ts["G_standard_hartree"] - rc["G_standard_hartree"]) * HARTREE_TO_KCAL_MOL
                dg_act_proc = None if ts is None else (ts["G_process_hartree"] - rc["G_process_hartree"]) * HARTREE_TO_KCAL_MOL
                ts_calc = combined_calc(str(ts_id)) if ts_id else None
                imag_freq = ts_calc.data.get("imag_freq_cm1") if ts_calc else None
                rec = {
                    "reaction_id": reaction_id,
                    "mol_id": rxn.data.get("mol_id"),
                    "site_id": rxn.data.get("site_id"),
                    "site_type": rxn.data.get("site_type"),
                    "conformer_id": rxn.data.get("conformer_id"),
                    "hf_n": rxn.data.get("hf_n"),
                    "reaction_type": rxn.data.get("reaction_type"),
                    "T_K": T,
                    "p_standard_bar": pstd,
                    "p_HF_bar": p_hf,
                    "p_candidate_bar": p_candidate,
                    "p_complex_bar": p_complex,
                    "G_candidate_standard_hartree": None if b is None else b["G_standard_hartree"],
                    "G_hf_cluster_standard_hartree": None if hfn is None else hfn["G_standard_hartree"],
                    "G_reactant_complex_standard_hartree": rc["G_standard_hartree"],
                    "G_product_ionpair_standard_hartree": ip["G_standard_hartree"],
                    "G_TS_standard_hartree": None if ts is None else ts["G_standard_hartree"],
                    "G_candidate_process_hartree": None if b is None else b["G_process_hartree"],
                    "G_hf_cluster_process_hartree": None if hfn is None else hfn["G_process_hartree"],
                    "G_reactant_complex_process_hartree": rc["G_process_hartree"],
                    "G_product_ionpair_process_hartree": ip["G_process_hartree"],
                    "G_TS_process_hartree": None if ts is None else ts["G_process_hartree"],
                    # Backward-compatible aliases expected by rank/report code.
                    "G_reactant_complex_hartree": rc["G_standard_hartree"],
                    "G_product_ionpair_hartree": ip["G_standard_hartree"],
                    "G_TS_hartree": None if ts is None else ts["G_standard_hartree"],
                    "delta_G_assoc_standard_kcal_mol": dg_assoc_std,
                    "delta_G_assoc_process_kcal_mol": dg_assoc_proc,
                    "delta_G_assoc_kcal_mol": dg_assoc_std,
                    "delta_G_assoc_pressure_corrected_kcal_mol": dg_assoc_proc,
                    "delta_G_ionpair_kcal_mol": dg_ion,
                    "delta_G_act_kcal_mol": dg_act,
                    "delta_G_act_process_kcal_mol": dg_act_proc,
                    "K_assoc_standard": equilibrium_constant(dg_assoc_std, T),
                    "K_assoc_process_adjusted": equilibrium_constant(dg_assoc_proc, T),
                    "K_ionpair": equilibrium_constant(dg_ion, T),
                    "imag_freq_cm1": imag_freq,
                    "quasi_rrho_applied": bool(any(x and x.get("quasi_rrho_applied") for x in [b, hfn, rc, ip, ts])),
                    "low_frequency_count_total": int(sum(int(x.get("low_frequency_count", 0)) for x in [b, hfn, rc, ip, ts] if x)),
                    "thermo_backend": getattr(thermo_engine, "name", str(engine_name)),
                    "temperature_model": "source_298K_or_backend_corrected",
                    "pressure_model": "ideal_gas_RTlnp_reactant_adjusted_for_assoc",
                    "irc_artifact_id": irc_artifact.artifact_id if irc_artifact else None,
                    "irc_validated": bool(irc_artifact.qc.get("irc_validated")) if irc_artifact else False,
                    "real_irc_executed": bool(irc_artifact.qc.get("real_irc_executed") or irc_artifact.qc.get("real_orca_executed") or irc_artifact.data.get("real_orca_executed")) if irc_artifact else False,
                    "quality_tier": q_tier,
                    "confidence_score_raw": raw_conf,
                    "confidence_score": conf,
                    "has_scientific_dft": bool(has_dft),
                    "has_validated_ts": bool(has_ts),
                    "has_validated_irc": bool(has_ts and has_irc),
                    "has_high_level_sp": bool(has_high_level_sp),
                    "main_values_are_dummy": bool(main_values_are_dummy),
                    "main_values_are_fallback": bool(main_values_are_fallback),
                    "energy_source_reactant": rc.get("energy_source_calc_id"),
                    "energy_source_product": ip.get("energy_source_calc_id"),
                    "energy_source_ts": None if ts is None else ts.get("energy_source_calc_id"),
                    "thermal_source_reactant": rc.get("energy_source_calc_id"),
                    "association_reference_status": "computed_from_separated_components" if dg_assoc_std is not None else "missing_separated_component_energy",
                    "production_thermo_ready": bool(all((x or {}).get("production_thermo_ready") for x in [rc, ip] if x) and (ts is None or ts.get("production_thermo_ready"))),
                    "external_goodvibes_executed": bool(any((x or {}).get("external_goodvibes_executed") for x in [b, hfn, rc, ip, ts] if x)),
                    "external_goodvibes_statuses": ";".join(sorted({str((x or {}).get("external_goodvibes_status")) for x in [b, hfn, rc, ip, ts] if x and (x or {}).get("external_goodvibes_status")})),
                    "connector_quality": "production_goodvibes" if bool(all((x or {}).get("production_thermo_ready") for x in [rc, ip] if x) and (ts is None or ts.get("production_thermo_ready"))) else "screening_or_internal_thermo",
                }
                rec["scientific_rank_eligible"] = scientific_rank_eligible(q_tier, has_dummy_or_fallback=rec["main_values_are_dummy"] or rec["main_values_are_fallback"])
                rec["production_rank_eligible"] = production_rank_eligible(
                    q_tier,
                    has_dummy_or_fallback=rec["main_values_are_dummy"] or rec["main_values_are_fallback"],
                    production_thermo_ready=rec["production_thermo_ready"],
                    real_irc_executed=rec["real_irc_executed"],
                )
                rec["science_gate_reason"] = (
                    "production_ready" if rec["production_rank_eligible"] else
                    "scientific_dft_ready" if rec["scientific_rank_eligible"] else
                    "screening_only_dummy_or_no_real_dft" if rec["main_values_are_dummy"] or not rec["has_scientific_dft"] else
                    "needs_ts_irc_or_production_thermo"
                )
                thermo_records.append(rec)
                out.add_artifact(Artifact(artifact_id=f"thermo_{reaction_id}_{int(T)}K", artifact_type="thermo", parents=[rxn.artifact_id], data=rec, qc={"thermo_status": "success", "standard_state": f"{pstd}bar", "quality_tier": q_tier, "quasi_rrho_applied": rec["quasi_rrho_applied"]}))

        write_jsonl(species_thermo_records, out_dir / "species_thermo_records.jsonl")
        write_jsonl(thermo_records, out_dir / "thermo_records.jsonl")
        if species_thermo_records:
            pd.DataFrame(species_thermo_records).to_csv(out_dir / "species_thermo.csv", index=False)
        if thermo_records:
            pd.DataFrame(thermo_records).to_csv(out_dir / "reaction_thermo.csv", index=False)
        out.add_artifact(Artifact(artifact_id="species_thermo_table", artifact_type="table", paths={"csv": str(out_dir / "species_thermo.csv"), "jsonl": str(out_dir / "species_thermo_records.jsonl")}, data={"n_rows": int(len(species_thermo_records)), "table_type": "species_thermo"}))
        out.add_artifact(Artifact(artifact_id="reaction_thermo_table", artifact_type="table", paths={"csv": str(out_dir / "reaction_thermo.csv"), "jsonl": str(out_dir / "thermo_records.jsonl")}, data={"n_rows": int(len(thermo_records)), "table_type": "reaction_thermo"}))
        return out

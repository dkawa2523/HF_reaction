from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.backends.qm.nwchem import backfill_nwchem_thermochemistry_frequencies
from hfauto.backends.registry import get_thermo_engine
from hfauto.chemistry.minima import is_accepted_optimized_minimum
from hfauto.chemistry.stoichiometry import normalize_terms, stoichiometric_sum
from hfauto.core.hashing import sha256_file
from hfauto.core.io import ensure_dir, write_jsonl
from hfauto.core.qc import (
    cap_confidence,
    confidence_score,
    production_rank_eligible,
    quality_tier,
    scientific_rank_eligible,
)
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.core.schemas.thermo import ThermoRecord
from hfauto.core.thermo_models import equilibrium_constant
from hfauto.core.units import HARTREE_TO_KCAL_MOL
from hfauto.stages.base import Stage, StageContext

THERMAL_FREQUENCY_TASKS = frozenset(
    {
        "opt_freq",
        "saddle_freq",
        "neb_ts",
        "optts",
        "ts_search",
        "neb_ts_optts_freq",
        "neb_ts_freq",
    }
)


def _activation_difference_kcal_mol(
    transition_state_hartree: float | None,
    reactants_hartree: float | None,
) -> float | None:
    if transition_state_hartree is None or reactants_hartree is None:
        return None
    return (
        float(transition_state_hartree) - float(reactants_hartree)
    ) * HARTREE_TO_KCAL_MOL


def _publishable_electronic_activation(
    transition_state_hartree: float | None,
    reactants_hartree: float | None,
    resolution_flag: bool | None,
) -> dict[str, float | bool | None]:
    """Keep a sub-resolution difference as evidence, not a published barrier."""

    value = _activation_difference_kcal_mol(
        transition_state_hartree, reactants_hartree
    )
    return {
        "delta_E_activation_kcal_mol": (
            None if resolution_flag is False else value
        ),
        "delta_E_activation_unresolved_kcal_mol": (
            value if resolution_flag is False else None
        ),
        "activation_energy_resolved": resolution_flag,
    }


class ThermoStage(Stage):
    """Build species and reaction thermochemistry records.

    Thermochemistry is separated into two levels:
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
        if engine in {"dummy", "xtb"}:
            return False
        if calc.qc.get("real_qm_executed", calc.qc.get("real_orca_executed")) is not True:
            return False
        if calc.qc.get("scf_converged") is not True or calc.qc.get("normal_termination") is False:
            return False
        task = ThermoStage._calc_task(calc)
        if task in {"opt_freq", "saddle_freq"}:
            state = calc.data.get("state")
            expected = 1 if state == "transition_state" or task == "saddle_freq" else 0
            if calc.data.get("n_imag") != expected:
                return False
            if task == "opt_freq" and state != "transition_state" and (
                (calc.method or {}).get("stage") == "dft-minima"
                or "minimum_accepted" in calc.qc
            ):
                if calc.qc.get("minimum_accepted") is not True:
                    return False
                if state in {
                    "reactant_complex",
                    "shared_proton",
                    "ion_pair",
                } and calc.qc.get("proton_state_accepted") is not True:
                    return False
        return True

    @staticmethod
    def _is_validated_minimum_calc(calc: Artifact) -> bool:
        state = calc.data.get("state")
        return bool(
            ThermoStage._calc_task(calc) == "opt_freq"
            and state != "transition_state"
            and ThermoStage._is_scientific_calc(calc)
            and calc.qc.get("minimum_accepted") is True
            and (
                state not in {"reactant_complex", "shared_proton", "ion_pair"}
                or calc.qc.get("proton_state_accepted") is True
            )
        )

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

    @classmethod
    def _endpoint_frequency_adapter(
        cls,
        manifest: Manifest,
        species: Artifact,
        *,
        source_claim_count: int = 1,
    ) -> tuple[Artifact | None, list[str]]:
        """Bind one accepted endpoint species to its exact frequency source.

        Endpoint calculations are commonly launched from seed species whose
        temporary IDs differ from the canonical reaction endpoint.  This
        adapter changes only the in-memory thermochemistry view: it keeps the
        original calculation ID and files while inheriting the published
        minimum/proton/state-method gates from the accepted endpoint artifact.
        """

        reasons: list[str] = []
        if not is_accepted_optimized_minimum(species):
            reasons.append("endpoint_species_is_not_an_accepted_minimum")
        if species.qc.get("state_method_evidence_validated") is not True:
            reasons.append("endpoint_state_method_evidence_is_not_validated")

        data_source_id = (species.data or {}).get(
            "source_frequency_calculation_id"
        )
        provenance_source_id = (species.provenance or {}).get(
            "source_frequency_calculation_id"
        )
        source_ids = {
            str(value)
            for value in (data_source_id, provenance_source_id)
            if isinstance(value, (str, int))
            and not isinstance(value, bool)
            and str(value)
        }
        if len(source_ids) != 1:
            reasons.append("source_frequency_calculation_id_missing_or_conflicting")
            return None, reasons
        source_id = next(iter(source_ids))
        if str(data_source_id or "") != source_id or str(
            provenance_source_id or ""
        ) != source_id:
            reasons.append("source_frequency_calculation_id_not_explicit_twice")
        if source_id not in {str(parent) for parent in species.parents}:
            reasons.append("source_frequency_calculation_is_not_endpoint_parent")
        if source_claim_count != 1:
            reasons.append("source_frequency_calculation_claimed_by_multiple_species")

        source = manifest.find(source_id)
        if source is None:
            reasons.append("source_frequency_calculation_missing")
            return None, list(dict.fromkeys(reasons))
        if source.artifact_type != "calculation":
            reasons.append("source_frequency_artifact_is_not_calculation")
        if source.status.status != "success":
            reasons.append("source_frequency_calculation_not_successful")
        for key in (
            "real_qm_executed",
            "scf_converged",
            "geometry_converged",
            "normal_termination",
            "dispersion_applied",
            "geometry_sane",
        ):
            if source.qc.get(key) is not True:
                reasons.append(f"source_frequency_{key}_is_not_true")
        if source.qc.get("fallback_dummy") is not False:
            reasons.append("source_frequency_fallback_dummy_is_not_false")
        if source.qc.get("engine_is_dummy") is not False:
            reasons.append("source_frequency_engine_is_dummy_is_not_false")
        if cls._calc_task(source) != "opt_freq":
            reasons.append("source_frequency_task_is_not_opt_freq")
        if source.data.get("frequency_analysis_present") is not True:
            reasons.append("source_frequency_analysis_missing")
        if source.data.get("frequency_count_complete") is not True:
            reasons.append("source_frequency_count_incomplete")
        if source.data.get("temperature_consistent") is not True:
            reasons.append("source_frequency_temperature_inconsistent")
        if source.data.get("n_imag") != 0 or source.qc.get("n_imag") != 0:
            reasons.append("source_frequency_is_not_a_minimum")
        if not cls._is_scientific_calc(source):
            reasons.append("source_frequency_calculation_is_not_scientific")

        for key in ("resolved_charge", "resolved_multiplicity", "electron_count"):
            expected = species.data.get(key)
            if expected is not None and source.data.get(key) != expected:
                reasons.append(f"source_frequency_{key}_mismatch")

        endpoint_xyz = (
            species.paths.get("final_xyz")
            or species.paths.get("xyz")
            or species.data.get("xyz_path")
        )
        source_xyz = source.paths.get("final_xyz")
        try:
            geometry_matches = bool(
                endpoint_xyz
                and source_xyz
                and Path(str(endpoint_xyz)).is_file()
                and Path(str(source_xyz)).is_file()
                and sha256_file(endpoint_xyz) == sha256_file(source_xyz)
            )
        except OSError:
            geometry_matches = False
        if not geometry_matches:
            reasons.append("source_frequency_final_geometry_mismatch")

        reasons = list(dict.fromkeys(reasons))
        if reasons:
            return None, reasons

        inherited_qc_keys = {
            "is_minimum",
            "minimum_accepted",
            "proton_state_accepted",
            "state_method_evidence_validated",
            "scf_converged",
            "geometry_converged",
            "normal_termination",
            "real_qm_executed",
            "engine_is_dummy",
            "dispersion_applied",
            "geometry_sane",
            "atom_count_ok",
            "collision_detected",
            "frequency_count_complete",
            "frequency_analysis_present",
            "temperature_consistent",
            "fallback_dummy",
            "n_imag",
        }
        adapter_qc = dict(source.qc or {})
        adapter_qc.update(
            {
                key: species.qc.get(key)
                for key in inherited_qc_keys
                if key in species.qc
            }
        )
        species_id = str(
            species.data.get("species_id")
            or species.data.get("source_species_id")
            or species.artifact_id
        )
        adapter_data = dict(source.data or {})
        if adapter_data.get("thermochemistry_frequencies_cm1") is None:
            migrated_frequencies, migration_reasons = (
                backfill_nwchem_thermochemistry_frequencies(
                    source,
                    species.provenance.get("source_frequency_file_sha256"),
                )
            )
            if migration_reasons or migrated_frequencies is None:
                return None, migration_reasons or [
                    "thermochemistry_frequency_backfill_failed"
                ]
            adapter_data["thermochemistry_frequencies_cm1"] = migrated_frequencies
            adapter_data["thermochemistry_frequency_cutoff_cm1"] = 1.0
            adapter_data[
                "thermochemistry_frequencies_migrated_from_raw"
            ] = True
        return (
            Artifact(
                artifact_id=f"thermo_endpoint_adapter_{species.artifact_id}",
                artifact_type="calculation",
                parents=[species.artifact_id, source.artifact_id],
                paths=dict(source.paths or {}),
                data={
                    **adapter_data,
                    "species_id": species_id,
                    "state": species.data.get("state"),
                    "endpoint_species_artifact_id": species.artifact_id,
                    "source_frequency_calculation_id": source.artifact_id,
                    "energy_source_calc_id": source.artifact_id,
                    "thermal_source_calc_id": source.artifact_id,
                    "thermal_source_lineage_accepted": True,
                },
                method=dict(source.method or {}),
                provenance={
                    **(source.provenance or {}),
                    "thermo_adapter": True,
                    "endpoint_species_artifact_id": species.artifact_id,
                    "source_frequency_calculation_id": source.artifact_id,
                },
                qc=adapter_qc,
            ),
            [],
        )

    def run(self, manifest: Manifest | None, config: dict[str, Any], context: StageContext) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        engine_name = config.get("engine", "goodvibes")
        engine_settings = {**(config.get("settings", {}) or {}), **{k: v for k, v in config.items() if k not in {"name", "enabled", "engine", "settings"}}}
        engine_settings.setdefault("work_root", str(out_dir / "backend_runs"))
        thermo_engine = get_thermo_engine(engine_name, **engine_settings)

        species_by_id: dict[str, Artifact] = {}
        for art_type in ["species", "species_preopt", "species_optimized"]:
            for sp in manifest.latest_artifacts(art_type):
                sid = (
                    sp.data.get("species_id")
                    or sp.data.get("source_species_id")
                    or sp.artifact_id
                )
                species_by_id[str(sid)] = sp

        explicit_endpoints_by_species: dict[str, list[Artifact]] = {}
        source_claims: dict[str, set[str]] = {}
        for sp in manifest.latest_artifacts("species_optimized"):
            source_values = (
                sp.data.get("source_frequency_calculation_id"),
                sp.provenance.get("source_frequency_calculation_id"),
            )
            if not any(value is not None and value != "" for value in source_values):
                continue
            sid = str(
                sp.data.get("species_id")
                or sp.data.get("source_species_id")
                or sp.artifact_id
            )
            explicit_endpoints_by_species.setdefault(sid, []).append(sp)
            for value in source_values:
                if isinstance(value, (str, int)) and not isinstance(value, bool):
                    source_claims.setdefault(str(value), set()).add(sp.artifact_id)

        endpoint_thermal_by_species: dict[str, Artifact] = {}
        endpoint_lineage_failures: dict[str, list[str]] = {}
        for sid, endpoints in explicit_endpoints_by_species.items():
            if len(endpoints) != 1:
                endpoint_lineage_failures[sid] = [
                    "multiple_explicit_endpoint_species_for_canonical_species_id"
                ]
                continue
            endpoint = endpoints[0]
            source_id = str(
                endpoint.data.get("source_frequency_calculation_id")
                or endpoint.provenance.get("source_frequency_calculation_id")
                or ""
            )
            adapter, reasons = self._endpoint_frequency_adapter(
                manifest,
                endpoint,
                source_claim_count=len(source_claims.get(source_id, set())),
            )
            if adapter is None:
                endpoint_lineage_failures[sid] = reasons
            else:
                endpoint_thermal_by_species[sid] = adapter

        for sid, reasons in endpoint_lineage_failures.items():
            endpoint = explicit_endpoints_by_species[sid][-1]
            out.add_artifact(
                Artifact.failure(
                    f"thermo_endpoint_lineage_failed_{endpoint.artifact_id}",
                    "thermo_lineage",
                    ";".join(reasons),
                    category="invalid_endpoint_frequency_lineage",
                    parents=[endpoint.artifact_id],
                    recoverable=False,
                    data={"species_id": sid, "reasons": reasons},
                )
            )

        calcs_by_species: dict[str, list[Artifact]] = {}
        for calc in manifest.iter_artifacts("calculation"):
            if calc.status.status != "success" or calc.data.get("electronic_energy_hartree") is None:
                continue
            sid = calc.data.get("species_id")
            if sid:
                calcs_by_species.setdefault(str(sid), []).append(calc)

        def calcs_for(species_id: str | None) -> list[Artifact]:
            return calcs_by_species.get(str(species_id or ""), [])

        claimed_endpoint_source_ids = set(source_claims)

        def validation_calcs_for(species_id: str | None) -> list[Artifact]:
            sid = str(species_id or "")
            adapter = endpoint_thermal_by_species.get(sid)
            return [*calcs_for(sid), *([adapter] if adapter is not None else [])]

        thermal_tasks = THERMAL_FREQUENCY_TASKS
        ts_tasks = set(thermal_tasks)

        def select_energy_calc(species_id: str | None) -> Artifact | None:
            sid = str(species_id or "")
            calcs = calcs_for(sid)
            endpoint_adapter = endpoint_thermal_by_species.get(sid)
            if not calcs:
                return endpoint_adapter
            scientific_sp = [c for c in calcs if self._calc_task(c) == "single_point" and self._is_scientific_calc(c)]
            if scientific_sp:
                return scientific_sp[-1]
            scientific_opt = [c for c in calcs if (self._calc_task(c) in thermal_tasks or (c.method or {}).get("stage") in {"dft-minima", "ts-search"}) and self._is_scientific_calc(c)]
            if scientific_opt:
                return scientific_opt[-1]
            sp = [c for c in calcs if self._calc_task(c) == "single_point"]
            if sp:
                return sp[-1]
            return endpoint_adapter or calcs[-1]

        def select_thermal_calc(species_id: str | None) -> Artifact | None:
            sid = str(species_id or "")
            if sid in explicit_endpoints_by_species:
                return endpoint_thermal_by_species.get(sid)
            calcs = [
                calc
                for calc in calcs_for(sid)
                if calc.artifact_id not in claimed_endpoint_source_ids
            ]
            thermal = [c for c in calcs if self._calc_task(c) in thermal_tasks or (c.method or {}).get("stage") in {"dft-minima", "ts-search"}]
            if thermal:
                return thermal[-1]
            energy_calc = select_energy_calc(species_id)
            if (
                energy_calc is not None
                and energy_calc.artifact_id in claimed_endpoint_source_ids
            ):
                return None
            return energy_calc

        def combined_calc(species_id: str | None) -> Artifact | None:
            ecalc = select_energy_calc(species_id)
            tcalc = select_thermal_calc(species_id)
            if ecalc is None or (
                str(species_id or "") in explicit_endpoints_by_species
                and tcalc is None
            ) or (
                tcalc is None
                and ecalc.artifact_id in claimed_endpoint_source_ids
            ):
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
                    "thermochemistry_frequencies_cm1",
                    "thermochemistry_frequency_cutoff_cm1",
                    "raw_frequencies_cm1",
                    "raw_frequency_count",
                    "expected_raw_frequency_count",
                    "frequency_count_complete",
                    "frequency_analysis_present",
                    "thermochemistry_temperature_K",
                    "requested_temperature_K",
                    "temperature_consistent",
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
                    "energy_source_calc_id": ecalc.data.get(
                        "energy_source_calc_id", ecalc.artifact_id
                    ),
                    "thermal_source_calc_id": (
                        tcalc.data.get("thermal_source_calc_id", tcalc.artifact_id)
                        if tcalc
                        else None
                    ),
                    "energy_source_is_scientific": self._is_scientific_calc(ecalc),
                    "thermal_source_is_scientific": self._is_scientific_calc(tcalc) if tcalc else False,
                    "thermal_source_lineage_accepted": bool(
                        tcalc is not None
                        and tcalc.data.get("thermal_source_lineage_accepted", True)
                    ),
                    "endpoint_species_artifact_id": (
                        tcalc.data.get("endpoint_species_artifact_id")
                        if tcalc is not None
                        else None
                    ),
                }
            )
            qc = {**(ecalc.qc or {})}
            if tcalc is not None:
                qc.update({f"thermal_{k}": v for k, v in (tcalc.qc or {}).items() if k in {"fallback_dummy", "engine_is_dummy", "real_qm_executed", "real_orca_executed"}})
            energy_parent = str(
                ecalc.data.get("energy_source_calc_id") or ecalc.artifact_id
            )
            thermal_parent = (
                str(tcalc.data.get("thermal_source_calc_id") or tcalc.artifact_id)
                if tcalc is not None
                else None
            )
            endpoint_parent = (
                str(tcalc.data.get("endpoint_species_artifact_id"))
                if tcalc is not None
                and tcalc.data.get("endpoint_species_artifact_id")
                else None
            )
            return Artifact(
                artifact_id=f"combined_thermo_source_{species_id}",
                artifact_type="calculation",
                parents=list(
                    dict.fromkeys(
                        x
                        for x in [energy_parent, thermal_parent, endpoint_parent]
                        if x
                    )
                ),
                paths={**((tcalc or ecalc).paths or {})},
                data=data,
                method={
                    "engine": (ecalc.method or {}).get("engine"),
                    "method_id": (ecalc.method or {}).get("method_id"),
                    "thermal_engine": (tcalc.method or {}).get("engine")
                    if tcalc is not None
                    else None,
                    "task": "combined_sp_plus_freq",
                },
                qc=qc,
            )

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
            if species_id in partial:
                return partial[species_id]
            if role in partial:
                return partial[role]
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
            rec["thermal_source_lineage_accepted"] = bool(
                calc.data.get("thermal_source_lineage_accepted", True)
            )
            rec["endpoint_species_artifact_id"] = calc.data.get(
                "endpoint_species_artifact_id"
            )
            if not rec["thermal_source_lineage_accepted"]:
                rec["production_thermo_ready"] = False
                rec["quantitative_thermo_ready"] = False
                reasons = list(rec.get("production_thermo_not_ready_reasons") or [])
                reasons.append("endpoint_frequency_lineage_not_accepted")
                rec["production_thermo_not_ready_reasons"] = list(
                    dict.fromkeys(reasons)
                )
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
                        "thermal_source_lineage_accepted": rec.get(
                            "thermal_source_lineage_accepted"
                        ),
                        "production_thermo_ready": rec.get("production_thermo_ready"),
                        "external_goodvibes_executed": rec.get("external_goodvibes_executed"),
                        "external_goodvibes_status": rec.get("external_goodvibes_status"),
                    },
                )
            )
            return rec

        thermo_records: list[dict[str, Any]] = []
        for rxn in reactions:
            stoichiometry = rxn.data.get("stoichiometry") or {}
            if stoichiometry:
                reactant_terms = normalize_terms(stoichiometry.get("reactants") or [])
                product_terms = normalize_terms(stoichiometry.get("products") or [])
                reaction_id = rxn.data.get("reaction_id", rxn.artifact_id)
                ts_id = rxn.data.get("ts_species_id")
                irc_artifact = irc_by_reaction.get(reaction_id)
                has_irc = self._has_real_irc(irc_artifact)
                minimum_ids = list(
                    dict.fromkeys(
                        term.species_id for term in [*reactant_terms, *product_terms]
                    )
                )
                calculation_ids = [*minimum_ids, *([str(ts_id)] if ts_id else [])]
                has_dft = bool(minimum_ids) and all(
                    any(
                        self._is_validated_minimum_calc(calc)
                        for calc in validation_calcs_for(species_id)
                    )
                    for species_id in minimum_ids
                )
                has_high_level_sp = bool(minimum_ids) and all(
                    any(
                        self._calc_task(calc) == "single_point"
                        and self._is_scientific_calc(calc)
                        for calc in calcs_for(species_id)
                    )
                    for species_id in calculation_ids
                )
                has_ts = bool(
                    ts_id
                    and any(
                        (
                            self._calc_task(calc) in ts_tasks
                            or (calc.method or {}).get("stage") == "ts-search"
                        )
                        and self._is_scientific_calc(calc)
                        and calc.qc.get("ts_validated_by_frequency") is True
                        for calc in calcs_for(str(ts_id))
                    )
                )
                q_tier = quality_tier(
                    has_dft=has_dft,
                    has_ts=has_ts,
                    has_irc=bool(has_ts and has_irc),
                    has_high_level_sp=has_high_level_sp,
                )
                raw_confidence = confidence_score(
                    1.0 if q_tier in {"Q4", "Q5"} else 0.7,
                    0.9,
                    1.0 if has_irc else 0.65,
                    0.5,
                    0.2,
                )
                main_values_are_dummy = not has_dft
                confidence = cap_confidence(
                    raw_confidence,
                    q_tier,
                    has_dummy_or_fallback=main_values_are_dummy,
                )

                for T in temps:
                    thermo_by_species: dict[str, dict[str, Any]] = {}
                    for term in [*reactant_terms, *product_terms]:
                        if term.species_id in thermo_by_species:
                            continue
                        species_thermo = make_species_thermo(
                            term.species_id, T, term.role or "reaction_participant"
                        )
                        if species_thermo is not None:
                            thermo_by_species[term.species_id] = species_thermo
                    ts_thermo = (
                        make_species_thermo(str(ts_id), T, "transition_state")
                        if ts_id
                        else None
                    )
                    missing = [
                        species_id
                        for species_id in minimum_ids
                        if species_id not in thermo_by_species
                    ]
                    if missing:
                        out.add_artifact(
                            Artifact.failure(
                                f"thermo_failed_{reaction_id}_{int(T)}K",
                                "thermo",
                                "missing thermochemistry for: " + ", ".join(missing),
                                category="missing_input",
                                parents=[rxn.artifact_id],
                                data={"reaction_id": reaction_id, "missing_species_ids": missing},
                            )
                        )
                        continue

                    def side_value(
                        terms: list[Any],
                        key: str,
                        values: dict[str, dict[str, Any]] = thermo_by_species,
                    ) -> float | None:
                        return stoichiometric_sum(
                            terms,
                            lambda term_species_id: values[term_species_id].get(key),
                        )

                    reactant_g = side_value(reactant_terms, "G_standard_hartree")
                    product_g = side_value(product_terms, "G_standard_hartree")
                    reactant_process_g = side_value(reactant_terms, "G_process_hartree")
                    product_process_g = side_value(product_terms, "G_process_hartree")
                    delta_g = (
                        None
                        if reactant_g is None or product_g is None
                        else (product_g - reactant_g) * HARTREE_TO_KCAL_MOL
                    )
                    delta_g_process = (
                        None
                        if reactant_process_g is None or product_process_g is None
                        else (product_process_g - reactant_process_g)
                        * HARTREE_TO_KCAL_MOL
                    )
                    delta_g_activation = (
                        None
                        if ts_thermo is None or reactant_g is None
                        else (ts_thermo["G_standard_hartree"] - reactant_g)
                        * HARTREE_TO_KCAL_MOL
                    )
                    delta_g_activation_process = (
                        None
                        if ts_thermo is None or reactant_process_g is None
                        else (ts_thermo["G_process_hartree"] - reactant_process_g)
                        * HARTREE_TO_KCAL_MOL
                    )
                    reactant_e = side_value(
                        reactant_terms, "electronic_energy_hartree"
                    )
                    reactant_e0 = side_value(
                        reactant_terms, "zero_point_corrected_energy_hartree"
                    )
                    reactant_h = side_value(
                        reactant_terms, "enthalpy_standard_hartree"
                    )
                    ts_e = (
                        None
                        if ts_thermo is None
                        else ts_thermo.get("electronic_energy_hartree")
                    )
                    ts_e0 = (
                        None
                        if ts_thermo is None
                        else ts_thermo.get(
                            "zero_point_corrected_energy_hartree"
                        )
                    )
                    ts_h = (
                        None
                        if ts_thermo is None
                        else ts_thermo.get("enthalpy_standard_hartree")
                    )

                    required_thermo = [
                        *thermo_by_species.values(),
                        *([ts_thermo] if ts_thermo is not None else []),
                    ]
                    quantitative_ready = all(
                        item.get("quantitative_thermo_ready") for item in required_thermo
                    )
                    production_ready = bool(
                        quantitative_ready
                        and all(
                            item.get("production_thermo_ready")
                            for item in required_thermo
                        )
                    )
                    not_ready_reasons = sorted(
                        {
                            str(reason)
                            for item in required_thermo
                            for reason in item.get(
                                "production_thermo_not_ready_reasons", []
                            )
                        }
                    )
                    ts_calc = combined_calc(str(ts_id)) if ts_id else None
                    activation_resolution = (
                        None
                        if ts_calc is None
                        else ts_calc.data.get(
                            "activation_energy_resolved",
                            ts_calc.qc.get("activation_energy_resolved"),
                        )
                    )
                    electronic_activation = _publishable_electronic_activation(
                        ts_e, reactant_e, activation_resolution
                    )
                    record = {
                        "reaction_id": reaction_id,
                        "reaction_type": rxn.data.get("reaction_type"),
                        "mechanism_family": rxn.data.get("mechanism_family"),
                        "T_K": T,
                        "p_standard_bar": pstd,
                        "stoichiometry": {
                            "reactants": [term.model_dump() for term in reactant_terms],
                            "products": [term.model_dump() for term in product_terms],
                        },
                        "G_reactants_standard_hartree": reactant_g,
                        "G_products_standard_hartree": product_g,
                        "G_TS_standard_hartree": (
                            None if ts_thermo is None else ts_thermo["G_standard_hartree"]
                        ),
                        "E_reactants_hartree": reactant_e,
                        "E_TS_hartree": ts_e,
                        "E0_reactants_hartree": reactant_e0,
                        "E0_TS_hartree": ts_e0,
                        "H_reactants_standard_hartree": reactant_h,
                        "H_TS_standard_hartree": ts_h,
                        **electronic_activation,
                        "delta_E0_activation_kcal_mol": _activation_difference_kcal_mol(
                            ts_e0, reactant_e0
                        ),
                        "delta_H_activation_standard_kcal_mol": (
                            _activation_difference_kcal_mol(ts_h, reactant_h)
                        ),
                        "delta_G_reaction_standard_kcal_mol": delta_g,
                        "delta_G_reaction_process_kcal_mol": delta_g_process,
                        "delta_G_activation_standard_kcal_mol": delta_g_activation,
                        "delta_G_activation_process_kcal_mol": delta_g_activation_process,
                        "delta_G_reaction_kcal_mol": delta_g,
                        "delta_G_act_kcal_mol": delta_g_activation,
                        "K_reaction_standard": equilibrium_constant(delta_g, T),
                        "imag_freq_cm1": (
                            None if ts_calc is None else ts_calc.data.get("imag_freq_cm1")
                        ),
                        "thermo_backend": getattr(thermo_engine, "name", str(engine_name)),
                        "pressure_model": "ideal_gas_RTlnp_per_stoichiometric_species",
                        "energy_source_ids": {
                            species_id: item.get("energy_source_calc_id")
                            for species_id, item in thermo_by_species.items()
                        },
                        "thermal_source_ids": {
                            species_id: item.get("thermal_source_calc_id")
                            for species_id, item in thermo_by_species.items()
                        },
                        "irc_artifact_id": (
                            None if irc_artifact is None else irc_artifact.artifact_id
                        ),
                        "has_scientific_dft": has_dft,
                        "has_validated_ts": has_ts,
                        "has_validated_irc": bool(has_ts and has_irc),
                        "activation_connectivity_validated": bool(
                            has_ts and has_irc
                        ),
                        "has_high_level_sp": has_high_level_sp,
                        "quality_tier": q_tier,
                        "confidence_score_raw": raw_confidence,
                        "confidence_score": confidence,
                        "main_values_are_dummy": main_values_are_dummy,
                        "main_values_are_fallback": main_values_are_dummy,
                        "quantitative_thermo_ready": quantitative_ready,
                        "production_thermo_ready": production_ready,
                        "production_thermo_not_ready_reasons": not_ready_reasons,
                    }
                    record["scientific_rank_eligible"] = scientific_rank_eligible(
                        q_tier, has_dummy_or_fallback=main_values_are_dummy
                    )
                    record["production_rank_eligible"] = production_rank_eligible(
                        q_tier,
                        has_dummy_or_fallback=main_values_are_dummy,
                        production_thermo_ready=production_ready,
                        real_irc_executed=has_irc,
                    )
                    record = ThermoRecord.model_validate(record).model_dump(
                        exclude_defaults=True, exclude_none=True
                    )
                    thermo_records.append(record)
                    out.add_artifact(
                        Artifact(
                            artifact_id=f"thermo_{reaction_id}_{int(T)}K",
                            artifact_type="thermo",
                            parents=[rxn.artifact_id],
                            data=record,
                            qc={
                                "thermo_status": "success",
                                "standard_state": f"{pstd}bar",
                                "quality_tier": q_tier,
                            },
                        )
                    )
                continue

            rc_id = rxn.data.get("reactant_species_id")
            ip_id = rxn.data.get("product_species_id")
            ts_id = rxn.data.get("ts_species_id")
            b_id = rxn.data.get("isolated_candidate_species_id") or rxn.data.get("candidate_species_id")
            hfn_id = rxn.data.get("hf_cluster_species_id")
            reaction_id = rxn.data.get("reaction_id", rxn.artifact_id)
            irc_artifact = irc_by_reaction.get(reaction_id)
            has_irc = self._has_real_irc(irc_artifact)

            required_minima_ids = [b_id, hfn_id, rc_id, ip_id]
            sp_species_ids = [sid for sid in [*required_minima_ids, ts_id] if sid]
            has_dft = all(required_minima_ids) and all(
                any(
                    self._is_validated_minimum_calc(c)
                    for c in validation_calcs_for(sid)
                )
                for sid in required_minima_ids
            )
            has_high_level_sp = all(required_minima_ids) and all(
                any(self._calc_task(c) == "single_point" and self._is_scientific_calc(c) for c in calcs_for(sid))
                for sid in sp_species_ids
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
                de0_act = _activation_difference_kcal_mol(
                    None
                    if ts is None
                    else ts.get("zero_point_corrected_energy_hartree"),
                    rc.get("zero_point_corrected_energy_hartree"),
                )
                dh_act = _activation_difference_kcal_mol(
                    None if ts is None else ts.get("enthalpy_standard_hartree"),
                    rc.get("enthalpy_standard_hartree"),
                )
                ts_calc = combined_calc(str(ts_id)) if ts_id else None
                activation_resolution = (
                    None
                    if ts_calc is None
                    else ts_calc.data.get(
                        "activation_energy_resolved",
                        ts_calc.qc.get("activation_energy_resolved"),
                    )
                )
                electronic_activation = _publishable_electronic_activation(
                    None if ts is None else ts.get("electronic_energy_hartree"),
                    rc.get("electronic_energy_hartree"),
                    activation_resolution,
                )
                imag_freq = ts_calc.data.get("imag_freq_cm1") if ts_calc else None
                required_thermo = [b, hfn, rc, ip] + ([ts] if ts_id else [])
                quantitative_thermo_ready = all(
                    x is not None and x.get("quantitative_thermo_ready")
                    for x in required_thermo
                )
                production_thermo_ready = bool(
                    quantitative_thermo_ready
                    and all(
                        x is not None and x.get("production_thermo_ready")
                        for x in required_thermo
                    )
                )
                thermo_not_ready_reasons = sorted(
                    {
                        str(reason)
                        for item in required_thermo
                        if item is not None
                        for reason in (
                            item.get("production_thermo_not_ready_reasons") or []
                        )
                    }
                )
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
                    "delta_G_activation_standard_kcal_mol": dg_act,
                    "delta_G_activation_process_kcal_mol": dg_act_proc,
                    **electronic_activation,
                    "delta_E0_activation_kcal_mol": de0_act,
                    "delta_H_activation_standard_kcal_mol": dh_act,
                    "K_assoc_standard": equilibrium_constant(dg_assoc_std, T),
                    "K_assoc_process_adjusted": equilibrium_constant(dg_assoc_proc, T),
                    "K_ionpair": equilibrium_constant(dg_ion, T),
                    "imag_freq_cm1": imag_freq,
                    "quasi_rrho_applied": bool(any(x and x.get("quasi_rrho_applied") for x in [b, hfn, rc, ip, ts])),
                    "low_frequency_count_total": int(sum(int(x.get("low_frequency_count", 0)) for x in [b, hfn, rc, ip, ts] if x)),
                    "thermo_backend": getattr(thermo_engine, "name", str(engine_name)),
                    "temperature_model": "per_species_fail_closed_temperature_gate",
                    "thermal_temperature_statuses": ";".join(
                        sorted(
                            {
                                str(item.get("thermal_temperature_status"))
                                for item in required_thermo
                                if item is not None
                                and item.get("thermal_temperature_status")
                            }
                        )
                    ),
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
                    "activation_connectivity_validated": bool(has_ts and has_irc),
                    "has_high_level_sp": bool(has_high_level_sp),
                    "main_values_are_dummy": bool(main_values_are_dummy),
                    "main_values_are_fallback": bool(main_values_are_fallback),
                    "energy_source_reactant": rc.get("energy_source_calc_id"),
                    "energy_source_product": ip.get("energy_source_calc_id"),
                    "energy_source_ts": None if ts is None else ts.get("energy_source_calc_id"),
                    "thermal_source_reactant": rc.get("thermal_source_calc_id"),
                    "association_reference_status": "computed_from_separated_components" if dg_assoc_std is not None else "missing_separated_component_energy",
                    "quantitative_thermo_ready": quantitative_thermo_ready,
                    "production_thermo_ready": production_thermo_ready,
                    "production_thermo_not_ready_reasons": thermo_not_ready_reasons,
                    "external_goodvibes_executed": bool(any((x or {}).get("external_goodvibes_executed") for x in [b, hfn, rc, ip, ts] if x)),
                    "external_goodvibes_statuses": ";".join(sorted({str((x or {}).get("external_goodvibes_status")) for x in [b, hfn, rc, ip, ts] if x and (x or {}).get("external_goodvibes_status")})),
                    "connector_quality": "production_goodvibes" if production_thermo_ready else "screening_or_internal_thermo",
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
        out.add_artifact(Artifact(artifact_id="species_thermo_table", artifact_type="table", paths={"csv": str(out_dir / "species_thermo.csv"), "jsonl": str(out_dir / "species_thermo_records.jsonl")}, data={"n_rows": len(species_thermo_records), "table_type": "species_thermo"}))
        out.add_artifact(Artifact(artifact_id="reaction_thermo_table", artifact_type="table", paths={"csv": str(out_dir / "reaction_thermo.csv"), "jsonl": str(out_dir / "thermo_records.jsonl")}, data={"n_rows": len(thermo_records), "table_type": "reaction_thermo"}))
        return out

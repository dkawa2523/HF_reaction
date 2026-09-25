from __future__ import annotations

import json
from typing import Any

import numpy as np
import pandas as pd

from hfauto.core.io import ensure_dir
from hfauto.core.public_data import process_penalty_from_public_data
from hfauto.core.qc import (
    cap_confidence,
    production_rank_eligible,
    scientific_rank_eligible,
    tier_numeric,
)
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext
from hfauto.stages.reaction_rank import ReactionRankStage


def _z(series: pd.Series, neutral: float = 0.0) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce")
    if s.notna().sum() == 0:
        return pd.Series([neutral] * len(series), index=series.index, dtype=float)
    s = s.fillna(s.median())
    std = float(s.std(ddof=0))
    if std < 1e-12:
        return pd.Series([neutral] * len(series), index=series.index, dtype=float)
    return (s - float(s.mean())) / std


def _moderate_binding_score(delta_g_assoc: pd.Series) -> pd.Series:
    x = pd.to_numeric(delta_g_assoc, errors="coerce").fillna(0.0)
    return -((x + 8.0).abs())


def _dict_value(obj: Any, *keys: str, default=None):
    cur = obj if isinstance(obj, dict) else {}
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


class RankStage(Stage):
    name = "rank"

    def run(self, manifest: Manifest | None, config: dict[str, Any], context: StageContext) -> Manifest:
        assert manifest is not None
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.name)
        thermo = [a.data for a in manifest.iter_artifacts("thermo")]
        desc = [a.data for a in manifest.iter_artifacts("descriptor")]
        kin = [a.data for a in manifest.iter_artifacts("kinetics")]
        mols = [a.data for a in list(manifest.iter_artifacts("molecule_enriched")) or list(manifest.iter_artifacts("molecule"))]
        failures = [a for a in manifest.artifacts if a.status.status == "failed"]
        if not thermo:
            out.add_artifact(Artifact.failure("rank_failed", "ranking", "no_thermo_records", category="missing_input"))
            return out
        generic_thermo = [record for record in thermo if record.get("stoichiometry")]
        if generic_thermo and len(generic_thermo) != len(thermo):
            out.add_artifact(
                Artifact.failure(
                    "rank_failed_mixed_domains",
                    "ranking",
                    "generic stoichiometric and legacy HF thermochemistry cannot share one ranking",
                    category="incomparable_reaction_domains",
                    recoverable=True,
                    recommended_fallback=(
                        "run reaction-rank for generic reactions and rank for the legacy HF study "
                        "in separate manifests"
                    ),
                )
            )
            return out
        if generic_thermo:
            return ReactionRankStage(output_stage_name=self.name).run(
                manifest, config, context
            )

        df = pd.DataFrame(thermo)
        dfd = pd.DataFrame(desc)
        if not dfd.empty:
            df = df.merge(dfd, on=["mol_id", "site_id", "site_type", "hf_n"], how="left")
        dfk = pd.DataFrame(kin)
        if not dfk.empty:
            kin_cols = [c for c in [
                "reaction_id", "T_K", "k_TST_s-1", "k_corrected_s-1", "transmission_coefficient",
                "tunneling_model", "K_assoc_standard", "K_assoc_process_adjusted", "K_ionpair",
                "arrhenius_A_s-1", "arrhenius_Ea_kcal_mol", "arrhenius_r2", "kinetics_quality"
            ] if c in dfk.columns]
            df = df.merge(dfk[kin_cols], on=["reaction_id", "T_K"], how="left")

        dfm = pd.DataFrame(mols)
        if not dfm.empty:
            keep = [c for c in ["mol_id", "name", "canonical_smiles", "inchikey", "identity", "public_data"] if c in dfm.columns]
            df = df.merge(dfm[keep], on="mol_id", how="left")

        # Public DB / validation support: a bonus and review field, never a hard filter.
        support_map: dict[str, float] = {}
        support_notes_map: dict[str, Any] = {}
        global_support = 0.0
        for art in list(manifest.iter_artifacts("method_validation")) + list(manifest.iter_artifacts("calibration")):
            data = art.data or {}
            support_map.update({str(k): float(v) for k, v in (data.get("db_calibration_support_by_mol") or {}).items()})
            support_notes_map.update(data.get("support_notes_by_mol") or {})
            if data.get("calibration_support_score") is not None:
                support_value = pd.to_numeric(
                    pd.Series([data.get("calibration_support_score")]), errors="coerce"
                ).iloc[0]
                if pd.notna(support_value):
                    global_support = max(global_support, float(support_value))
        df["db_calibration_support"] = df["mol_id"].map(lambda x: float(support_map.get(str(x), global_support))) if "mol_id" in df.columns else global_support
        df["db_calibration_notes"] = df["mol_id"].map(lambda x: ";".join(support_notes_map.get(str(x), [])) if isinstance(support_notes_map.get(str(x), []), list) else str(support_notes_map.get(str(x), ""))) if "mol_id" in df.columns else ""

        if "delta_G_assoc_pressure_corrected_kcal_mol" not in df.columns:
            df["delta_G_assoc_pressure_corrected_kcal_mol"] = df.get("delta_G_assoc_kcal_mol", pd.Series([np.nan] * len(df)))
        for required in ["delta_G_assoc_kcal_mol", "delta_G_ionpair_kcal_mol", "delta_G_act_kcal_mol"]:
            if required not in df.columns:
                df[required] = np.nan
        df["delta_G_act_effective_kcal_mol"] = df.get("delta_G_act_pressure_corrected_kcal_mol", df["delta_G_act_kcal_mol"])
        df["delta_G_act_effective_kcal_mol"] = pd.to_numeric(df["delta_G_act_effective_kcal_mol"], errors="coerce").fillna(pd.to_numeric(df["delta_G_act_kcal_mol"], errors="coerce"))
        if "quality_tier" not in df.columns:
            df["quality_tier"] = "Q0"
        df["quality_numeric"] = df["quality_tier"].map(tier_numeric).astype(float)

        # Simple science gate: development/dummy values are allowed in screening
        # outputs, but are clearly excluded from scientific/production rankings.
        for col in ["main_values_are_dummy", "main_values_are_fallback", "production_thermo_ready", "real_irc_executed"]:
            if col not in df.columns:
                df[col] = False
            df[col] = df[col].fillna(False).astype(bool)
        if "has_scientific_dft" not in df.columns:
            df["has_scientific_dft"] = df["quality_numeric"] >= tier_numeric("Q2")
        df["main_values_are_dummy"] = df["main_values_are_dummy"] | (df["quality_tier"].astype(str).isin(["Q0", "Q1"]) & ~df["has_scientific_dft"].astype(bool))
        if "confidence_score_raw" not in df.columns:
            df["confidence_score_raw"] = df.get("confidence_score", pd.Series([np.nan] * len(df)))
        df["confidence_score"] = df.apply(
            lambda r: cap_confidence(r.get("confidence_score_raw"), r.get("quality_tier"), bool(r.get("main_values_are_dummy") or r.get("main_values_are_fallback"))),
            axis=1,
        )
        if "scientific_rank_eligible" not in df.columns:
            df["scientific_rank_eligible"] = df.apply(lambda r: scientific_rank_eligible(r.get("quality_tier"), bool(r.get("main_values_are_dummy") or r.get("main_values_are_fallback"))), axis=1)
        else:
            df["scientific_rank_eligible"] = df["scientific_rank_eligible"].fillna(False).astype(bool)
        if "production_rank_eligible" not in df.columns:
            df["production_rank_eligible"] = df.apply(lambda r: production_rank_eligible(r.get("quality_tier"), bool(r.get("main_values_are_dummy") or r.get("main_values_are_fallback")), bool(r.get("production_thermo_ready")), bool(r.get("real_irc_executed"))), axis=1)
        else:
            df["production_rank_eligible"] = df["production_rank_eligible"].fillna(False).astype(bool)
        if "science_gate_reason" not in df.columns:
            df["science_gate_reason"] = df.apply(
                lambda r: "production_ready" if r.get("production_rank_eligible") else ("scientific_dft_ready" if r.get("scientific_rank_eligible") else "screening_only"),
                axis=1,
            )

        assoc_for_risk = -pd.to_numeric(df.get("delta_G_assoc_pressure_corrected_kcal_mol", df.get("delta_G_assoc_kcal_mol")), errors="coerce").fillna(0.0)
        lowfreq_for_risk = pd.to_numeric(df.get("low_frequency_count_reactant", pd.Series([0] * len(df))), errors="coerce").fillna(0.0)
        df["cluster_growth_score"] = (pd.to_numeric(df["hf_n"], errors="coerce").fillna(1.0) * 0.20) + 0.05 * assoc_for_risk + 0.03 * lowfreq_for_risk

        public_penalty = df.get("public_data", pd.Series([{}] * len(df))).map(lambda x: process_penalty_from_public_data(x if isinstance(x, dict) else {}))
        identity_penalty = df.get("identity", pd.Series([{}] * len(df))).map(lambda x: 0.4 if isinstance(x, dict) and x.get("identity_conflict") else 0.0)
        df["process_penalty"] = pd.to_numeric(public_penalty, errors="coerce").fillna(0.0) + pd.to_numeric(identity_penalty, errors="coerce").fillna(0.0)
        df["irreversible_trapping_flag"] = pd.to_numeric(df["delta_G_assoc_kcal_mol"], errors="coerce").fillna(0.0) < -18.0
        df.loc[df["irreversible_trapping_flag"], "process_penalty"] += 0.5
        df["gas_process_feasibility"] = df.get("public_data", pd.Series([{}] * len(df))).map(lambda x: _dict_value(x, "gas_process", "gas_process_feasibility") or _dict_value(x, "comptox", "gas_process_feasibility"))
        df["ehs_review_flag"] = df.get("public_data", pd.Series([{}] * len(df))).map(lambda x: _dict_value(x, "gas_process", "ehs_review_flag") or _dict_value(x, "comptox", "ehs_review_flag"))

        df["scavenger_score"] = (
            0.25 * _z(-pd.to_numeric(df["delta_G_assoc_pressure_corrected_kcal_mol"], errors="coerce"))
            + 0.25 * _z(-pd.to_numeric(df["delta_G_ionpair_kcal_mol"], errors="coerce"))
            + 0.20 * _z(-pd.to_numeric(df["delta_G_act_effective_kcal_mol"], errors="coerce"))
            - 0.15 * _z(df["cluster_growth_score"])
            - 0.10 * _z(df["process_penalty"])
            + 0.03 * _z(df["quality_numeric"])
            + 0.02 * _z(df["db_calibration_support"])
        )
        df["activation_score"] = (
            0.25 * _z(df.get("delta_r_HF_A", pd.Series([np.nan] * len(df))))
            + 0.20 * _z(-pd.to_numeric(df.get("delta_nu_HF_cm1", pd.Series([np.nan] * len(df))), errors="coerce"))
            + 0.25 * _z(-pd.to_numeric(df["delta_G_act_effective_kcal_mol"], errors="coerce"))
            + 0.15 * _z(_moderate_binding_score(df["delta_G_assoc_kcal_mol"]))
            - 0.10 * _z(df["process_penalty"])
            + 0.03 * _z(df["quality_numeric"])
            + 0.02 * _z(df["db_calibration_support"])
        )

        def _fmt(v, digits=2):
            try:
                return f"{float(v):.{digits}f}" if pd.notna(v) else "NA"
            except (TypeError, ValueError):
                return "NA"

        df["key_evidence"] = df.apply(
            lambda r: (
                f"Q={r.get('quality_tier')}; ΔGassoc(p)={_fmt(r.get('delta_G_assoc_pressure_corrected_kcal_mol'))}; "
                f"ΔG‡={_fmt(r.get('delta_G_act_effective_kcal_mol'))}; "
                f"k={_fmt(r.get('k_corrected_s-1'), 3)}; "
                f"ΔrHF={_fmt(r.get('delta_r_HF_A'), 3)}; "
                f"DBsupport={_fmt(r.get('db_calibration_support'), 2)}"
            ),
            axis=1,
        )

        df["overall_score"] = df[["scavenger_score", "activation_score"]].max(axis=1)
        scav = df.sort_values(["scavenger_score", "quality_numeric"], ascending=False).copy()
        act = df.sort_values(["activation_score", "quality_numeric"], ascending=False).copy()
        cluster = df.sort_values(["cluster_growth_score", "quality_numeric"], ascending=False).copy()
        screening = df.sort_values(["overall_score", "quality_numeric"], ascending=False).copy()
        scientific = screening[screening["scientific_rank_eligible"]].copy()
        production = screening[screening["production_rank_eligible"]].copy()
        scav["rank"] = range(1, len(scav) + 1)
        act["rank"] = range(1, len(act) + 1)
        cluster["risk_rank"] = range(1, len(cluster) + 1)
        screening["rank"] = range(1, len(screening) + 1)
        scientific["rank"] = range(1, len(scientific) + 1)
        production["rank"] = range(1, len(production) + 1)
        best = df.sort_values(["production_rank_eligible", "scientific_rank_eligible", "overall_score", "quality_numeric"], ascending=False).groupby("mol_id", as_index=False).first()
        candidate_summary = best[[c for c in [
            "mol_id", "name", "canonical_smiles", "inchikey", "site_id", "site_type", "hf_n",
            "quality_tier", "confidence_score", "confidence_score_raw", "scientific_rank_eligible", "production_rank_eligible",
            "main_values_are_dummy", "main_values_are_fallback", "science_gate_reason", "scavenger_score", "activation_score",
            "delta_G_assoc_pressure_corrected_kcal_mol", "delta_G_act_kcal_mol", "k_corrected_s-1",
            "K_assoc_process_adjusted", "db_calibration_support", "db_calibration_notes", "process_penalty",
            "gas_process_feasibility", "ehs_review_flag", "key_evidence"
        ] if c in best.columns]].copy()

        def _science_action(row):
            if row.get("production_rank_eligible"):
                return "review_production_candidate"
            if row.get("scientific_rank_eligible"):
                return "run_ts_irc_or_high_level_sp"
            if row.get("main_values_are_dummy") or row.get("quality_tier") in {"Q0", "Q1"}:
                return "run_real_dft_minima"
            return "review_science_gate"

        def _process_action(row):
            if row.get("process_penalty", 0.0) and float(row.get("process_penalty", 0.0)) >= 0.5:
                return "manual_process_EHS_review"
            return "check_vapor_pressure_and_tool_compatibility"

        def _ops_action(row):
            if row.get("production_rank_eligible"):
                return "archive_or_compare_high_level_result"
            if row.get("scientific_rank_eligible"):
                return "submit_ts_or_sp_jobs"
            return "submit_dft_screening_jobs"

        candidate_summary["science_next_action"] = candidate_summary.apply(_science_action, axis=1)
        candidate_summary["process_next_action"] = candidate_summary.apply(_process_action, axis=1)
        candidate_summary["ops_next_action"] = candidate_summary.apply(_ops_action, axis=1)
        candidate_summary["recommended_next_action"] = candidate_summary["science_next_action"]

        scav_path = out_dir / "rank_scavenger.csv"
        act_path = out_dir / "rank_activation.csv"
        screening_path = out_dir / "rank_screening.csv"
        scientific_path = out_dir / "rank_scientific.csv"
        production_path = out_dir / "rank_production.csv"
        summary_path = out_dir / "candidate_summary.csv"
        cluster_path = out_dir / "cluster_risk.csv"
        failure_path = out_dir / "failure_report.csv"
        metadata_path = out_dir / "ranking_summary.json"
        parquet_path = out_dir / "reaction_results.parquet"
        scav.to_csv(scav_path, index=False)
        act.to_csv(act_path, index=False)
        screening.to_csv(screening_path, index=False)
        scientific.to_csv(scientific_path, index=False)
        production.to_csv(production_path, index=False)
        candidate_summary.to_csv(summary_path, index=False)
        cluster.to_csv(cluster_path, index=False)
        failure_columns = ["artifact_id", "artifact_type", "category", "reason", "recommended_fallback"]
        failure_rows = [
            {
                "artifact_id": f.artifact_id,
                "artifact_type": f.artifact_type,
                "category": f.status.category,
                "reason": f.status.reason,
                "recommended_fallback": f.status.recommended_fallback,
            }
            for f in failures
        ]
        pd.DataFrame(failure_rows, columns=failure_columns).to_csv(failure_path, index=False)
        ranking_summary = {
            "mode": context.global_config.get("mode", "development"),
            "n_rows": len(df),
            "n_screening_rows": len(screening),
            "n_scientific_rows": len(scientific),
            "n_production_rows": len(production),
            "contains_dummy_or_fallback": bool((df["main_values_are_dummy"] | df["main_values_are_fallback"]).any()),
            "warning": "screening rankings may include dummy/fallback values; use rank_production.csv for production decisions",
        }
        metadata_path.write_text(json.dumps(ranking_summary, ensure_ascii=False, indent=2), encoding="utf-8")
        try:
            df.to_parquet(parquet_path, index=False)
            table_path = parquet_path
            table_key = "parquet"
        except (ImportError, OSError, ValueError):
            table_path = out_dir / "reaction_results.csv"
            table_key = "csv"
            df.to_csv(table_path, index=False)
        for artifact_id, ranking_type, path, rows in [
            ("rank_scavenger", "scavenger", scav_path, len(scav)),
            ("rank_activation", "activation", act_path, len(act)),
            ("rank_screening", "screening", screening_path, len(screening)),
            ("rank_scientific", "scientific", scientific_path, len(scientific)),
            ("rank_production", "production", production_path, len(production)),
            ("cluster_risk", "cluster_risk", cluster_path, len(cluster)),
        ]:
            out.add_artifact(Artifact(artifact_id=artifact_id, artifact_type="ranking", paths={"csv": str(path)}, data={"ranking_type": ranking_type, "n_rows": int(rows)}))
        out.add_artifact(Artifact(artifact_id="ranking_summary", artifact_type="table", paths={"json": str(metadata_path)}, data=ranking_summary))
        out.add_artifact(Artifact(artifact_id="candidate_summary", artifact_type="table", paths={"csv": str(summary_path)}, data={"n_rows": len(candidate_summary)}))
        out.add_artifact(Artifact(artifact_id="failure_report", artifact_type="table", paths={"csv": str(failure_path)}, data={"n_rows": len(failures)}))
        out.add_artifact(Artifact(artifact_id="reaction_results", artifact_type="table", paths={table_key: str(table_path)}, data={"n_rows": len(df)}))
        return out

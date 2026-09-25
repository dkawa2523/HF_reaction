"""Evidence-first ranking for chemistry-neutral reaction workflows."""

from __future__ import annotations

import json
from math import log
from typing import Any

import pandas as pd

from hfauto.chemistry.populations import R_KCAL_MOL_K
from hfauto.chemistry.stoichiometry import reaction_side_terms
from hfauto.core.io import ensure_dir
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import Stage, StageContext


def _reaction_metadata(manifest: Manifest) -> dict[str, dict[str, Any]]:
    metadata: dict[str, dict[str, Any]] = {}
    for artifact_type in ("reaction", "reaction_validated"):
        for artifact in manifest.latest_artifacts(artifact_type):
            reaction_id = str(artifact.data.get("reaction_id") or artifact.artifact_id)
            metadata[reaction_id] = {
                key: artifact.data.get(key)
                for key in (
                    "candidate_id",
                    "comparison_group",
                    "mechanism_family",
                    "reaction_type",
                )
            }
    return metadata


def _method_uncertainty(manifest: Manifest) -> dict[str, dict[str, Any]]:
    return {
        str(artifact.data["reaction_id"]): dict(artifact.data)
        for artifact in manifest.latest_artifacts("method_uncertainty")
        if artifact.data.get("reaction_id")
    }


def _reaction_artifacts(manifest: Manifest) -> dict[str, Artifact]:
    reactions: dict[str, Artifact] = {}
    for artifact_type in ("reaction", "reaction_validated", "reaction_path_validated"):
        for artifact in manifest.latest_artifacts(artifact_type):
            reaction_id = str(artifact.data.get("reaction_id") or artifact.artifact_id)
            previous = reactions.get(reaction_id)
            if previous is None:
                reactions[reaction_id] = artifact
                continue
            merged = artifact.model_copy(deep=True)
            merged.data = {**previous.data, **artifact.data}
            reactions[reaction_id] = merged
    return reactions


def _basin_population_evidence(
    manifest: Manifest,
    reaction: Artifact | None,
    temperature_K: float,
) -> dict[str, Any]:
    if reaction is None:
        return {"population_evidence_available": False}
    populations = {
        (str(artifact.data.get("species_id")), float(artifact.data.get("T_K"))): artifact
        for artifact in manifest.latest_artifacts("basin_population")
        if artifact.data.get("species_id") and artifact.data.get("T_K") is not None
    }
    terms = reaction_side_terms(reaction, "reactants")
    selected = [
        (term, populations.get((term.species_id, float(temperature_K))))
        for term in terms
    ]
    if not selected or any(artifact is None for _term, artifact in selected):
        return {"population_evidence_available": False}
    population = 1.0
    for term, artifact in selected:
        assert artifact is not None
        population *= float(artifact.data["conditional_population"]) ** float(
            term.coefficient
        )
    adjustment = (
        -R_KCAL_MOL_K * float(temperature_K) * log(population)
        if population > 0.0
        else None
    )
    artifacts = [artifact for _term, artifact in selected if artifact is not None]
    return {
        "population_evidence_available": adjustment is not None,
        "reactant_conditional_population": population,
        "population_adjustment_kcal_mol": adjustment,
        "population_ensemble_complete": all(
            artifact.data.get("ensemble_complete") is True for artifact in artifacts
        ),
        "population_scope": "observed_frequency_validated_basins_only",
    }


class ReactionRankStage(Stage):
    """Rank like-for-like reactions by one declared metric and evidence tier."""

    name = "reaction-rank"

    def __init__(self, *, output_stage_name: str | None = None):
        self.output_stage_name = output_stage_name or self.name

    def run(
        self,
        manifest: Manifest | None,
        config: dict[str, Any],
        context: StageContext,
    ) -> Manifest:
        if manifest is None:
            raise ValueError("reaction-rank requires an input manifest")
        out_dir = ensure_dir(context.out_dir)
        out = manifest.carry_forward(self.output_stage_name)
        rows = [
            dict(artifact.data)
            for artifact in manifest.latest_artifacts("thermo")
            if artifact.data.get("stoichiometry")
        ]
        if not rows:
            out.add_artifact(
                Artifact.failure(
                    "reaction_rank_failed",
                    "ranking",
                    "no stoichiometric reaction thermochemistry records",
                    category="missing_input",
                )
            )
            return out

        metadata = _reaction_metadata(manifest)
        uncertainty = _method_uncertainty(manifest)
        reactions = _reaction_artifacts(manifest)
        for row in rows:
            reaction_id = str(row.get("reaction_id"))
            for key, value in metadata.get(reaction_id, {}).items():
                row.setdefault(key, value)
            for key, value in uncertainty.get(reaction_id, {}).items():
                if key != "reaction_id":
                    row[f"method_{key}"] = value
            row.update(
                _basin_population_evidence(
                    manifest,
                    reactions.get(reaction_id),
                    float(row.get("T_K", 298.15)),
                )
            )
        frame = pd.DataFrame(rows)
        metric = str(
            config.get("metric", "delta_G_activation_standard_kcal_mol")
        )
        if metric not in frame.columns:
            fallback = "delta_G_reaction_standard_kcal_mol"
            if fallback not in frame.columns:
                raise ValueError(f"reaction-rank metric is unavailable: {metric!r}")
            metric = fallback
        objective = str(config.get("objective", "lower")).lower()
        if objective not in {"lower", "higher"}:
            raise ValueError("reaction-rank objective must be 'lower' or 'higher'")
        default_group = frame.get(
            "mechanism_family", pd.Series(["unspecified"] * len(frame))
        ).fillna("unspecified")
        if "comparison_group" not in frame.columns:
            frame["comparison_group"] = default_group
        else:
            frame["comparison_group"] = frame["comparison_group"].fillna(default_group)
        frame["ranking_metric"] = metric
        frame["ranking_objective"] = objective
        frame["ranking_point_value"] = pd.to_numeric(frame[metric], errors="coerce")
        uncertainty_policy = str(
            config.get("uncertainty_policy", "report_only")
        ).lower()
        if uncertainty_policy not in {
            "report_only",
            "conservative_bound",
            "require_characterized",
        }:
            raise ValueError(
                "reaction-rank uncertainty_policy must be 'report_only', "
                "'conservative_bound', or 'require_characterized'"
            )
        uncertainty_column = str(
            config.get(
                "uncertainty_column",
                "method_activation_max_abs_shift_from_reference_kcal_mol",
            )
        )
        fallback_column = "method_activation_half_range_kcal_mol"
        if uncertainty_column in frame.columns:
            uncertainty_values = pd.to_numeric(
                frame[uncertainty_column], errors="coerce"
            )
        else:
            uncertainty_values = pd.Series(float("nan"), index=frame.index)
        if fallback_column in frame.columns:
            uncertainty_values = uncertainty_values.fillna(
                pd.to_numeric(frame[fallback_column], errors="coerce")
            )
        frame["ranking_uncertainty_kcal_mol"] = uncertainty_values
        frame["uncertainty_policy"] = uncertainty_policy
        if uncertainty_policy == "conservative_bound":
            direction = 1.0 if objective == "lower" else -1.0
            frame["ranking_value"] = frame["ranking_point_value"] + direction * frame[
                "ranking_uncertainty_kcal_mol"
            ].fillna(0.0)
        else:
            frame["ranking_value"] = frame["ranking_point_value"]
        population_policy = str(
            config.get("population_policy", "report_only")
        ).lower()
        if population_policy not in {
            "report_only",
            "conditional_adjustment",
            "require_complete",
        }:
            raise ValueError(
                "reaction-rank population_policy must be 'report_only', "
                "'conditional_adjustment', or 'require_complete'"
            )
        frame["population_policy"] = population_policy
        if population_policy == "conditional_adjustment":
            if "activation" not in metric.lower():
                raise ValueError(
                    "conditional population adjustment is only defined for an "
                    "activation free-energy metric"
                )
            available = frame.get(
                "population_evidence_available",
                pd.Series(False, index=frame.index),
            ).fillna(False).astype(bool)
            adjustment = pd.to_numeric(
                frame.get(
                    "population_adjustment_kcal_mol",
                    pd.Series(float("nan"), index=frame.index),
                ),
                errors="coerce",
            )
            frame["ranking_value"] = frame["ranking_value"] + adjustment
            frame["population_adjustment_applied"] = available & adjustment.notna()
        else:
            frame["population_adjustment_applied"] = False
        evidence_policy = str(config.get("evidence_policy", "scientific")).lower()
        eligibility_columns = {
            "all": None,
            "scientific": "scientific_rank_eligible",
            "production": "production_rank_eligible",
        }
        if evidence_policy not in eligibility_columns:
            raise ValueError(
                "reaction-rank evidence_policy must be 'all', 'scientific', or 'production'"
            )
        eligibility_column = eligibility_columns[evidence_policy]
        if eligibility_column is None:
            evidence_eligible = pd.Series(True, index=frame.index)
        elif eligibility_column in frame.columns:
            evidence_eligible = frame[eligibility_column].fillna(False).astype(bool)
        else:
            evidence_eligible = pd.Series(False, index=frame.index)
        frame["evidence_policy"] = evidence_policy
        frame["evidence_eligible"] = evidence_eligible
        if uncertainty_policy == "require_characterized":
            characterized_column = "method_method_uncertainty_characterized"
            characterized = (
                frame[characterized_column].fillna(False).astype(bool)
                if characterized_column in frame.columns
                else pd.Series(False, index=frame.index)
            )
            frame["evidence_eligible"] &= characterized
        if population_policy == "conditional_adjustment":
            frame["evidence_eligible"] &= frame[
                "population_adjustment_applied"
            ].astype(bool)
        elif population_policy == "require_complete":
            complete = frame.get(
                "population_ensemble_complete",
                pd.Series(False, index=frame.index),
            ).fillna(False).astype(bool)
            frame["evidence_eligible"] &= complete
        frame["rankable"] = frame["ranking_value"].notna() & frame["evidence_eligible"]
        frame["rank_within_comparison_group"] = frame["ranking_value"].where(
            frame["rankable"]
        ).groupby(
            [frame["comparison_group"], frame["T_K"]], dropna=False
        ).rank(
            method="min", ascending=objective == "lower", na_option="bottom"
        )
        frame = frame.sort_values(
            ["comparison_group", "T_K", "rankable", "rank_within_comparison_group"],
            ascending=[True, True, False, True],
        )

        csv_path = out_dir / "reaction_ranking.csv"
        json_path = out_dir / "reaction_ranking_summary.json"
        frame.to_csv(csv_path, index=False)
        summary = {
            "metric": metric,
            "objective": objective,
            "evidence_policy": evidence_policy,
            "uncertainty_policy": uncertainty_policy,
            "population_policy": population_policy,
            "uncertainty_column": uncertainty_column,
            "n_rows": len(frame),
            "n_rankable": int(frame["rankable"].sum()),
            "comparison_groups": sorted(
                str(value) for value in frame["comparison_group"].dropna().unique()
            ),
            "warning": (
                "Ranks are valid only within comparison_group and temperature; "
                "ineligible or dummy/fallback rows are retained for audit but are not ranked."
            ),
        }
        json_path.write_text(
            json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        out.add_artifact(
            Artifact(
                artifact_id="reaction_ranking",
                artifact_type="ranking",
                paths={"csv": str(csv_path), "summary_json": str(json_path)},
                data=summary,
                qc={"like_for_like_grouping_required": True},
            )
        )
        return out

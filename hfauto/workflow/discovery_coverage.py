"""Evidence accounting for bounded reaction-discovery campaigns.

Coverage is an audit of what a finite search budget actually exercised. It is
not an estimator of undiscovered reactions and never turns a negative bounded
search into a proof that no reaction exists.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Any

from hfauto.core.artifacts import preferred_species_by_id
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest


def _successful(manifest: Manifest, artifact_type: str) -> list[Artifact]:
    return [
        artifact
        for artifact in manifest.latest_artifacts(artifact_type)
        if artifact.status.status == "success"
    ]


def _composition_key(species: Artifact | None) -> str:
    if species is None:
        return "unknown"
    counts: dict[str, int] = defaultdict(int)
    for component in species.data.get("components", []) or []:
        for element, count in (component.get("element_counts") or {}).items():
            counts[str(element)] += int(count)
    if not counts:
        for element, count in (species.data.get("element_counts") or {}).items():
            counts[str(element)] += int(count)
    return "".join(
        f"{element}{counts[element]}" for element in sorted(counts)
    ) or "unknown"


def _scientific_conclusion(row: dict[str, Any]) -> str:
    if row["irc_validated_reaction_count"]:
        return "validated_elementary_connection_observed"
    if row["frequency_validated_ts_count"]:
        return "frequency_validated_ts_without_irc"
    if row["promoted_reaction_count"]:
        return "distinct_dft_basins_without_validated_ts"
    if row["candidate_count"]:
        return "candidates_without_distinct_dft_basins"
    if row["discovery_attempt_count"]:
        return "no_distinct_product_observed_within_budget"
    return "trials_not_executed"


def build_discovery_coverage(
    manifest: Manifest,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Measure trial-to-IRC attrition for every searched source geometry."""

    trials = _successful(manifest, "reaction_trial")
    attempts = manifest.latest_artifacts("reaction_discovery_attempt")
    candidates = _successful(manifest, "reaction_candidate")
    reactions = _successful(manifest, "reaction")
    validated_ts = _successful(manifest, "reaction_validated")
    validated_paths = [
        *_successful(manifest, "reaction_path_validated"),
        *[
            artifact
            for artifact in _successful(manifest, "irc")
            if artifact.qc.get("irc_validated") is True
        ],
    ]
    species = preferred_species_by_id(
        manifest,
        artifact_types=("species", "species_preopt", "species_optimized"),
    )

    attempts_by_trial: dict[str, list[Artifact]] = defaultdict(list)
    for attempt in attempts:
        trial_id = attempt.data.get("trial_id")
        if trial_id:
            attempts_by_trial[str(trial_id)].append(attempt)
    candidates_by_trial: dict[str, list[Artifact]] = defaultdict(list)
    for candidate in candidates:
        trial_id = candidate.data.get("trial_id")
        if trial_id:
            candidates_by_trial[str(trial_id)].append(candidate)

    reactions_by_candidate: dict[str, list[Artifact]] = defaultdict(list)
    for reaction in [*reactions, *validated_ts]:
        candidate_id = reaction.data.get("candidate_id")
        if candidate_id:
            reactions_by_candidate[str(candidate_id)].append(reaction)
    ts_reactions = {
        str(item.data.get("reaction_id") or item.artifact_id)
        for item in validated_ts
        if item.qc.get("ts_validated_by_frequency") is True
    }
    irc_reactions = {
        str(item.data.get("reaction_id") or item.artifact_id)
        for item in validated_paths
        if item.qc.get("irc_validated") is True
    }

    registry_items = manifest.latest_artifacts("minimum_registry")
    species_to_basin = (
        dict(registry_items[-1].data.get("species_to_basin") or {})
        if registry_items
        else {}
    )
    trial_summary_items = [
        item
        for item in manifest.latest_artifacts("table")
        if item.artifact_id == "reaction_trial_summary"
    ]
    trial_summary = trial_summary_items[-1].data if trial_summary_items else {}
    per_source_budget = trial_summary.get("max_trials_per_species")

    trials_by_source: dict[str, list[Artifact]] = defaultdict(list)
    for trial in trials:
        source_id = trial.data.get("source_species_id")
        if source_id:
            trials_by_source[str(source_id)].append(trial)

    rows: list[dict[str, Any]] = []
    for source_id in sorted(trials_by_source):
        source_trials = trials_by_source[source_id]
        source_attempts = [
            attempt
            for trial in source_trials
            for attempt in attempts_by_trial.get(
                str(trial.data.get("trial_id") or trial.artifact_id), []
            )
        ]
        source_candidates = [
            candidate
            for trial in source_trials
            for candidate in candidates_by_trial.get(
                str(trial.data.get("trial_id") or trial.artifact_id), []
            )
        ]
        source_reactions = [
            reaction
            for candidate in source_candidates
            for reaction in reactions_by_candidate.get(
                str(candidate.data.get("candidate_id") or candidate.artifact_id),
                [],
            )
        ]
        reaction_ids = {
            str(reaction.data.get("reaction_id") or reaction.artifact_id)
            for reaction in source_reactions
        }
        low_level_candidates = [
            candidate
            for candidate in source_candidates
            if candidate.data.get("endpoint_evidence_level", "low_level")
            == "low_level"
        ]
        dft_ensemble_candidates = [
            candidate
            for candidate in source_candidates
            if candidate.data.get("endpoint_evidence_level")
            == "dft_minimum_ensemble"
        ]
        expected_attempts = sum(
            min(
                int(trial.data.get("max_attempts", 1)),
                len(trial.data.get("driver_order") or []),
            )
            for trial in source_trials
        )
        requested_drivers = sorted(
            {
                str(driver)
                for trial in source_trials
                for driver in trial.data.get("driver_order", []) or []
            }
        )
        attempted_drivers = sorted(
            {
                str(attempt.data.get("driver"))
                for attempt in source_attempts
                if attempt.data.get("driver")
            }
        )
        basins = {
            str(basin_id)
            for candidate in source_candidates
            for species_id in (
                candidate.data.get("reactant_species_id"),
                candidate.data.get("product_species_id"),
            )
            if species_id
            and (basin_id := species_to_basin.get(str(species_id)))
        }
        source = species.get(source_id)
        row: dict[str, Any] = {
            "source_species_id": source_id,
            "source_artifact_id": source.artifact_id if source else None,
            "composition_key": _composition_key(source),
            "component_count": len(source.data.get("components") or []) if source else 0,
            "trial_count": len(source_trials),
            "expected_attempt_count": expected_attempts,
            "discovery_attempt_count": len(source_attempts),
            "attempt_completion_fraction": (
                min(1.0, len(source_attempts) / expected_attempts)
                if expected_attempts
                else 0.0
            ),
            "requested_drivers": ";".join(requested_drivers),
            "attempted_drivers": ";".join(attempted_drivers),
            "driver_coverage_fraction": (
                len(set(attempted_drivers) & set(requested_drivers))
                / len(requested_drivers)
                if requested_drivers
                else 0.0
            ),
            "candidate_count": len(source_candidates),
            "low_level_candidate_count": len(low_level_candidates),
            "dft_minimum_ensemble_candidate_count": len(dft_ensemble_candidates),
            "candidate_yield_per_attempt": (
                len(low_level_candidates) / len(source_attempts)
                if source_attempts
                else 0.0
            ),
            "registered_basin_count_touched": len(basins),
            "promoted_reaction_count": len(source_reactions),
            "frequency_validated_ts_count": len(reaction_ids & ts_reactions),
            "irc_validated_reaction_count": len(reaction_ids & irc_reactions),
            "trial_budget_saturated": bool(
                per_source_budget is not None
                and len(source_trials) >= int(per_source_budget)
            ),
            "finite_search_only": True,
            "exhaustive_claim_allowed": False,
        }
        row["scientific_conclusion"] = _scientific_conclusion(row)
        rows.append(row)

    source_count_declared = int(
        trial_summary.get("source_species_count", len(rows))
    )
    summary = {
        "source_species_count_declared": source_count_declared,
        "source_species_with_trials": len(rows),
        "source_species_without_trials": max(0, source_count_declared - len(rows)),
        "trial_count": len(trials),
        "discovery_attempt_count": len(attempts),
        "candidate_count": len(candidates),
        "low_level_candidate_count": sum(
            candidate.data.get("endpoint_evidence_level", "low_level")
            == "low_level"
            for candidate in candidates
        ),
        "dft_minimum_ensemble_candidate_count": sum(
            candidate.data.get("endpoint_evidence_level")
            == "dft_minimum_ensemble"
            for candidate in candidates
        ),
        "minimum_basin_count": len(_successful(manifest, "minimum_basin")),
        "promoted_reaction_count": len(
            {
                str(reaction.data.get("reaction_id") or reaction.artifact_id)
                for reaction in [*reactions, *validated_ts]
                if reaction.data.get("candidate_id")
                and (reaction.data.get("basin_assessment") or {}).get(
                    "distinct_basin"
                )
                is True
            }
        ),
        "frequency_validated_ts_count": len(ts_reactions),
        "irc_validated_reaction_count": len(irc_reactions),
        "trial_budget_saturated": bool(
            trial_summary.get("max_total_trials") is not None
            and len(trials) >= int(trial_summary["max_total_trials"])
        ),
        "finite_search_only": True,
        "exhaustive_claim_allowed": False,
        "interpretation": (
            "Counts describe observed attrition under declared finite budgets; "
            "zero candidates is not proof that no reaction exists."
        ),
    }
    return rows, summary


def build_discovery_saturation(
    campaigns: list[tuple[str, Manifest]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Compare observed yield as declared campaign budgets increase.

    A flat observed count is reported as a plateau, never as proof that the
    underlying reaction network is exhausted.
    """

    series: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for campaign_id, manifest in campaigns:
        rows, _summary = build_discovery_coverage(manifest)
        for row in rows:
            key = (
                str(row["source_species_id"]),
                str(row["composition_key"]),
            )
            series[key].append(
                {
                    "campaign_id": str(campaign_id),
                    **row,
                }
            )

    output: list[dict[str, Any]] = []
    plateau_series = 0
    for key in sorted(series):
        ordered = sorted(
            series[key],
            key=lambda row: (
                int(row["discovery_attempt_count"]),
                int(row["trial_count"]),
                str(row["campaign_id"]),
            ),
        )
        previous: dict[str, Any] | None = None
        series_has_plateau = False
        for level, row in enumerate(ordered, start=1):
            compared = dict(row)
            compared["budget_level"] = level
            compared["previous_campaign_id"] = (
                previous["campaign_id"] if previous else None
            )
            budget_expanded = bool(
                previous
                and (
                    int(row["discovery_attempt_count"])
                    > int(previous["discovery_attempt_count"])
                    or int(row["trial_count"]) > int(previous["trial_count"])
                )
            )
            compared["budget_expanded"] = budget_expanded
            for metric in (
                "trial_count",
                "discovery_attempt_count",
                "candidate_count",
                "registered_basin_count_touched",
                "promoted_reaction_count",
                "frequency_validated_ts_count",
                "irc_validated_reaction_count",
            ):
                compared[f"delta_{metric}"] = (
                    int(row[metric]) - int(previous[metric])
                    if previous
                    else None
                )
            plateau = bool(
                budget_expanded
                and all(
                    compared[f"delta_{metric}"] == 0
                    for metric in (
                        "candidate_count",
                        "registered_basin_count_touched",
                        "promoted_reaction_count",
                        "frequency_validated_ts_count",
                        "irc_validated_reaction_count",
                    )
                )
            )
            compared.update(
                {
                    "observed_yield_plateau": plateau,
                    "finite_campaign_comparison_only": True,
                    "exhaustive_claim_allowed": False,
                }
            )
            series_has_plateau = series_has_plateau or plateau
            output.append(compared)
            previous = row
        plateau_series += int(series_has_plateau)

    summary = {
        "campaign_count": len(campaigns),
        "source_series_count": len(series),
        "source_series_with_observed_plateau": plateau_series,
        "comparison_row_count": len(output),
        "finite_campaign_comparison_only": True,
        "discovery_saturation_established": False,
        "exhaustive_claim_allowed": False,
        "interpretation": (
            "A plateau means only that the observed counts did not increase "
            "between the supplied finite budgets."
        ),
    }
    return output, summary

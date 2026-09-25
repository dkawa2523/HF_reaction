from hfauto.core.schemas.artifact import Artifact, ArtifactStatus
from hfauto.core.schemas.manifest import Manifest
from hfauto.workflow.discovery_coverage import (
    build_discovery_coverage,
    build_discovery_saturation,
)


def test_discovery_coverage_reports_budgeted_attrition_without_exhaustive_claim() -> None:
    manifest = Manifest.new(run_id="coverage", stage="reaction-plan")
    manifest.extend(
        [
            Artifact(
                artifact_id="source",
                artifact_type="species_preopt",
                data={
                    "species_id": "source",
                    "components": [
                        {"component_id": "a", "element_counts": {"H": 2}},
                        {"component_id": "b", "element_counts": {"F": 1}},
                    ],
                },
            ),
            Artifact(
                artifact_id="trial",
                artifact_type="reaction_trial",
                data={
                    "trial_id": "trial",
                    "source_species_id": "source",
                    "driver_order": ["nt2", "afir"],
                    "max_attempts": 2,
                },
            ),
            Artifact(
                artifact_id="attempt",
                artifact_type="reaction_discovery_attempt",
                data={"trial_id": "trial", "driver": "nt2"},
                status=ArtifactStatus(status="partial"),
            ),
            Artifact(
                artifact_id="reaction_trial_summary",
                artifact_type="table",
                data={
                    "source_species_count": 1,
                    "max_trials_per_species": 1,
                    "max_total_trials": 1,
                },
            ),
        ]
    )

    rows, summary = build_discovery_coverage(manifest)

    assert rows[0]["composition_key"] == "F1H2"
    assert rows[0]["attempt_completion_fraction"] == 0.5
    assert rows[0]["driver_coverage_fraction"] == 0.5
    assert rows[0]["scientific_conclusion"] == (
        "no_distinct_product_observed_within_budget"
    )
    assert rows[0]["trial_budget_saturated"] is True
    assert rows[0]["exhaustive_claim_allowed"] is False
    assert summary["trial_budget_saturated"] is True
    assert summary["exhaustive_claim_allowed"] is False


def test_campaign_comparison_reports_plateau_without_exhaustive_claim() -> None:
    def campaign(run_id: str, attempt_count: int) -> Manifest:
        manifest = Manifest.new(run_id=run_id, stage="discovery-audit")
        artifacts = [
            Artifact(
                artifact_id="source",
                artifact_type="species_preopt",
                data={
                    "species_id": "source",
                    "components": [
                        {"component_id": "a", "element_counts": {"H": 1}}
                    ],
                },
            ),
            Artifact(
                artifact_id="trial",
                artifact_type="reaction_trial",
                data={
                    "trial_id": "trial",
                    "source_species_id": "source",
                    "driver_order": ["nt2", "afir"],
                    "max_attempts": 2,
                },
            ),
        ]
        artifacts.extend(
            Artifact(
                artifact_id=f"attempt_{index}",
                artifact_type="reaction_discovery_attempt",
                data={"trial_id": "trial", "driver": f"driver_{index}"},
                status=ArtifactStatus(status="partial"),
            )
            for index in range(attempt_count)
        )
        manifest.extend(artifacts)
        return manifest

    rows, summary = build_discovery_saturation(
        [("small", campaign("small", 1)), ("large", campaign("large", 2))]
    )

    assert rows[-1]["budget_expanded"] is True
    assert rows[-1]["observed_yield_plateau"] is True
    assert rows[-1]["delta_discovery_attempt_count"] == 1
    assert rows[-1]["delta_candidate_count"] == 0
    assert summary["source_series_with_observed_plateau"] == 1
    assert summary["discovery_saturation_established"] is False
    assert summary["exhaustive_claim_allowed"] is False

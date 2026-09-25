from pathlib import Path

from hfauto.core.io import write_manifest
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import StageContext
from hfauto.stages.ts_search import TSSearchStage, _next_ts_workdir


def test_interrupted_ts_attempt_is_preserved_on_retry(tmp_path: Path) -> None:
    base = tmp_path / "attempt_00_nwchem_neb"
    assert _next_ts_workdir(base) == base
    base.mkdir()
    (base / "partial.raw").write_text("partial", encoding="utf-8")
    assert _next_ts_workdir(base) == tmp_path / "attempt_00_nwchem_neb_retry_02"
    retry = tmp_path / "attempt_00_nwchem_neb_retry_02"
    retry.mkdir()
    (retry / "partial.raw").write_text("second", encoding="utf-8")
    assert _next_ts_workdir(base) == tmp_path / "attempt_00_nwchem_neb_retry_03"


def _reaction(identifier: str) -> Artifact:
    return Artifact(
        artifact_id=identifier,
        artifact_type="reaction",
        data={
            "reaction_id": identifier,
            "reactant_species_id": f"{identifier}_reactant",
            "product_species_id": f"{identifier}_product",
        },
    )


def test_ts_stage_records_resource_deferral_as_nonchemical_result(
    tmp_path: Path,
) -> None:
    manifest = Manifest.new(run_id="budget", stage="reaction-plan")
    manifest.extend([_reaction("r1"), _reaction("r2")])

    result = TSSearchStage().run(
        manifest,
        {"max_reactions": 1, "max_total_attempts": 2},
        StageContext(out_dir=tmp_path, run_id="budget", global_config={}),
    )

    failures = result.latest_artifacts("ts_result")
    assert failures[0].status.category == "missing_input"
    assert failures[1].status.category == "execution_budget_exhausted"
    assert failures[1].data["budget_reason"] == "max_reactions_reached"
    assert "not evidence about reaction existence" in str(
        failures[1].status.recommended_fallback
    )


def test_ts_stage_resumes_frequency_validated_checkpoint(tmp_path: Path) -> None:
    source = Manifest.new(run_id="resume", stage="reaction-plan")
    source.add_artifact(_reaction("r1"))
    checkpoint = source.carry_forward("ts-search")
    checkpoint.add_artifact(
        Artifact(
            artifact_id="r1_with_ts",
            artifact_type="reaction_validated",
            data={"reaction_id": "r1"},
            qc={"ts_validated_by_frequency": True},
        )
    )
    write_manifest(checkpoint, tmp_path)

    result = TSSearchStage().run(
        source,
        {"resume_completed_reactions": True},
        StageContext(out_dir=tmp_path, run_id="resume", global_config={}),
    )

    assert result.metadata["ts_reactions_resumed"] == 1
    assert not result.latest_artifacts("ts_result")

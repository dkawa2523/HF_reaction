from pathlib import Path

from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto_ops.core.artifact_jobs import ARTIFACT_STAGE_HINTS
from hfauto_ops.core.resource import DEFAULT_STAGE_RESOURCES
from hfauto_ops.retry.policy import DEFAULT_RETRY_POLICY
from hfauto_ops.scheduler.resources import (
    STAGE_INPUT_TYPES,
    make_stage_job_plan,
)


def test_ops_maps_use_current_complex_and_irc_stages() -> None:
    assert ARTIFACT_STAGE_HINTS["conformer"] == "build-complexes"
    assert "build-complexes" in DEFAULT_STAGE_RESOURCES
    assert "build-hf" not in DEFAULT_STAGE_RESOURCES
    assert DEFAULT_RETRY_POLICY["endpoint_mismatch"]["stage"] == "irc"
    assert STAGE_INPUT_TYPES["build-complexes"] == ["conformer", "species"]


def test_scheduler_job_uses_selected_pipeline_config(tmp_path: Path) -> None:
    manifest = Manifest.new(run_id="run", stage="conformers")
    manifest.extend(
        [
            Artifact(
                artifact_id="conf_1",
                artifact_type="conformer",
                data={"species_id": "species_1"},
            )
        ]
    )

    jobs = make_stage_job_plan(
        manifest,
        tmp_path,
        ["build-complexes"],
        {"pipeline_config": "configs/pipelines/current.yaml"},
    )

    assert len(jobs) == 1
    assert "configs/pipelines/current.yaml" in jobs[0].command
    assert "phase9_hpc_template" not in jobs[0].command

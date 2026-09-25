from pathlib import Path

from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.sp import _single_point_inputs


def test_method_panel_uses_only_frequency_validated_transition_state(
    tmp_path: Path,
) -> None:
    seed_path = tmp_path / "seed.xyz"
    final_path = tmp_path / "final.xyz"
    seed_path.write_text("1\nseed\nH 0 0 0\n", encoding="utf-8")
    final_path.write_text("1\nfinal\nH 0 0 0\n", encoding="utf-8")
    manifest = Manifest.new(run_id="sp", stage="ts-search")
    manifest.extend(
        [
            Artifact(
                artifact_id="ts_seed",
                artifact_type="species",
                paths={"xyz": str(seed_path)},
                data={
                    "species_id": "ts",
                    "state": "transition_state",
                    "xyz_path": str(seed_path),
                },
            ),
            Artifact(
                artifact_id="ts_frequency",
                artifact_type="calculation",
                paths={"final_xyz": str(final_path)},
                data={"species_id": "ts", "state": "transition_state"},
                qc={
                    "ts_validated_by_frequency": True,
                    "real_qm_executed": True,
                },
            ),
        ]
    )

    inputs = _single_point_inputs(
        manifest,
        {"transition_state"},
        include_validated_transition_states=True,
    )

    assert [artifact.artifact_id for artifact in inputs] == ["ts_frequency"]

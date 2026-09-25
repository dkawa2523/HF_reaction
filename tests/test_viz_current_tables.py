from pathlib import Path

from hfauto.core.io import write_manifest
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto_viz.data.loaders import RunData


def test_visualization_uses_latest_reaction_ranking(tmp_path: Path) -> None:
    old_csv = tmp_path / "old.csv"
    new_csv = tmp_path / "new.csv"
    old_csv.write_text("reaction_id,value\nr1,1\n", encoding="utf-8")
    new_csv.write_text("reaction_id,value\nr1,2\n", encoding="utf-8")
    manifest = Manifest.new(run_id="run", stage="reaction-rank")
    manifest.extend(
        [
            Artifact(
                artifact_id="reaction_ranking",
                artifact_type="ranking",
                paths={"csv": str(old_csv)},
            ),
            Artifact(
                artifact_id="reaction_ranking",
                artifact_type="ranking",
                paths={"csv": str(new_csv)},
            ),
        ]
    )
    stage_dir = tmp_path / "00_reaction-rank"
    stage_dir.mkdir()
    write_manifest(manifest, stage_dir)

    table = RunData(tmp_path).read_table("reaction_results")

    assert table.iloc[0]["value"] == 2

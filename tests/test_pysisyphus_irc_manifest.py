from pathlib import Path

from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from scripts.run_pysisyphus_irc import assemble_result_manifest


def test_irc_result_manifest_preserves_upstream_calculations() -> None:
    source = Manifest.new(run_id="source", stage="ts-search")
    calculation = Artifact(
        artifact_id="calc_ts",
        artifact_type="calculation",
        data={"task": "saddle_freq"},
    )
    irc = Artifact(
        artifact_id="irc_rxn",
        artifact_type="irc",
        data={"reaction_id": "rxn"},
    )
    source.add_artifact(calculation)

    result = assemble_result_manifest(
        source,
        run_id="irc-run",
        source_path=Path("source/manifest.json"),
        reaction_id="rxn",
        validated_reaction_id="rxn_with_ts",
        method={"program": "nwchem"},
        artifacts=[irc],
    )

    assert result.run_id == "irc-run"
    assert result.stage == "bidirectional-irc"
    assert result.find("calc_ts") is calculation
    assert result.find("irc_rxn") is irc
    assert result.metadata["reaction_id"] == "rxn"

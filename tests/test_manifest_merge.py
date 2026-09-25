import pytest

from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest


def test_manifest_merge_preserves_branches_and_prefers_later_revision() -> None:
    minima = Manifest.new(run_id="minima", stage="minimum-registry")
    minima.extend(
        [
            Artifact(
                artifact_id="calc_minimum",
                artifact_type="calculation",
                data={"task": "opt_freq"},
            ),
            Artifact(
                artifact_id="endpoint",
                artifact_type="species_optimized",
                data={"revision": "minimum"},
            ),
        ]
    )
    path = Manifest.new(run_id="path", stage="irc")
    path.extend(
        [
            Artifact(
                artifact_id="endpoint",
                artifact_type="species_optimized",
                data={"revision": "path"},
            ),
            Artifact(
                artifact_id="irc",
                artifact_type="irc",
                data={"reaction_id": "rxn"},
            ),
        ]
    )

    merged = Manifest.merge(
        [minima, path], run_id="complete", stage="evidence-merge"
    )

    assert merged.parents == [minima.manifest_id, path.manifest_id]
    assert merged.find("calc_minimum") is not None
    assert merged.find("irc") is not None
    assert merged.find("endpoint").data["revision"] == "path"


def test_manifest_merge_rejects_empty_input() -> None:
    with pytest.raises(ValueError, match="at least one"):
        Manifest.merge([], run_id="empty", stage="evidence-merge")

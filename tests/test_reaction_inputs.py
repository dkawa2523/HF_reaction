import pytest

from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.workflow.reaction_inputs import (
    basin_assessment,
    reaction_and_endpoints,
    saddle_seed_candidate,
    validated_ts_and_endpoints,
)


def _manifest() -> Manifest:
    manifest = Manifest.new(run_id="inputs", stage="test")
    manifest.extend(
        [
            Artifact(
                artifact_id="rxn",
                artifact_type="reaction",
                data={
                    "reaction_id": "rxn",
                    "reactant_species_id": "r",
                    "product_species_id": "p",
                    "basin_assessment": {
                        "status": "distinct_basin",
                        "accepted": True,
                    },
                },
                qc={"distinct_registry_basins": True},
            ),
            Artifact(
                artifact_id="opt_r",
                artifact_type="species_optimized",
                data={"species_id": "r"},
            ),
            Artifact(
                artifact_id="opt_p",
                artifact_type="species_optimized",
                data={"species_id": "p"},
            ),
            Artifact(
                artifact_id="attempt",
                artifact_type="saddle_attempt",
                data={
                    "reaction_id": "rxn",
                    "diagnosis": "resolved_saddle_candidate",
                    "candidate": {"xyz_path": "seed.xyz"},
                },
                qc={
                    "saddle_seed_resolved": True,
                    "method_evidence_validated": True,
                },
            ),
            Artifact(
                artifact_id="ts",
                artifact_type="species",
                data={"species_id": "ts"},
                qc={"method_evidence_validated": True},
            ),
            Artifact(
                artifact_id="rxn_with_ts",
                artifact_type="reaction_validated",
                data={"reaction_id": "rxn", "ts_species_id": "ts"},
                qc={"ts_validated_by_frequency": True, "n_imag": 1},
            ),
        ]
    )
    return manifest


def test_reaction_input_selection_uses_validated_artifacts() -> None:
    manifest = _manifest()

    reaction, reactant, product = reaction_and_endpoints(manifest, "rxn")
    candidate, attempt = saddle_seed_candidate(manifest, "rxn")

    assert reaction.artifact_id == "rxn"
    assert reactant.artifact_id == "opt_r"
    assert product.artifact_id == "opt_p"
    assert basin_assessment(manifest, "rxn")["accepted"] is True
    assert candidate["xyz_path"] == "seed.xyz"
    assert attempt.artifact_id == "attempt"
    assert validated_ts_and_endpoints(manifest, "rxn")[1].artifact_id == "ts"


def test_saddle_seed_selection_rejects_unevidenced_candidate() -> None:
    manifest = _manifest()
    attempt = next(
        artifact
        for artifact in manifest.artifacts
        if artifact.artifact_type == "saddle_attempt"
    )
    attempt.qc["method_evidence_validated"] = False

    with pytest.raises(KeyError, match="Validated saddle seed"):
        saddle_seed_candidate(manifest, "rxn")

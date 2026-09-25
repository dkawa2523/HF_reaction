from pathlib import Path

import numpy as np

from hfauto.chemistry.xyz import XYZ, write_xyz
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import StageContext
from hfauto.stages.reaction_segments import ReactionSegmentsStage


def _minimum(
    artifact_id: str,
    xyz_path: Path,
    state: str,
    image_index: int | None = None,
) -> Artifact:
    data = {
        "species_id": artifact_id,
        "state": state,
        "xyz_path": str(xyz_path),
        "charge": 0,
        "multiplicity": 1,
    }
    if image_index is not None:
        data.update(
            {"source_image_index": image_index, "source_path_attempt_id": "path"}
        )
    return Artifact(
        artifact_id=artifact_id,
        artifact_type="species_optimized",
        paths={"xyz": str(xyz_path)},
        data=data,
        qc={
            "minimum_accepted": True,
            "is_minimum": True,
            "n_imag": 0,
            "scf_converged": True,
            "geometry_converged": True,
            "geometry_sane": True,
            "fallback_dummy": False,
        },
    )


def test_validated_intermediate_creates_independent_segment_hypotheses(
    tmp_path: Path,
) -> None:
    symbols = ["N", "H", "F"]
    geometries = [
        [[0.0, 0.0, 0.0], [1.70, 0.0, 0.0], [2.62, 0.0, 0.0]],
        [[0.0, 0.0, 0.0], [1.35, 0.7, 0.0], [2.45, 0.0, 0.0]],
        [[0.0, 0.0, 0.0], [1.05, 0.0, 0.0], [2.35, 0.0, 0.0]],
    ]
    paths = [
        write_xyz(XYZ(symbols, np.asarray(coords)), tmp_path / f"node_{i}.xyz")
        for i, coords in enumerate(geometries)
    ]
    reactant = _minimum("reactant", paths[0], "reactant")
    intermediate = _minimum(
        "intermediate", paths[1], "reaction_intermediate_candidate", 4
    )
    product = _minimum("product", paths[2], "product")
    reaction = Artifact(
        artifact_id="rxn",
        artifact_type="reaction",
        data={
            "reaction_id": "rxn",
            "reactant_species_id": "reactant",
            "product_species_id": "product",
            "reaction_coordinate": {"terms": []},
            "bond_changes": [],
            "basin_assessment": {
                "status": "distinct_basin",
                "accepted": True,
            },
        },
        qc={"distinct_registry_basins": True},
    )
    classification = Artifact(
        artifact_id="classification",
        artifact_type="reaction_classification",
        data={
            "reaction_id": "rxn",
            "classification": "multistep_with_intermediate",
            "scientific_conclusion_supported": True,
            "evidence": {
                "validated_intermediates": [
                    {
                        "artifact_id": "intermediate",
                        "source_image_index": 4,
                        "source_path_attempt_id": "path",
                    }
                ],
                "multistep_path_attempt_ids": ["path"],
            },
        },
    )
    manifest = Manifest.new(run_id="segments", stage="reaction-classify")
    manifest.extend([reaction, reactant, intermediate, product, classification])

    output = ReactionSegmentsStage().run(
        manifest,
        {},
        StageContext(
            out_dir=tmp_path / "segments",
            run_id="segments",
            global_config={},
        ),
    )
    plan = output.latest_artifacts("reaction_segmentation_plan")[-1]
    children = [
        item
        for item in output.latest_artifacts("reaction")
        if item.data.get("parent_reaction_id") == "rxn"
    ]

    assert plan.status.status == "success"
    assert plan.data["segment_count"] == 2
    assert len(children) == 2
    assert all(item.data["bond_changes"] == [] for item in children)
    assert all(item.data["reaction_coordinate"]["terms"] for item in children)
    assert all(item.qc["elementary_step_claimed"] is False for item in children)


def test_unresolved_multistep_candidate_is_not_split(tmp_path: Path) -> None:
    manifest = Manifest.new(run_id="segments", stage="reaction-classify")
    manifest.add_artifact(
        Artifact(
            artifact_id="classification",
            artifact_type="reaction_classification",
            data={
                "reaction_id": "rxn",
                "classification": "unresolved",
                "candidate_classification": "multistep_with_intermediate",
                "scientific_conclusion_supported": False,
            },
        )
    )

    output = ReactionSegmentsStage().run(
        manifest,
        {},
        StageContext(
            out_dir=tmp_path / "unresolved",
            run_id="segments",
            global_config={},
        ),
    )

    assert output.latest_artifacts("reaction_segmentation_plan") == []

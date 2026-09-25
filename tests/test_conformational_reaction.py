from pathlib import Path

from hfauto.chemistry.reactions import (
    validate_reaction_coordinate_between_geometries,
    validate_reaction_endpoint_pair,
)
from hfauto.core.schemas.artifact import Artifact
from hfauto.stages.base import StageContext
from hfauto.stages.build_complexes import BuildComplexesStage
from hfauto.workflow.runner import run_pipeline

ROOT = Path(__file__).resolve().parents[1]
REACTANT = ROOT / "configs" / "systems" / "xyz" / "nh3_inversion" / "reactant.xyz"
PRODUCT = ROOT / "configs" / "systems" / "xyz" / "nh3_inversion" / "product.xyz"
COORDINATE = {
    "bond_changes": [],
    "reaction_coordinate": {
        "min_change": 0.20,
        "terms": [
            {
                "kind": "dihedral",
                "atoms": [0, 1, 2, 3],
                "coefficient": 1.0,
                "label": "ammonia_umbrella_improper",
            }
        ],
    },
}


def test_conformational_coordinate_distinguishes_inversion_endpoints() -> None:
    qc = validate_reaction_coordinate_between_geometries(
        COORDINATE, REACTANT, PRODUCT
    )

    assert qc["accepted"] is True
    assert qc["progress"] >= 0.20


def test_endpoint_pair_accepts_coordinate_without_bond_changes() -> None:
    reactant = Artifact(
        artifact_id="reactant",
        artifact_type="species",
        paths={"xyz": str(REACTANT)},
        data={"xyz_path": str(REACTANT), "charge": 0, "multiplicity": 1},
    )
    product = Artifact(
        artifact_id="product",
        artifact_type="species",
        paths={"xyz": str(PRODUCT)},
        data={"xyz_path": str(PRODUCT), "charge": 0, "multiplicity": 1},
    )

    qc = validate_reaction_endpoint_pair(COORDINATE, reactant, product)

    assert qc["accepted"] is True
    assert qc["validation_scope"].endswith("declared_reaction_coordinate")


def test_build_complexes_can_start_from_reviewed_xyz(tmp_path: Path) -> None:
    out = BuildComplexesStage().run(
        None,
        {
            "systems": [
                {
                    "system_id": "ammonia",
                    "state": "reactant",
                    "xyz_path": str(REACTANT),
                    "charge": 0,
                    "multiplicity": 1,
                    "atom_order_key": "N_H1_H2_H3",
                }
            ]
        },
        StageContext(out_dir=tmp_path, run_id="m3_test", global_config={}),
    )

    species = out.latest_artifacts("species")
    assert len(species) == 1
    assert species[0].data["source_xyz_sha256"]
    assert Path(species[0].paths["xyz"]).is_file()


def test_pipeline_can_start_from_source_capable_stage(tmp_path: Path) -> None:
    run_dir = run_pipeline(
        {
            "run_root": str(tmp_path / "runs"),
            "stages": [
                {
                    "name": "build-complexes",
                    "systems": [
                        {
                            "system_id": "ammonia",
                            "state": "reactant",
                            "xyz_path": str(REACTANT),
                        }
                    ],
                }
            ],
        },
        "source_pipeline",
    )

    assert (run_dir / "00_build-complexes" / "manifest.json").is_file()

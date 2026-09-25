from pathlib import Path

from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import StageContext
from hfauto.stages.recover_path import RecoverPathStage


def _xyz(distance: float, energy: float) -> str:
    return f"2\nenergy={energy}\nH 0 0 0\nF 0 0 {distance}\n"


def test_interrupted_path_is_recovered_only_as_saddle_seed(tmp_path: Path) -> None:
    reactant_path = tmp_path / "reactant.xyz"
    product_path = tmp_path / "product.xyz"
    path_xyz = tmp_path / "path.xyz"
    reactant_path.write_text(_xyz(1.0, -10.0), encoding="utf-8")
    product_path.write_text(_xyz(1.4, -10.1), encoding="utf-8")
    path_xyz.write_text(
        _xyz(1.0, -10.0)
        + _xyz(1.2, -9.8)
        + _xyz(1.4, -10.1),
        encoding="utf-8",
    )
    manifest = Manifest.new(run_id="recover", stage="reaction-plan")
    manifest.extend(
        [
            Artifact(
                artifact_id="reactant",
                artifact_type="species_optimized",
                paths={"xyz": str(reactant_path)},
                data={"species_id": "reactant", "xyz_path": str(reactant_path)},
            ),
            Artifact(
                artifact_id="product",
                artifact_type="species_optimized",
                paths={"xyz": str(product_path)},
                data={"species_id": "product", "xyz_path": str(product_path)},
            ),
            Artifact(
                artifact_id="reaction",
                artifact_type="reaction",
                data={
                    "reaction_id": "reaction",
                    "reactant_species_id": "reactant",
                    "product_species_id": "product",
                    "basin_assessment": {"status": "distinct_basin", "accepted": True},
                },
                qc={"distinct_registry_basins": True},
            ),
        ]
    )

    result = RecoverPathStage().run(
        manifest,
        {"paths": [{"reaction_id": "reaction", "path_xyz": str(path_xyz)}]},
        StageContext(out_dir=tmp_path / "recovered", run_id="recover", global_config={}),
    )

    recovered_path = result.latest_artifacts("reaction_path")[-1]
    saddle_attempt = result.latest_artifacts("saddle_attempt")[-1]
    case = result.latest_artifacts("reaction_case")[-1]
    assert recovered_path.qc["path_energy_publishable"] is False
    assert recovered_path.status.status == "partial"
    assert saddle_attempt.data["diagnosis"] == "resolved_saddle_candidate"
    assert saddle_attempt.data["candidate"]["path_image_index"] == 1
    assert case.data["next_strategy"] == "saddle_refinement"
    assert case.data["ts_search_allowed"] is True


def test_recovered_string_uses_string_optimization_history(tmp_path: Path) -> None:
    reactant_path = tmp_path / "reactant.xyz"
    product_path = tmp_path / "product.xyz"
    path_xyz = tmp_path / "path.xyz"
    history = tmp_path / "nwchem_string.out"
    reactant_path.write_text(_xyz(1.0, -10.0), encoding="utf-8")
    product_path.write_text(_xyz(1.4, -10.1), encoding="utf-8")
    path_xyz.write_text(
        _xyz(1.0, -10.0)
        + _xyz(1.2, -9.8)
        + _xyz(1.4, -10.1),
        encoding="utf-8",
    )
    history.write_text(
        "@zts  4 0.003600 0.090000 -10 -9 -10 -9 100\n",
        encoding="utf-8",
    )
    manifest = Manifest.new(run_id="recover_string", stage="reaction-plan")
    manifest.extend(
        [
            Artifact(
                artifact_id="reactant",
                artifact_type="species_optimized",
                paths={"xyz": str(reactant_path)},
                data={"species_id": "reactant", "xyz_path": str(reactant_path)},
            ),
            Artifact(
                artifact_id="product",
                artifact_type="species_optimized",
                paths={"xyz": str(product_path)},
                data={"species_id": "product", "xyz_path": str(product_path)},
            ),
            Artifact(
                artifact_id="reaction",
                artifact_type="reaction",
                data={
                    "reaction_id": "reaction",
                    "reactant_species_id": "reactant",
                    "product_species_id": "product",
                    "basin_assessment": {
                        "status": "distinct_basin",
                        "accepted": True,
                    },
                },
                qc={"distinct_registry_basins": True},
            ),
        ]
    )

    result = RecoverPathStage().run(
        manifest,
        {
            "paths": [
                {
                    "reaction_id": "reaction",
                    "engine": "nwchem_string",
                    "path_xyz": str(path_xyz),
                    "optimization_history": str(history),
                }
            ]
        },
        StageContext(
            out_dir=tmp_path / "recovered_string",
            run_id="recover_string",
            global_config={},
        ),
    )

    attempt = result.latest_artifacts("path_attempt")[-1]
    assert attempt.data["engine"] == "nwchem_string"
    assert attempt.data["optimization"]["iterations"] == 4
    assert attempt.data["optimization"]["rms_displacement"] == 0.0036

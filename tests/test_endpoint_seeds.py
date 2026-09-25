from pathlib import Path

import numpy as np

from hfauto.chemistry.complex_seed_ensemble import generate_rigid_fragment_seeds
from hfauto.chemistry.xyz import XYZ, read_xyz, write_xyz
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import StageContext
from hfauto.stages.endpoint_seed_screen import EndpointSeedScreenStage
from hfauto.stages.endpoint_seeds import EndpointSeedsStage


def _complex() -> XYZ:
    return XYZ(
        symbols=["N", "H", "H", "H", "H", "F"],
        coords=np.asarray(
            [
                [0.0, 0.0, 0.0],
                [0.9, 0.0, 0.0],
                [-0.3, 0.85, 0.0],
                [-0.3, -0.4, 0.75],
                [0.0, 0.0, 2.0],
                [0.0, 0.0, 2.92],
            ]
        ),
    )


def _reaction() -> Artifact:
    return Artifact(
        artifact_id="rxn",
        artifact_type="reaction",
        data={
            "reaction_id": "rxn",
            "reactant_species_id": "reactant",
            "product_species_id": "product",
            "mechanism_family": "intermolecular_proton_transfer",
            "bond_changes": [],
            "reaction_coordinate": {
                "terms": [
                    {
                        "kind": "distance",
                        "atoms": [0, 4],
                        "coefficient": -1.0,
                    },
                    {
                        "kind": "distance",
                        "atoms": [4, 5],
                        "coefficient": 1.0,
                    },
                ]
            },
        },
    )


def test_rigid_fragment_seeds_preserve_internal_geometry() -> None:
    source = _complex()
    seeds = generate_rigid_fragment_seeds(
        source,
        _reaction().data,
        maximum_seeds=4,
        minimum_interfragment_distance_A=0.60,
    )

    assert len(seeds) == 4
    source_hf = np.linalg.norm(source.coords[4] - source.coords[5])
    assert all(
        np.isclose(np.linalg.norm(seed.coords[4] - seed.coords[5]), source_hf)
        for seed, _evidence in seeds
    )
    assert all(evidence["minimum_claimed"] is False for _, evidence in seeds)
    assert len(
        {
            tuple(np.round(seed.coords.reshape(-1), 6))
            for seed, _evidence in seeds
        }
    ) == 4


def test_endpoint_seed_stage_is_bounded_and_keeps_root_lineage(
    tmp_path: Path,
) -> None:
    xyz_path = write_xyz(_complex(), tmp_path / "reactant.xyz")
    reaction = _reaction()
    source = Artifact(
        artifact_id="reactant",
        artifact_type="species",
        paths={"xyz": str(xyz_path)},
        data={
            "species_id": "reactant",
            "state": "reactant",
            "charge": 0,
            "multiplicity": 1,
            "xyz_path": str(xyz_path),
        },
    )
    calculation = Artifact(
        artifact_id="calc_reactant",
        artifact_type="calculation",
        paths={"final_xyz": str(xyz_path)},
        data={"species_id": "reactant", "task": "opt_freq", "n_imag": 1},
        method={"stage": "dft-minima"},
    )
    selection = Artifact(
        artifact_id="endpoint_pair_selection_rxn",
        artifact_type="endpoint_pair_selection",
        data={
            "reaction_id": "rxn",
            "reactant_root_species_id": "reactant",
            "product_root_species_id": "product",
            "reasons": ["all_endpoint_candidates_collapse_to_same_basin"],
        },
        qc={"endpoint_pair_ready": False},
    )
    manifest = Manifest.new(run_id="test", stage="endpoint-seeds")
    manifest.extend([reaction, source, calculation, selection])

    output = EndpointSeedsStage().run(
        manifest,
        {
            "roles": ["reactant"],
            "maximum_seeds_per_role": 2,
            "minimum_interfragment_distance_A": 0.60,
        },
        StageContext(
            out_dir=tmp_path / "seeds",
            run_id="test",
            global_config={"mode": "production"},
        ),
    )
    seeds = [
        artifact
        for artifact in output.latest_artifacts("species")
        if artifact.data.get("endpoint_seed_generation") == 1
    ]
    plan = output.latest_artifacts("endpoint_seed_ensemble_plan")[-1]

    assert len(seeds) == 2
    assert plan.status.status == "success"
    assert plan.qc["bounded_generation"] is True
    assert all(seed.data["source_species_id"] == "reactant" for seed in seeds)
    assert all(read_xyz(seed.paths["xyz"]).symbols == _complex().symbols for seed in seeds)


def test_endpoint_seed_stage_does_not_recurse(tmp_path: Path) -> None:
    first = test_manifest = Manifest.new(run_id="test", stage="endpoint-seeds")
    plan = Artifact(
        artifact_id="endpoint_seed_plan_rxn_reactant_existing",
        artifact_type="endpoint_seed_ensemble_plan",
        data={"reaction_id": "rxn"},
    )
    test_manifest.add_artifact(plan)
    # The public generation contract rejects a second generation before any
    # chemistry or filesystem work is attempted.
    try:
        EndpointSeedsStage().run(
            first,
            {"generation": 2},
            StageContext(
                out_dir=tmp_path / "second",
                run_id="test",
                global_config={},
            ),
        )
    except ValueError as exc:
        assert "exactly one bounded generation" in str(exc)
    else:
        raise AssertionError("a recursive endpoint seed generation was accepted")


def test_endpoint_seed_screen_uses_energy_window_and_geometry_diversity(
    tmp_path: Path,
) -> None:
    manifest = Manifest.new(run_id="screen", stage="preopt")
    geometries = [_complex(), _complex(), _complex()]
    geometries[1].coords[4:] += np.asarray([0.001, 0.0, 0.0])
    geometries[2].coords[4:] += np.asarray([0.4, 0.2, 0.0])
    energies = [-10.0000, -9.9999, -9.9998]
    for index, (geometry, energy) in enumerate(zip(geometries, energies)):
        species_id = f"seed_{index}"
        xyz_path = write_xyz(geometry, tmp_path / f"{species_id}.xyz")
        manifest.add_artifact(
            Artifact(
                artifact_id=species_id,
                artifact_type="species",
                paths={"xyz": str(xyz_path)},
                data={
                    "species_id": species_id,
                    "state": "endpoint_minimum_candidate",
                    "xyz_path": str(xyz_path),
                    "endpoint_seed_generation": 1,
                    "endpoint_role": "reactant",
                    "source_species_id": "root",
                    "reaction_ids": ["rxn"],
                },
                method={"stage": "preopt"},
            )
        )
        manifest.add_artifact(
            Artifact(
                artifact_id=f"calc_{species_id}",
                artifact_type="calculation",
                data={
                    "species_id": species_id,
                    "electronic_energy_hartree": energy,
                },
                method={"stage": "preopt", "engine": "xtb"},
                qc={
                    "geometry_converged": True,
                    "geometry_sane": True,
                    "fallback_dummy": False,
                    "preopt_promotion_accepted": True,
                },
            )
        )

    output = EndpointSeedScreenStage().run(
        manifest,
        {
            "maximum_selected_per_role": 2,
            "energy_window_kcal_mol": 1.0,
            "minimum_distance_rmsd_A": 0.03,
        },
        StageContext(
            out_dir=tmp_path / "screen",
            run_id="screen",
            global_config={},
        ),
    )
    selection = output.latest_artifacts("endpoint_seed_selection")[-1]

    assert selection.status.status == "success"
    assert selection.data["selected_species_ids"] == ["seed_0", "seed_2"]
    assert selection.data["preopt_energy_is_final_chemistry"] is False

from __future__ import annotations

from pathlib import Path

import numpy as np

from hfauto.backends.conformer.crest import CRESTConformerBackend
from hfauto.chemistry.reaction_trials import (
    automatic_reaction_trials,
    deduplicate_trials,
    explicit_reaction_trials,
)
from hfauto.chemistry.xyz import XYZ, read_xyz, write_xyz
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import StageContext
from hfauto.stages.build_complexes import BuildComplexesStage
from hfauto.stages.generate_reactions import GenerateReactionsStage


def test_explicit_conformational_coordinate_is_a_product_free_trial() -> None:
    trials = explicit_reaction_trials(
        "ammonia",
        4,
        [
            {
                "system_id": "ammonia",
                "mechanism_hint": "umbrella_inversion",
                "reaction_coordinate": {
                    "min_change": 0.2,
                    "terms": [
                        {
                            "kind": "dihedral",
                            "atoms": [0, 1, 2, 3],
                            "coefficient": 1.0,
                        }
                    ],
                },
            }
        ],
        component_ids=["ammonia"],
        charge=0,
        multiplicity=1,
        default_driver_order=["dft_minimum_ensemble"],
        default_max_attempts=1,
    )

    assert len(trials) == 1
    assert trials[0].associations == []
    assert trials[0].dissociations == []
    assert trials[0].driver_order == ["dft_minimum_ensemble"]
    assert trials[0].reaction_coordinate.terms[0].kind == "dihedral"


def test_crest_species_identity_cannot_overwrite_a_molecular_conformer() -> None:
    complex_species = Artifact(
        artifact_id="spc_tma_hf2_seed00",
        artifact_type="species",
        data={
            "species_id": "spc_tma_hf2_seed00",
            "mol_id": "mol00001_state000",
        },
    )
    molecule_state = Artifact(
        artifact_id="mol00001_state000",
        artifact_type="molecule_state",
        data={"mol_id": "mol00001_state000"},
    )

    assert (
        CRESTConformerBackend._calculation_identity(complex_species)
        == "spc_tma_hf2_seed00"
    )
    assert (
        CRESTConformerBackend._calculation_identity(molecule_state)
        == "mol00001_state000"
    )


def test_explicit_trial_has_no_product_geometry() -> None:
    trials = explicit_reaction_trials(
        "spc_hcn",
        3,
        [
            {
                "system_id": "hcn",
                "associations": [[1, 2]],
                "dissociations": [[0, 2]],
                "driver_order": ["nt2"],
            }
        ],
        component_ids=["hcn"],
        charge=0,
        multiplicity=1,
        default_driver_order=["nt2", "afir"],
        default_max_attempts=2,
    )

    assert len(trials) == 1
    record = trials[0].model_dump(mode="json")
    assert record["source_species_id"] == "spc_hcn"
    assert record["associations"][0]["atoms"] == [1, 2]
    assert record["dissociations"][0]["atoms"] == [0, 2]
    assert "product" not in record
    assert "product_xyz_path" not in record


def test_automatic_trials_include_transfer_relay_and_are_deduplicated() -> None:
    xyz = XYZ(
        symbols=["N", "H", "H", "H", "H", "F", "H", "F"],
        coords=np.asarray(
            [
                [0.00, 0.00, 0.00],
                [0.95, 0.00, 0.00],
                [-0.32, 0.90, 0.00],
                [-0.32, -0.45, 0.78],
                [1.65, 0.10, 0.00],
                [2.57, 0.10, 0.00],
                [3.05, 0.10, 0.00],
                [3.97, 0.10, 0.00],
            ]
        ),
    )
    components = [
        {"component_id": "amine", "atom_indices": [0, 1, 2, 3]},
        {"component_id": "hf_a", "atom_indices": [4, 5]},
        {"component_id": "hf_b", "atom_indices": [6, 7]},
    ]

    trials = automatic_reaction_trials(
        "spc_amine_hf2",
        xyz,
        components,
        charge=0,
        multiplicity=1,
        driver_order=["nt2", "afir"],
        max_attempts=2,
        max_relay_steps=2,
        relay_contact_cutoff_A=3.2,
    )
    mechanisms = {trial.mechanism_hint for trial in trials}
    assert "proton_transfer" in mechanisms
    assert "proton_relay" in mechanisms
    assert len(deduplicate_trials([*trials, *trials])) == len(
        deduplicate_trials(trials)
    )


def _source(
    tmp_path: Path,
    artifact_id: str,
    symbols: list[str],
    coords: list[list[float]],
) -> Artifact:
    path = write_xyz(XYZ(symbols, np.asarray(coords)), tmp_path / f"{artifact_id}.xyz")
    return Artifact(
        artifact_id=artifact_id,
        artifact_type="conformer",
        paths={"xyz": str(path)},
        data={
            "conformer_id": artifact_id,
            "mol_id": artifact_id,
            "xyz_path": str(path),
            "charge": 0,
            "multiplicity": 1,
        },
    )


def test_composition_count_builds_three_component_complex(tmp_path: Path) -> None:
    manifest = Manifest.new(run_id="test", stage="conformers")
    manifest.extend(
        [
            _source(
                tmp_path,
                "amine",
                ["N", "H", "H", "H"],
                [
                    [0.0, 0.0, 0.0],
                    [0.95, 0.0, 0.0],
                    [-0.3, 0.9, 0.0],
                    [-0.3, -0.45, 0.78],
                ],
            ),
            _source(
                tmp_path,
                "hf",
                ["H", "F"],
                [[0.0, 0.0, 0.0], [0.92, 0.0, 0.0]],
            ),
        ]
    )
    result = BuildComplexesStage().run(
        manifest,
        {
            "compositions": [
                {
                    "composition_id": "amine_hf2",
                    "components": [
                        {"selector": {"artifact_id": "amine"}, "count": 1},
                        {"selector": {"artifact_id": "hf"}, "count": 2},
                    ],
                    "max_seeds": 1,
                    "crest_nci": False,
                }
            ]
        },
        StageContext(tmp_path / "out", "test", {}),
    )

    complexes = [
        item
        for item in result.latest_artifacts("species")
        if item.data.get("composition_id") == "amine_hf2"
    ]
    assert len(complexes) == 1
    assert complexes[0].qc["component_count"] == 3
    assert len(complexes[0].data["components"]) == 3
    assert complexes[0].data["element_counts"] == {"F": 2, "H": 5, "N": 1}


def test_multicomponent_species_does_not_inherit_one_fragment_molecule_id(
    tmp_path: Path,
) -> None:
    amine = _source(
        tmp_path,
        "amine",
        ["N", "H", "H", "H"],
        [[0, 0, 0], [1, 0, 0], [-0.3, 0.9, 0], [-0.3, -0.45, 0.78]],
    )
    hf = _source(tmp_path, "hf", ["H", "F"], [[0, 0, 0], [0.92, 0, 0]])
    hf.data.pop("mol_id")
    manifest = Manifest.new(run_id="identity", stage="conformers")
    manifest.extend([amine, hf])

    result = BuildComplexesStage().run(
        manifest,
        {
            "compositions": [
                {
                    "composition_id": "amine_hf",
                    "components": [
                        {"selector": {"artifact_id": "amine"}},
                        {"selector": {"artifact_id": "hf"}},
                    ],
                    "max_seeds": 1,
                    "crest_nci": False,
                }
            ]
        },
        StageContext(tmp_path / "out", "identity", {}),
    )
    complex_species = next(
        item
        for item in result.latest_artifacts("species")
        if item.data.get("composition_id") == "amine_hf"
    )

    assert "mol_id" not in complex_species.data


def test_hf_stoichiometry_series_keeps_the_original_base_geometry(
    tmp_path: Path, monkeypatch
) -> None:
    def fake_generate(backend, species, _config, workdir):
        identity = backend._calculation_identity(species)
        xyz_path = write_xyz(
            read_xyz(species.paths["xyz"]), Path(workdir) / f"{identity}.xyz"
        )
        return [
            Artifact(
                artifact_id=f"{identity}_conf0000",
                artifact_type="conformer",
                parents=[species.artifact_id],
                paths={"xyz": str(xyz_path)},
                data={"relative_energy_kcal_mol": 0.0},
            )
        ]

    monkeypatch.setattr(
        "hfauto.stages.build_complexes.CRESTConformerBackend.generate",
        fake_generate,
    )
    amine = _source(
        tmp_path,
        "amine",
        ["N", "H", "H", "H"],
        [[0, 0, 0], [1, 0, 0], [-0.3, 0.9, 0], [-0.3, -0.45, 0.78]],
    )
    hf = _source(tmp_path, "hf", ["H", "F"], [[0, 0, 0], [0.92, 0, 0]])
    hf.data.pop("mol_id")
    manifest = Manifest.new(run_id="stoichiometry", stage="conformers")
    manifest.extend([amine, hf])
    compositions = [
        {
            "composition_id": f"NH3_HF{count}",
            "components": [
                {"selector": {"artifact_id": "amine"}},
                {"selector": {"artifact_id": "hf"}, "count": count},
            ],
            "max_seeds": 1,
            "max_nci_complexes": 1,
        }
        for count in (1, 2, 3)
    ]

    result = BuildComplexesStage().run(
        manifest,
        {"compositions": compositions, "max_components": 4},
        StageContext(tmp_path / "out", "stoichiometry", {}),
    )

    for count, atom_count in ((1, 6), (2, 8), (3, 10)):
        seed = result.find(f"spc_NH3_HF{count}_seed00")
        assert seed is not None
        assert len(read_xyz(seed.paths["xyz"]).symbols) == atom_count
    assert result.find("amine") is amine


def test_build_complexes_rejects_cross_lineage_artifact_id_collision(
    tmp_path: Path, monkeypatch
) -> None:
    def colliding_generate(_backend, species, _config, workdir):
        xyz_path = write_xyz(
            read_xyz(species.paths["xyz"]), Path(workdir) / "collision.xyz"
        )
        return [
            Artifact(
                artifact_id="amine",
                artifact_type="conformer",
                parents=[species.artifact_id],
                paths={"xyz": str(xyz_path)},
                data={"relative_energy_kcal_mol": 0.0},
            )
        ]

    monkeypatch.setattr(
        "hfauto.stages.build_complexes.CRESTConformerBackend.generate",
        colliding_generate,
    )
    amine = _source(
        tmp_path,
        "amine",
        ["N", "H", "H", "H"],
        [[0, 0, 0], [1, 0, 0], [-0.3, 0.9, 0], [-0.3, -0.45, 0.78]],
    )
    hf = _source(tmp_path, "hf", ["H", "F"], [[0, 0, 0], [0.92, 0, 0]])
    manifest = Manifest.new(run_id="collision", stage="conformers")
    manifest.extend([amine, hf])

    result = BuildComplexesStage().run(
        manifest,
        {
            "compositions": [
                {
                    "composition_id": "amine_hf",
                    "components": [
                        {"selector": {"artifact_id": "amine"}},
                        {"selector": {"artifact_id": "hf"}},
                    ],
                    "max_seeds": 1,
                    "max_nci_complexes": 1,
                }
            ]
        },
        StageContext(tmp_path / "out", "collision", {}),
    )

    assert result.find("amine") is amine
    failure = result.find("build_complex_identity_collision_amine_hf_seed00")
    assert failure is not None
    assert failure.status.category == "artifact_identity_collision"
    assert result.metadata["build_complex_failed_definition_ids"] == [
        "amine_hf_seed00"
    ]


def test_configured_limit_allows_base_plus_three_hf(tmp_path: Path) -> None:
    manifest = Manifest.new(run_id="test", stage="conformers")
    manifest.extend(
        [
            _source(tmp_path, "base", ["N", "H"], [[0, 0, 0], [1, 0, 0]]),
            _source(tmp_path, "hf", ["H", "F"], [[0, 0, 0], [0.92, 0, 0]]),
        ]
    )

    result = BuildComplexesStage().run(
        manifest,
        {
            "max_components": 4,
            "compositions": [
                {
                    "composition_id": "base_hf3",
                    "components": [
                        {"selector": {"artifact_id": "base"}},
                        {"selector": {"artifact_id": "hf"}, "count": 3},
                    ],
                    "max_seeds": 1,
                    "crest_nci": False,
                }
            ],
        },
        StageContext(tmp_path / "out", "test", {}),
    )

    complex_species = next(
        item
        for item in result.latest_artifacts("species")
        if item.data.get("composition_id") == "base_hf3"
    )
    assert complex_species.qc["component_count"] == 4
    assert complex_species.data["element_counts"] == {"F": 3, "H": 4, "N": 1}


def test_build_complexes_checkpoints_and_resumes_each_definition(
    tmp_path: Path, monkeypatch
) -> None:
    calls: list[Path] = []

    def fake_generate(_backend, species, _config, workdir):
        root = Path(workdir)
        calls.append(root)
        xyz_path = write_xyz(
            read_xyz(species.paths["xyz"]),
            root / "fake_nci.xyz",
        )
        return [
            Artifact(
                artifact_id=f"{species.artifact_id}_fake_conf",
                artifact_type="conformer",
                parents=[species.artifact_id],
                paths={"xyz": str(xyz_path)},
                data={"relative_energy_kcal_mol": 0.0},
                method={"engine": "crest"},
                qc={"nci_search": True},
            )
        ]

    monkeypatch.setattr(
        "hfauto.stages.build_complexes.CRESTConformerBackend.generate",
        fake_generate,
    )
    base_path = write_xyz(
        XYZ(
            ["N", "H", "H", "H"],
            np.asarray(
                [[0, 0, 0], [1, 0, 0], [-0.3, 0.9, 0], [-0.3, -0.45, 0.78]],
                dtype=float,
            ),
        ),
        tmp_path / "base.xyz",
    )
    hf_path = write_xyz(
        XYZ(["H", "F"], np.asarray([[0, 0, 0], [0.92, 0, 0]], dtype=float)),
        tmp_path / "hf.xyz",
    )
    settings = {
        "systems": [
            {"system_id": "base", "xyz_path": str(base_path)},
            {"system_id": "hf", "xyz_path": str(hf_path)},
        ],
        "compositions": [
            {
                "composition_id": "base_hf",
                "components": [
                    {"selector": {"artifact_id": "spc_base"}},
                    {"selector": {"artifact_id": "spc_hf"}},
                ],
                "max_seeds": 1,
                "max_nci_complexes": 1,
            }
        ],
    }
    context = StageContext(tmp_path / "out", "resume", {"mode": "production"})

    first = BuildComplexesStage().run(None, settings, context)
    second = BuildComplexesStage().run(None, settings, context)

    assert len(calls) == 1
    assert calls[0].name == "attempt_00"
    assert (tmp_path / "out" / "manifest.json").is_file()
    assert second.metadata["build_complex_resumed_definition_count"] == 3
    assert "base_hf_seed00" in second.metadata[
        "build_complex_completed_definition_ids"
    ]
    assert any(
        item.qc.get("crest_nci_searched")
        for item in first.latest_artifacts("species")
    )

    changed_settings = {**settings, "crest_nci_settings": {"ewin": 6.0}}
    BuildComplexesStage().run(None, changed_settings, context)
    assert [path.name for path in calls] == ["attempt_00", "attempt_01"]


def test_failed_nci_retry_uses_a_new_raw_attempt_directory(
    tmp_path: Path, monkeypatch
) -> None:
    calls: list[Path] = []

    def fail_then_succeed(_backend, species, _config, workdir):
        root = Path(workdir)
        calls.append(root)
        if len(calls) == 1:
            return [
                Artifact.failure(
                    f"failed_{species.artifact_id}",
                    "conformer",
                    "synthetic interruption",
                    category="crest_failed",
                )
            ]
        xyz_path = write_xyz(read_xyz(species.paths["xyz"]), root / "retry.xyz")
        return [
            Artifact(
                artifact_id=f"{species.artifact_id}_retry_conf",
                artifact_type="conformer",
                parents=[species.artifact_id],
                paths={"xyz": str(xyz_path)},
                data={"relative_energy_kcal_mol": 0.0},
            )
        ]

    monkeypatch.setattr(
        "hfauto.stages.build_complexes.CRESTConformerBackend.generate",
        fail_then_succeed,
    )
    source = Manifest.new(run_id="retry", stage="conformers")
    source.extend(
        [
            _source(
                tmp_path,
                "amine",
                ["N", "H", "H", "H"],
                [[0, 0, 0], [1, 0, 0], [-0.3, 0.9, 0], [-0.3, -0.45, 0.78]],
            ),
            _source(tmp_path, "hf", ["H", "F"], [[0, 0, 0], [0.92, 0, 0]]),
        ]
    )
    settings = {
        "compositions": [
            {
                "composition_id": "amine_hf",
                "components": [
                    {"selector": {"artifact_id": "amine"}},
                    {"selector": {"artifact_id": "hf"}},
                ],
                "max_seeds": 1,
                "max_nci_complexes": 1,
            }
        ]
    }
    context = StageContext(tmp_path / "out", "retry", {"mode": "production"})

    BuildComplexesStage().run(source, settings, context)
    result = BuildComplexesStage().run(source, settings, context)
    BuildComplexesStage().run(source, settings, context)

    assert [path.name for path in calls] == ["attempt_00", "attempt_01"]
    assert result.metadata["build_complex_failed_definition_ids"] == []
    assert result.metadata["build_complex_nci_counts"] == {"amine_hf": 1}


def test_realistic_ammonia_plus_three_hf_uses_collision_alternatives(
    tmp_path: Path,
) -> None:
    root = Path(__file__).resolve().parents[1]
    result = BuildComplexesStage().run(
        None,
        {
            "systems": [
                {
                    "system_id": "ammonia",
                    "xyz_path": str(
                        root / "examples/m3_ammonia_inversion/reactant.xyz"
                    ),
                },
                {
                    "system_id": "hf",
                    "xyz_path": str(root / "examples/reaction_discovery/hf.xyz"),
                },
            ],
            "max_components": 4,
            "compositions": [
                {
                    "composition_id": "ammonia_hf3",
                    "components": [
                        {"selector": {"artifact_id": "spc_ammonia"}},
                        {"selector": {"artifact_id": "spc_hf"}, "count": 3},
                    ],
                    "max_seeds": 1,
                    "crest_nci": False,
                }
            ],
        },
        StageContext(tmp_path / "out", "test", {}),
    )

    complex_species = next(
        item
        for item in result.latest_artifacts("species")
        if item.data.get("composition_id") == "ammonia_hf3"
    )
    assert complex_species.qc["component_count"] == 4
    assert all(
        placement["minimum_interfragment_distance_A"] >= 0.55
        for placement in complex_species.qc["placements"]
    )


def test_composition_can_select_a_system_created_in_the_same_stage(
    tmp_path: Path,
) -> None:
    base_path = write_xyz(
        XYZ(["N", "H"], np.asarray([[0, 0, 0], [1, 0, 0]], dtype=float)),
        tmp_path / "base.xyz",
    )
    hf_path = write_xyz(
        XYZ(["H", "F"], np.asarray([[0, 0, 0], [0.92, 0, 0]], dtype=float)),
        tmp_path / "hf.xyz",
    )

    result = BuildComplexesStage().run(
        None,
        {
            "systems": [
                {"system_id": "base", "xyz_path": str(base_path)},
                {"system_id": "hf", "xyz_path": str(hf_path)},
            ],
            "compositions": [
                {
                    "composition_id": "base_hf",
                    "components": [
                        {"selector": {"artifact_id": "spc_base"}},
                        {"selector": {"artifact_id": "spc_hf"}},
                    ],
                    "max_seeds": 1,
                    "crest_nci": False,
                }
            ],
        },
        StageContext(tmp_path / "out", "test", {}),
    )

    assert result.find("spc_base") is not None
    assert result.find("spc_hf") is not None
    assert any(
        artifact.data.get("composition_id") == "base_hf"
        for artifact in result.latest_artifacts("species")
    )


def test_chain_seed_avoids_collisions_for_bulky_host_and_two_donors(
    tmp_path: Path,
) -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = Manifest.new(run_id="test", stage="build-complexes")
    manifest.extend(
        [
            Artifact(
                artifact_id="tma",
                artifact_type="species",
                paths={"xyz": str(root / "examples/reaction_discovery/tma.xyz")},
                data={"species_id": "tma", "charge": 0, "multiplicity": 1},
            ),
            Artifact(
                artifact_id="hf",
                artifact_type="species",
                paths={"xyz": str(root / "examples/reaction_discovery/hf.xyz")},
                data={"species_id": "hf", "charge": 0, "multiplicity": 1},
            ),
        ]
    )

    result = BuildComplexesStage().run(
        manifest,
        {
            "compositions": [
                {
                    "composition_id": "bulky_host_hf2",
                    "components": [
                        {"selector": {"artifact_id": "tma"}, "count": 1},
                        {"selector": {"artifact_id": "hf"}, "count": 2},
                    ],
                    "max_seeds": 1,
                    "placement_topologies": ["chain"],
                    "crest_nci": False,
                }
            ]
        },
        StageContext(tmp_path / "out", "test", {}),
    )

    complex_species = next(
        item
        for item in result.latest_artifacts("species")
        if item.data.get("composition_id") == "bulky_host_hf2"
    )
    placements = complex_species.qc["placements"]
    assert complex_species.qc["component_count"] == 3
    assert len(placements) == 2
    assert all(item["minimum_interfragment_distance_A"] >= 0.55 for item in placements)


def test_global_trial_budget_is_distributed_across_nci_sources(
    tmp_path: Path,
) -> None:
    manifest = Manifest.new(run_id="round_robin", stage="preopt")
    for index, energy in enumerate((0.0, 0.5)):
        xyz = write_xyz(
            XYZ(
                ["H", "C", "N"],
                np.asarray([[0.0, 0.0, 0.0], [1.1, 0.0, 0.0], [2.2, 0.0, 0.0]]),
            ),
            tmp_path / f"nci_{index}.xyz",
        )
        manifest.add_artifact(
            Artifact(
                artifact_id=f"preopt_nci_{index}",
                artifact_type="species_preopt",
                paths={"xyz": str(xyz)},
                data={
                    "species_id": f"nci_{index}",
                    "state": "encounter_complex",
                    "xyz_path": str(xyz),
                    "crest_nci_conformer_id": f"conf_{index}",
                    "relative_energy_kcal_mol": energy,
                },
            )
        )

    result = GenerateReactionsStage().run(
        manifest,
        {
            "states": ["encounter_complex"],
            "required_source_data_keys": ["crest_nci_conformer_id"],
            "automatic": False,
            "explicit_trials": [
                {"associations": [[0, 2]], "dissociations": [[0, 1]]}
            ],
            "max_trials_per_species": 4,
            "max_total_trials": 2,
        },
        StageContext(tmp_path / "trials", "round_robin", {}),
    )
    sources = {
        trial.data["source_species_id"]
        for trial in result.latest_artifacts("reaction_trial")
        if trial.status.status == "success"
    }

    assert sources == {"nci_0", "nci_1"}

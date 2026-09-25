from pathlib import Path

from hfauto.chemistry.stoichiometry import molecular_surface_key
from hfauto.core.hashing import sha256_file
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import StageContext
from hfauto.stages.dft_minima import (
    DFTMinimaStage,
    _anchor_minimum_species_ids,
    _candidate_endpoint_species_ids,
    _next_minimum_workdir,
    _preferred_species_inputs,
    _reusable_minimum_calculation,
    _validated_minimum_provenance,
)


def test_dft_minima_can_explicitly_refine_an_optimized_endpoint() -> None:
    endpoint = Artifact(
        artifact_id="opt_reactant",
        artifact_type="species_optimized",
        data={"species_id": "reactant", "state": "reactant"},
    )
    manifest = Manifest.new(run_id="test", stage="dft-minima")
    manifest.add_artifact(endpoint)

    assert _preferred_species_inputs(manifest, {"reactant"}) == []
    assert _preferred_species_inputs(
        manifest,
        {"reactant"},
        include_optimized=True,
    ) == [endpoint]


def test_molecular_surface_key_includes_composition_charge_and_multiplicity(
    tmp_path: Path,
) -> None:
    xyz = tmp_path / "nh.xyz"
    xyz.write_text("2\nsurface\nN 0 0 0\nH 1 0 0\n", encoding="utf-8")

    def species(charge: int, multiplicity: int) -> Artifact:
        return Artifact(
            artifact_id=f"nh_{charge}_{multiplicity}",
            artifact_type="species",
            paths={"xyz": str(xyz)},
            data={
                "species_id": f"nh_{charge}_{multiplicity}",
                "xyz_path": str(xyz),
                "charge": charge,
                "multiplicity": multiplicity,
            },
        )

    assert molecular_surface_key(species(0, 1)) == "H1N1|q=0|m=1"
    assert molecular_surface_key(species(1, 1)) != molecular_surface_key(
        species(0, 1)
    )
    assert molecular_surface_key(species(0, 2)) != molecular_surface_key(
        species(0, 1)
    )


def test_dft_anchor_selection_is_bounded_per_state() -> None:
    species = [
        Artifact(
            artifact_id=f"nci_{index}",
            artifact_type="species_preopt",
            data={
                "species_id": f"nci_{index}",
                "state": "encounter_complex",
            },
        )
        for index in range(3)
    ]

    selected = _anchor_minimum_species_ids(
        species, {"encounter_complex": 2}
    )

    assert selected == {"nci_0", "nci_1"}


def test_dft_anchor_selection_skips_clear_geometry_duplicates(tmp_path: Path) -> None:
    species = []
    for index, distance in enumerate((0.75, 0.751, 0.95)):
        xyz = tmp_path / f"nci_{index}.xyz"
        xyz.write_text(
            f"2\nanchor\nH 0 0 0\nH {distance} 0 0\n",
            encoding="utf-8",
        )
        species.append(
            Artifact(
                artifact_id=f"nci_{index}",
                artifact_type="species_preopt",
                paths={"xyz": str(xyz)},
                data={
                    "species_id": f"nci_{index}",
                    "state": "encounter_complex",
                    "xyz_path": str(xyz),
                },
            )
        )

    selected = _anchor_minimum_species_ids(
        species,
        {"encounter_complex": 2},
        duplicate_rmsd_A=0.01,
    )

    assert selected == {"nci_0", "nci_2"}


def test_dft_anchor_budget_is_balanced_across_molecular_surfaces(
    tmp_path: Path,
) -> None:
    species = []
    for species_id, atoms in (
        ("h2_first", "H 0 0 0\nH 0.75 0 0"),
        ("h2_second", "H 0 0 0\nH 1.20 0 0"),
        ("h3_first", "H 0 0 0\nH 0.75 0 0\nH 1.50 0 0"),
    ):
        xyz = tmp_path / f"{species_id}.xyz"
        xyz.write_text(
            f"{len(atoms.splitlines())}\nanchor\n{atoms}\n",
            encoding="utf-8",
        )
        species.append(
            Artifact(
                artifact_id=species_id,
                artifact_type="species_preopt",
                paths={"xyz": str(xyz)},
                data={
                    "species_id": species_id,
                    "state": "encounter_complex",
                    "xyz_path": str(xyz),
                },
            )
        )

    selected = _anchor_minimum_species_ids(
        species,
        {"encounter_complex": 2},
    )

    assert selected == {"h2_first", "h3_first"}


def test_nwchem_minimum_provenance_is_bound_to_the_rendered_input(
    tmp_path,
) -> None:
    input_path = tmp_path / "nwchem.nw"
    input_path.write_text(
        "charge 0\n"
        "basis spherical\n  * library def2-svpd\nend\n"
        "dft\n"
        "  xc b3lyp\n"
        "  mult 1\n"
        "  grid xfine\n"
        "  convergence energy 1.0e-08\n"
        "  disp vdw 3\n"
        "end\n",
        encoding="utf-8",
    )
    calculation = Artifact(
        artifact_id="calc",
        artifact_type="calculation",
        paths={"input": str(input_path)},
        method={
            "engine": "nwchem",
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
            "required_program_version": "7.2.3",
            "grid": "xfine",
            "scf_energy_tolerance": 1.0e-8,
        },
        data={
            "program_version": "7.2.3",
            "resolved_charge": 0,
            "resolved_multiplicity": 1,
            "electron_count": 10,
        },
        qc={
            "real_qm_executed": True,
            "fallback_dummy": False,
            "minimum_accepted": True,
            "dispersion_applied": True,
        },
    )

    provenance, method_validated, numerical_validated = (
        _validated_minimum_provenance(calculation)
    )

    assert method_validated is False
    assert numerical_validated is True
    assert provenance["validated_input_method"]["accepted"] is False
    assert provenance["validated_input_method"]["reasons"] == [
        "nwchem_input_functional_does_not_match_calculation"
    ]


def test_dft_minima_can_be_limited_by_a_seed_selection(tmp_path: Path) -> None:
    xyz = tmp_path / "seed.xyz"
    xyz.write_text("1\nseed\nH 0.0 0.0 0.0\n", encoding="utf-8")
    manifest = Manifest.new(run_id="selected", stage="endpoint-seed-screen")
    for species_id in ("seed_a", "seed_b"):
        manifest.add_artifact(
            Artifact(
                artifact_id=species_id,
                artifact_type="species",
                paths={"xyz": str(xyz)},
                data={
                    "species_id": species_id,
                    "state": "endpoint_minimum_candidate",
                    "xyz_path": str(xyz),
                    "endpoint_seed_generation": 1,
                },
            )
        )
    manifest.add_artifact(
        Artifact(
            artifact_id="selection",
            artifact_type="endpoint_seed_selection",
            data={"selected_species_ids": ["seed_b"]},
        )
    )

    output = DFTMinimaStage().run(
        manifest,
        {
            "engine": "dummy",
            "states": ["endpoint_minimum_candidate"],
            "selection_artifact_type": "endpoint_seed_selection",
            "data_filters": {"endpoint_seed_generation": 1},
            "require_real_qm": False,
        },
        StageContext(
            out_dir=tmp_path / "minima",
            run_id="selected",
            global_config={},
        ),
    )
    calculations = output.latest_artifacts("calculation")

    assert [item.data["species_id"] for item in calculations] == ["seed_b"]
    assert output.metadata["dft_minima_selected_species_count"] == 1


def test_dft_minima_reuses_only_an_integrity_checked_exact_checkpoint(
    tmp_path: Path,
) -> None:
    xyz = tmp_path / "reactant.xyz"
    xyz.write_text("1\nreactant\nH 0.0 0.0 0.0\n", encoding="utf-8")
    output = tmp_path / "nwchem.out"
    output.write_text("complete", encoding="utf-8")
    species = Artifact(
        artifact_id="source_geometry",
        artifact_type="species",
        paths={"xyz": str(xyz)},
        data={
            "species_id": "reactant",
            "state": "reactant",
            "xyz_path": str(xyz),
        },
    )
    method = {"functional": "pbe0", "basis": "def2-svpd"}
    calculation = Artifact(
        artifact_id="calc_reactant",
        artifact_type="calculation",
        method={"engine": "nwchem", **method},
        data={
            "species_id": "reactant",
            "source_geometry_artifact_id": "source_geometry",
        },
        qc={
            "minimum_accepted": True,
            "real_qm_executed": True,
            "fallback_dummy": False,
        },
        provenance={
            "input_xyz_sha256": sha256_file(xyz),
            "command": {"returncode": 0, "timed_out": False},
        },
    )
    optimized = Artifact(
        artifact_id="opt_reactant",
        artifact_type="species_optimized",
        qc={"dft_calc_id": calculation.artifact_id},
        provenance={
            "validated_source_integrity": {
                "accepted": True,
                "reasons": [],
                "files": {
                    "output": {
                        "path": str(output),
                        "sha256": sha256_file(output),
                    }
                },
            }
        },
    )
    checkpoint = Manifest.new(run_id="resume", stage="dft-minima")
    checkpoint.extend([calculation, optimized])

    reused = _reusable_minimum_calculation(
        checkpoint,
        species,
        engine="nwchem",
        method=method,
        require_real_qm=True,
    )

    assert reused is not None
    assert reused.provenance["reused_from_checkpoint"] is True

    output.write_text("tampered", encoding="utf-8")
    assert (
        _reusable_minimum_calculation(
            checkpoint,
            species,
            engine="nwchem",
            method=method,
            require_real_qm=True,
        )
        is None
    )


def test_dft_minima_retry_uses_a_new_directory(tmp_path: Path) -> None:
    base = tmp_path / "minimum"
    assert _next_minimum_workdir(base) == base
    base.mkdir()
    (base / "nwchem.nw").write_text("incomplete", encoding="utf-8")
    assert _next_minimum_workdir(base) == base / "attempt_02"
    (base / "attempt_02").mkdir()
    (base / "attempt_02" / "nwchem.nw").write_text(
        "failed retry",
        encoding="utf-8",
    )
    assert _next_minimum_workdir(base) == base / "attempt_03"


def test_candidate_endpoint_selection_keeps_the_exact_trial_source(
    tmp_path: Path,
) -> None:
    manifest = Manifest.new(run_id="pairs", stage="explore-reactions")
    xyz = tmp_path / "endpoint.xyz"
    xyz.write_text(
        "3\nendpoint\nH 0 0 0\nC 1.1 0 0\nN 2.2 0 0\n",
        encoding="utf-8",
    )
    for species_id, state, energy in (
        ("unrelated_low_nci", "encounter_complex", -20.0),
        ("trial_source", "encounter_complex", -19.0),
        ("product_low", "product_candidate", -21.0),
        ("product_high", "product_candidate", -18.0),
    ):
        manifest.add_artifact(
            Artifact(
                artifact_id=species_id,
                artifact_type="species",
                paths={"xyz": str(xyz)},
                data={
                    "species_id": species_id,
                    "state": state,
                    "xyz_path": str(xyz),
                    "discovery_electronic_energy_hartree": energy,
                },
            )
        )
    manifest.extend(
        [
            Artifact(
                artifact_id="candidate_low",
                artifact_type="reaction_candidate",
                data={
                    "candidate_id": "candidate_low",
                    "reactant_species_id": "trial_source",
                    "product_species_id": "product_low",
                },
                qc={"structural_change_detected": True},
            ),
            Artifact(
                artifact_id="candidate_high",
                artifact_type="reaction_candidate",
                data={
                    "candidate_id": "candidate_high",
                    "reactant_species_id": "unrelated_low_nci",
                    "product_species_id": "product_high",
                },
                qc={"structural_change_detected": True},
            ),
        ]
    )

    species_ids, candidate_ids = _candidate_endpoint_species_ids(manifest, 1)

    assert candidate_ids == ["candidate_low"]
    assert species_ids == {"trial_source", "product_low"}


def test_candidate_endpoint_selection_removes_only_clear_low_level_duplicates(
    tmp_path: Path,
) -> None:
    manifest = Manifest.new(run_id="diverse", stage="explore-reactions")
    paths = []
    for index, terminal_x in enumerate((2.20, 2.201, 4.0)):
        path = tmp_path / f"product_{index}.xyz"
        path.write_text(
            f"3\nproduct\nH 0 0 0\nC 1.1 0 0\nN {terminal_x} 0 0\n",
            encoding="utf-8",
        )
        paths.append(path)
        manifest.add_artifact(
            Artifact(
                artifact_id=f"product_{index}",
                artifact_type="species",
                paths={"xyz": str(path)},
                data={
                    "species_id": f"product_{index}",
                    "state": "product_candidate",
                    "xyz_path": str(path),
                    "discovery_electronic_energy_hartree": -10.0 + index * 0.01,
                },
            )
        )
        manifest.add_artifact(
            Artifact(
                artifact_id=f"source_{index}",
                artifact_type="species",
                paths={"xyz": str(path)},
                data={
                    "species_id": f"source_{index}",
                    "state": "encounter_complex",
                    "xyz_path": str(path),
                },
            )
        )
        manifest.add_artifact(
            Artifact(
                artifact_id=f"candidate_{index}",
                artifact_type="reaction_candidate",
                data={
                    "candidate_id": f"candidate_{index}",
                    "reactant_species_id": f"source_{index}",
                    "product_species_id": f"product_{index}",
                },
                qc={"structural_change_detected": True},
            )
        )

    species_ids, candidate_ids = _candidate_endpoint_species_ids(
        manifest,
        3,
        duplicate_rmsd_A=0.05,
    )

    assert candidate_ids == ["candidate_0", "candidate_2"]
    assert species_ids == {
        "source_0",
        "product_0",
        "source_2",
        "product_2",
    }


def test_candidate_budget_never_compares_absolute_energy_across_formulas(
    tmp_path: Path,
) -> None:
    manifest = Manifest.new(run_id="balanced", stage="explore-reactions")
    definitions = (
        ("h2_low", "H 0 0 0\nH 0.75 0 0", -100.0),
        ("h2_high", "H 0 0 0\nH 1.20 0 0", -90.0),
        ("h3", "H 0 0 0\nH 0.75 0 0\nH 1.50 0 0", -1.0),
    )
    for species_id, atoms, energy in definitions:
        xyz = tmp_path / f"{species_id}.xyz"
        xyz.write_text(
            f"{len(atoms.splitlines())}\nproduct\n{atoms}\n",
            encoding="utf-8",
        )
        for role in ("source", "product"):
            artifact_id = f"{role}_{species_id}"
            manifest.add_artifact(
                Artifact(
                    artifact_id=artifact_id,
                    artifact_type="species",
                    paths={"xyz": str(xyz)},
                    data={
                        "species_id": artifact_id,
                        "state": (
                            "encounter_complex"
                            if role == "source"
                            else "product_candidate"
                        ),
                        "xyz_path": str(xyz),
                        "discovery_electronic_energy_hartree": energy,
                    },
                )
            )
        manifest.add_artifact(
            Artifact(
                artifact_id=f"candidate_{species_id}",
                artifact_type="reaction_candidate",
                data={
                    "candidate_id": f"candidate_{species_id}",
                    "reactant_species_id": f"source_{species_id}",
                    "product_species_id": f"product_{species_id}",
                },
                qc={"structural_change_detected": True},
            )
        )

    species_ids, candidate_ids = _candidate_endpoint_species_ids(manifest, 2)

    assert candidate_ids == ["candidate_h2_low", "candidate_h3"]
    assert species_ids == {
        "source_h2_low",
        "product_h2_low",
        "source_h3",
        "product_h3",
    }

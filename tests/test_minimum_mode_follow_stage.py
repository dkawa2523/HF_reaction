from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from hfauto.chemistry.xyz import XYZ, read_xyz, write_xyz
from hfauto.core.hashing import sha256_file
from hfauto.core.schemas.artifact import Artifact, ArtifactStatus
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import StageContext
from hfauto.stages.minimum_mode_assess import MinimumModeAssessStage
from hfauto.stages.minimum_mode_follow import MinimumModeFollowStage


def _manifest_with_one_mode_source(tmp_path: Path) -> Manifest:
    source_xyz = tmp_path / "source.xyz"
    write_xyz(
        XYZ(
            symbols=["N", "H", "F"],
            coords=np.array(
                [[0.0, 0.0, 0.0], [1.0, 0.1, 0.0], [2.0, 0.0, 0.1]]
            ),
            comment="one-mode source",
        ),
        source_xyz,
    )
    input_path = tmp_path / "nwchem.nw"
    input_path.write_text(
        """geometry units angstrom noautoz
 N 0.0 0.0 0.0
 H 1.0 0.1 0.0
 F 2.0 0.0 0.1
end
charge 0
basis
 * library def2-svpd
end
dft
 mult 2
 xc pbe0
 grid xfine
 convergence energy 1.0e-8
 disp vdw 3
end
driver
 tight
 maxiter 100
end
task dft optimize
task dft frequencies
""",
        encoding="utf-8",
    )
    output_path = tmp_path / "nwchem.out"
    stderr_path = tmp_path / "nwchem.err"
    output_path.write_text("real output bound by hash\n", encoding="utf-8")
    stderr_path.write_text("", encoding="utf-8")
    mode = [[0.1, 0.0, 0.0], [0.0, -0.2, 0.0], [-0.1, 0.1, 0.0]]

    source = Artifact(
        artifact_id="opt_spc_source",
        artifact_type="species_optimized",
        paths={"xyz": str(source_xyz)},
        data={
            "species_id": "spc_source",
            "state": "reactant",
            "charge": 0,
            "multiplicity": 2,
            "electron_count": 17,
            "xyz_path": str(source_xyz),
        },
    )
    calculation = Artifact(
        artifact_id="calc_one_mode",
        artifact_type="calculation",
        parents=[source.artifact_id],
        paths={
            "input": str(input_path),
            "output": str(output_path),
            "stderr": str(stderr_path),
            "final_xyz": str(source_xyz),
        },
        data={
            "species_id": "spc_source",
            "source_geometry_artifact_id": source.artifact_id,
            "task": "opt_freq",
            "calculation_level": "nwchem_real",
            "resolved_charge": 0,
            "resolved_multiplicity": 2,
            "electron_count": 17,
            "program_version": "7.2.3",
            "dft_d3_applied": True,
            "frequency_count_complete": True,
            "raw_frequency_count": 9,
            "expected_raw_frequency_count": 9,
            "n_imag": 1,
            "imag_freq_cm1": -42.0,
            "imaginary_mode_displacements": mode,
            "projected_imaginary_mode": {
                "mode_number": 1,
                "frequency_cm1": -42.0,
                "cartesian_displacements": mode,
                "coordinate_convention": (
                    "NWChem projected normal-mode eigenvector in Cartesian "
                    "coordinates"
                ),
                "component_units": "amu^-1/2",
                "mass_weighting_handling": (
                    "already transformed to Cartesian coordinates"
                ),
            },
        },
        method={
            "stage": "dft-minima",
            "engine": "nwchem",
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
            "grid": "xfine",
            "scf_energy_tolerance": 1.0e-8,
            "optimization_convergence": "tight",
            "geometry_maxiter": 100,
            "charge": 0,
            "multiplicity": 2,
        },
        qc={
            "real_qm_executed": True,
            "fallback_dummy": False,
            "scf_converged": True,
            "geometry_converged": True,
            "normal_termination": True,
            "n_imag": 1,
            "minimum_accepted": False,
        },
        provenance={
            "created_by": "NWChemEngine",
            "electronic_state": {
                "charge": 0,
                "multiplicity": 2,
                "electron_count": 17,
            },
            "command": {
                "returncode": 0,
                "timed_out": False,
                "stdout_path": str(output_path),
                "stderr_path": str(stderr_path),
            },
        },
        status=ArtifactStatus(status="success"),
    )
    reaction = Artifact(
        artifact_id="rxn_source_to_product",
        artifact_type="reaction",
        data={
            "reaction_id": "rxn_source_to_product",
            "reactant_species_id": "spc_source",
            "product_species_id": "spc_product",
        },
    )
    manifest = Manifest.new(run_id="mode_follow_test", stage="dft-minima")
    manifest.extend([source, reaction, calculation])
    return manifest


def _context(tmp_path: Path) -> StageContext:
    return StageContext(
        out_dir=tmp_path / "mode_follow",
        run_id="mode_follow_test",
        global_config={"mode": "production"},
    )


def test_stage_generates_only_two_bounded_unpromoted_candidates(
    tmp_path: Path,
) -> None:
    manifest = _manifest_with_one_mode_source(tmp_path)
    result = MinimumModeFollowStage().run(
        manifest,
        {
            "source_calculation_ids": ["calc_one_mode"],
            "maximum_atom_displacement_A": 0.10,
        },
        _context(tmp_path),
    )

    plan = result.latest_artifacts("minimum_mode_following_plan")[0]
    seeds = [
        artifact
        for artifact in result.latest_artifacts("species")
        if artifact.data.get("mode_following_generation") == 1
    ]
    assert plan.status.status == "success"
    assert plan.data["accepted"] is True
    assert plan.qc["source_promoted_to_minimum"] is False
    assert plan.qc["source_promoted_to_transition_state"] is False
    assert len(seeds) == 2
    assert {seed.data["mode_following_direction"] for seed in seeds} == {
        "plus",
        "minus",
    }
    assert len(result.latest_artifacts("species_optimized")) == 1

    calculation = manifest.find("calc_one_mode")
    assert calculation is not None
    source = read_xyz(calculation.paths["final_xyz"])
    geometries = {
        seed.data["mode_following_direction"]: read_xyz(seed.paths["xyz"])
        for seed in seeds
    }
    for geometry in geometries.values():
        atom_displacements = np.linalg.norm(
            geometry.coords - source.coords, axis=1
        )
        np.testing.assert_allclose(
            atom_displacements.max(), 0.10, atol=1.0e-8
        )
    np.testing.assert_allclose(
        (geometries["plus"].coords + geometries["minus"].coords) / 2.0,
        source.coords,
        atol=1.0e-8,
    )
    assert geometries["plus"].coords.tolist() != geometries[
        "minus"
    ].coords.tolist()


def test_stage_fails_closed_when_raw_command_evidence_is_invalid(
    tmp_path: Path,
) -> None:
    manifest = _manifest_with_one_mode_source(tmp_path)
    calculation = manifest.find("calc_one_mode")
    assert calculation is not None
    calculation.provenance["command"]["returncode"] = 1

    result = MinimumModeFollowStage().run(
        manifest,
        {"source_calculation_ids": ["calc_one_mode"]},
        _context(tmp_path),
    )

    plan = result.latest_artifacts("minimum_mode_following_plan")[0]
    assert plan.status.status == "failed"
    assert plan.data["seed_artifact_ids"] == []
    assert "source_command_returncode_is_not_zero" in plan.data[
        "source_eligibility"
    ]["reasons"]
    assert result.metadata["minimum_mode_following_seed_count"] == 0


def test_assessment_resolves_local_same_basin_without_claiming_reaction(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    manifest = _manifest_with_one_mode_source(tmp_path)
    generated = MinimumModeFollowStage().run(
        manifest,
        {"source_calculation_ids": ["calc_one_mode"]},
        _context(tmp_path),
    )
    seeds = [
        artifact
        for artifact in generated.latest_artifacts("species")
        if artifact.data.get("mode_following_generation") == 1
    ]
    common_xyz = tmp_path / "common_minimum.xyz"
    source_calculation = generated.find("calc_one_mode")
    assert source_calculation is not None
    write_xyz(
        read_xyz(source_calculation.paths["final_xyz"]),
        common_xyz,
    )
    for index, seed in enumerate(seeds):
        calc_id = f"calc_successor_{index}"
        calculation = Artifact(
            artifact_id=calc_id,
            artifact_type="calculation",
            parents=[seed.artifact_id],
            paths={"final_xyz": str(common_xyz)},
            data={
                "source_geometry_artifact_id": seed.artifact_id,
                "input_xyz_sha256": sha256_file(seed.paths["xyz"]),
                "electronic_energy_hartree": -100.0,
                "frequency_count_complete": True,
                "n_imag": 0,
            },
            qc={"real_qm_executed": True, "minimum_accepted": True},
        )
        minimum = Artifact(
            artifact_id=f"opt_{seed.artifact_id}",
            artifact_type="species_optimized",
            parents=[seed.artifact_id, calc_id],
            paths={"xyz": str(common_xyz)},
            data={
                **seed.data,
                "xyz_path": str(common_xyz),
                "resolved_charge": 0,
                "resolved_multiplicity": 2,
                "electron_count": 17,
                "dft_minima": {"calc_id": calc_id},
            },
            qc={
                "dft_calc_id": calc_id,
                "minimum_accepted": True,
                "is_minimum": True,
                "n_imag": 0,
                "real_qm_executed": True,
                "state_method_evidence_validated": True,
            },
        )
        generated.extend([calculation, minimum])

    monkeypatch.setattr(
        "hfauto.stages.minimum_mode_assess.evaluate_endpoint_pair_lineage",
        lambda _plus, _minus: {"accepted": True, "reasons": []},
    )
    result = MinimumModeAssessStage().run(
        generated,
        {},
        StageContext(
            out_dir=tmp_path / "assess",
            run_id="mode_follow_test",
            global_config={"mode": "production"},
        ),
    )

    assessment = result.latest_artifacts(
        "minimum_mode_following_assessment"
    )[0]
    assert assessment.status.status == "success"
    assert assessment.data["outcome"] == "same_basin_descents"
    assert assessment.data["local_basin_conclusion_supported"] is True
    assert assessment.data["reaction_classification_supported"] is False
    assert assessment.qc["source_promoted_to_transition_state"] is False
    assert assessment.qc["bidirectional_irc_validated"] is False


def test_assessment_rejects_source_mutation_after_seed_generation(
    tmp_path: Path,
) -> None:
    manifest = _manifest_with_one_mode_source(tmp_path)
    generated = MinimumModeFollowStage().run(
        manifest,
        {"source_calculation_ids": ["calc_one_mode"]},
        _context(tmp_path),
    )
    source = generated.find("calc_one_mode")
    assert source is not None
    Path(source.paths["output"]).write_text(
        "output changed after seed generation\n", encoding="utf-8"
    )

    result = MinimumModeAssessStage().run(
        generated,
        {},
        StageContext(
            out_dir=tmp_path / "changed_assess",
            run_id="mode_follow_test",
            global_config={"mode": "production"},
        ),
    )

    assessment = result.latest_artifacts(
        "minimum_mode_following_assessment"
    )[0]
    assert assessment.status.status == "partial"
    assert "source_evidence_changed_after_seed_generation" in assessment.data[
        "source_stationary_point_reaudit"
    ]["reasons"]

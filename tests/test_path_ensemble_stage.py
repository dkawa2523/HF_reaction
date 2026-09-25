from pathlib import Path

from hfauto.chemistry.method_lineage import (
    evaluate_endpoint_pair_lineage,
    evaluate_endpoint_path_numerical_lineage,
    evaluate_endpoint_source_integrity,
)
from hfauto.chemistry.path_convergence import path_method_signature
from hfauto.chemistry.reaction_profile import analyze_reaction_path
from hfauto.core.hashing import sha256_file
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import StageContext
from hfauto.stages.path_ensemble_assess import PathEnsembleAssessStage


def _endpoint(
    identifier: str,
    state: str,
    evidence_dir: Path | None = None,
) -> Artifact:
    validated_method = {
        "engine": "nwchem",
        "functional": "pbe0",
        "basis": "def2-svpd",
        "disp_vdw": 3,
        "program_version": "7.2.3",
    }
    electronic_state = {"charge": 0, "multiplicity": 1, "electron_count": 10}
    numerical_input = {}
    source_integrity = {}
    if evidence_dir is not None:
        input_path = evidence_dir / f"{identifier}.nw"
        input_path.write_text(
            "grid xfine\nconvergence energy 1.0e-08\n",
            encoding="utf-8",
        )
        numerical_input = {
            "path": str(input_path),
            "sha256": sha256_file(input_path),
            "grid": "xfine",
            "scf_energy_tolerance": 1.0e-8,
        }
        source_paths = {
            "input": input_path,
            "output": evidence_dir / f"{identifier}.out",
            "stderr": evidence_dir / f"{identifier}.err",
            "final_xyz": evidence_dir / f"{identifier}.xyz",
        }
        source_paths["output"].write_text("completed\n", encoding="utf-8")
        source_paths["stderr"].write_text("", encoding="utf-8")
        source_paths["final_xyz"].write_text(
            "1\nminimum\nH 0.0 0.0 0.0\n", encoding="utf-8"
        )
        source_integrity = {
            "accepted": True,
            "reasons": [],
            "files": {
                key: {"path": str(path), "sha256": sha256_file(path)}
                for key, path in source_paths.items()
            },
        }
    return Artifact(
        artifact_id=identifier,
        artifact_type="species_optimized",
        data={
            "species_id": identifier,
            "state": state,
            "resolved_charge": 0,
            "resolved_multiplicity": 1,
            "electron_count": 10,
        },
        method={
            **validated_method,
            "required_program_version": "7.2.3",
            "grid": "xfine",
            "scf_energy_tolerance": 1.0e-8,
        },
        qc={
            "state_method_evidence_validated": True,
            "numerical_method_evidence_validated": True,
        },
        provenance={
            "validated_qm_method": validated_method,
            "validated_electronic_state": electronic_state,
            "validated_numerical_settings": {
                "grid": "xfine",
                "scf_energy_tolerance": 1.0e-8,
            },
            "validated_numerical_input": numerical_input,
            "validated_source_integrity": source_integrity,
        },
    )


def _source_path(tmp_path: Path, identifier: str, image_count: int) -> Artifact:
    input_path = tmp_path / f"{identifier}.nw"
    output_path = tmp_path / f"{identifier}.out"
    path_xyz = tmp_path / f"{identifier}.xyz"
    input_path.write_text(
        "charge 0\n"
        "basis spherical\n"
        "  * library def2-svpd\n"
        "end\n"
        "dft\n"
        "  xc pbe0\n"
        "  mult 1\n"
        "  grid xfine\n"
        "  convergence energy 1.0e-08\n"
        "  disp vdw 3\n"
        "end\n"
        f"string\n  nbeads {image_count}\n"
        "  maxiter 150\n"
        "  stepsize 0.01\n"
        "  nhist 10\n"
        "  interpol 3\n"
        "  tol 0.0002\n"
        "end\n",
        encoding="utf-8",
    )
    output_path.write_text(f"raw-output-{image_count}\n", encoding="utf-8")
    path_xyz.write_text(
        "".join(
            "2\n"
            f"image={index}\n"
            "H 0.0 0.0 0.0\n"
            f"F 0.0 0.0 {1.0 + index / (image_count - 1):.8f}\n"
            for index in range(image_count)
        ),
        encoding="utf-8",
    )
    return Artifact(
        artifact_id=identifier,
        artifact_type="reaction_path",
        paths={
            "input": str(input_path),
            "output": str(output_path),
            "path_xyz": str(path_xyz),
        },
        data={"reaction_id": "rxn", "program_version": "7.2.3"},
        method={
            "engine": "nwchem",
            "backend": "nwchem_string",
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
            "program_version": "7.2.3",
            "required_program_version": "7.2.3",
            "charge": 0,
            "multiplicity": 1,
            "settings": {
                "grid": "xfine",
                "scf_energy_tolerance": 1.0e-8,
                "string_tolerance": 0.0002,
                "string_stepsize": 0.01,
                "string_maxiter": 150,
                "string_nhist": 10,
                "string_interpol": 3,
                "path_image_count": image_count,
            },
        },
        qc={
            "real_path_executed": True,
            "path_geometry_validated": True,
            "method_evidence_validated": True,
            "fallback_dummy": False,
        },
        provenance={
            "input_sha256": sha256_file(input_path),
            "output_sha256": sha256_file(output_path),
            "path_sha256": sha256_file(path_xyz),
        },
    )


def test_endpoint_numerical_settings_must_match_the_path_resolution() -> None:
    reactant = _endpoint("reactant", "reactant")
    product = _endpoint("product", "product")
    reactant.method["grid"] = "fine"
    reactant.method["scf_energy_tolerance"] = 1.0e-7
    reactant.provenance["validated_numerical_settings"] = {
        "grid": "fine",
        "scf_energy_tolerance": 1.0e-7,
    }

    result = evaluate_endpoint_path_numerical_lineage(
        reactant,
        product,
        {
            "engine": "nwchem",
            "grid": "xfine",
            "scf_energy_tolerance": 1.0e-8,
        },
    )

    assert result["accepted"] is False
    assert "reactant_grid_does_not_match_path_setting" in result["reasons"]


def test_endpoint_raw_evidence_is_rehashed_before_path_execution(
    tmp_path: Path,
) -> None:
    endpoint = _endpoint("reactant", "reactant", tmp_path)

    assert evaluate_endpoint_source_integrity(endpoint)["accepted"] is True
    output = endpoint.provenance["validated_source_integrity"]["files"][
        "output"
    ]["path"]
    Path(output).write_text("changed after validation\n", encoding="utf-8")

    audit = evaluate_endpoint_source_integrity(endpoint)
    assert audit["accepted"] is False
    assert "validated_source_output_hash_mismatch" in audit["reasons"]


def test_endpoint_pair_lineage_rejects_mixed_numerical_surfaces(
    tmp_path: Path,
) -> None:
    reactant = _endpoint("reactant", "reactant", tmp_path)
    product = _endpoint("product", "product", tmp_path)
    product.provenance["validated_numerical_settings"] = {
        "grid": "fine",
        "scf_energy_tolerance": 1.0e-7,
    }

    lineage = evaluate_endpoint_pair_lineage(reactant, product)

    assert lineage["accepted"] is False
    assert "endpoint_numerical_settings_do_not_match" in lineage["reasons"]


def _member(source: Artifact, identifier: str, image_count: int) -> Artifact:
    energies = [
        -10.0 - 0.4 * index / (image_count - 1)
        for index in range(image_count)
    ]
    path = analyze_reaction_path(
        reaction_id="rxn",
        engine="nwchem_string",
        comments=[f"energy_hartree={energy}" for energy in energies],
        converged=True,
        barrier_threshold_kcal_mol=0.05,
    )
    return Artifact(
        artifact_id=identifier,
        artifact_type="path_ensemble_member",
        parents=[source.artifact_id],
        data={
            "reaction_id": "rxn",
            "study_id": "study",
            "source_path_artifact_id": source.artifact_id,
            "refines_path_attempt_id": "coarse_attempt",
            "image_count": image_count,
            "path": path.model_dump(),
            "method_signature": path_method_signature(source),
            "input_sha256": source.provenance["input_sha256"],
            "output_sha256": source.provenance["output_sha256"],
        },
        qc={
            "include_in_assessment": True,
            "real_path_executed": True,
            "path_geometry_validated": True,
            "method_evidence_validated": True,
            "fallback_dummy": False,
        },
    )


def test_assessment_reaudits_raw_inputs_and_publishes_lineage(
    tmp_path: Path,
) -> None:
    source_5 = _source_path(tmp_path, "path_5", 5)
    source_9 = _source_path(tmp_path, "path_9", 9)
    member_5 = _member(source_5, "member_5", 5)
    member_9 = _member(source_9, "member_9", 9)
    member_5.data["input_sha256"] = "untrusted-serialized-hash"
    manifest = Manifest.new(run_id="run", stage="path-ensemble")
    reaction = Artifact(
        artifact_id="rxn",
        artifact_type="reaction",
        data={
            "reaction_id": "rxn",
            "reactant_species_id": "reactant",
            "product_species_id": "product",
            "basin_assessment": {
                "status": "distinct_basin",
                "accepted": True,
            },
        },
        qc={"distinct_registry_basins": True},
    )
    manifest.extend(
        [
            reaction,
            _endpoint("reactant", "reactant", tmp_path),
            _endpoint("product", "product", tmp_path),
            source_5,
            source_9,
            member_5,
            member_9,
        ]
    )

    output = PathEnsembleAssessStage().run(
        manifest,
        {
            "minimum_resolutions": 2,
            "profile_energy_tolerance_kcal_mol": 0.05,
            "supersedes_path_attempt_ids": ["older_attempt"],
        },
        StageContext(
            out_dir=tmp_path / "assessment",
            run_id="run",
            global_config={"mode": "production"},
        ),
    )
    assessment = output.latest_artifacts("path_ensemble_assessment")[-1]

    assert assessment.status.status == "success"
    assert assessment.data["barrierless_at_resolution"] is True
    assert assessment.data["path_barrier_resolution_kcal_mol"] == 0.05
    assert assessment.data["next_action"] == (
        "classify_effectively_barrierless_at_resolution"
    )
    assert assessment.data["supersedes_path_attempt_ids"] == [
        "coarse_attempt",
        "older_attempt",
    ]
    assert all(
        audit["rendered_input"]["accepted"] is True
        and audit["endpoint_method_lineage"]["accepted"] is True
        for audit in assessment.data["member_audits"].values()
    )
    assert output.find("member_5").data["input_sha256"] == (
        source_5.provenance["input_sha256"]
    )

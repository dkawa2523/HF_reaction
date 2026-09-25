from hfauto.chemistry.path_convergence import assess_path_ensemble
from hfauto.chemistry.reaction_profile import analyze_reaction_path
from hfauto.core.schemas.artifact import Artifact


def _member(
    identifier: str,
    energies: list[float],
    *,
    input_hash: str,
    output_hash: str,
    path_hash: str | None = None,
    method_signature: dict | None = None,
    coordinate: list[float] | None = None,
) -> Artifact:
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
        data={
            "reaction_id": "rxn",
            "path": path.model_dump(),
            "method_signature": method_signature
            or {"engine": "nwchem", "functional": "pbe0", "basis": "def2-svpd"},
            "input_sha256": input_hash,
            "output_sha256": output_hash,
            "path_sha256": path_hash or f"path-{identifier}",
            "normalized_path_coordinate": coordinate
            or [index / (len(energies) - 1) for index in range(len(energies))],
        },
        qc={
            "real_path_executed": True,
            "path_geometry_validated": True,
            "method_evidence_validated": True,
            "input_method_evidence_validated": True,
            "endpoint_method_lineage_validated": True,
            "path_coordinate_validated": True,
            "fallback_dummy": False,
        },
    )


def test_two_matching_discretizations_support_barrierless_profile() -> None:
    members = [
        _member(
            "images_5",
            [-10.0, -10.1, -10.2, -10.3, -10.4],
            input_hash="input5",
            output_hash="output5",
        ),
        _member(
            "images_9",
            [-10.0, -10.05, -10.1, -10.15, -10.2, -10.25, -10.3, -10.35, -10.4],
            input_hash="input9",
            output_hash="output9",
        ),
    ]

    result = assess_path_ensemble(members)

    assert result["accepted"] is True
    assert result["barrierless_at_resolution"] is True
    assert result["image_counts"] == [5, 9]


def test_same_calculation_cannot_count_as_two_resolutions() -> None:
    members = [
        _member(
            "copy_a",
            [-10.0, -10.1, -10.2],
            input_hash="same-input",
            output_hash="same-output",
        ),
        _member(
            "copy_b",
            [-10.0, -10.1, -10.2],
            input_hash="same-input",
            output_hash="same-output",
        ),
    ]

    result = assess_path_ensemble(members)

    assert result["accepted"] is False
    assert "independent_image_discretizations_insufficient" in result["reasons"]


def test_different_quantum_methods_cannot_be_mixed() -> None:
    members = [
        _member(
            "pbe0",
            [-10.0, -10.1, -10.2],
            input_hash="input3",
            output_hash="output3",
        ),
        _member(
            "b3lyp",
            [-10.0, -10.05, -10.1, -10.15, -10.2],
            input_hash="input5",
            output_hash="output5",
            method_signature={
                "engine": "nwchem",
                "functional": "b3lyp",
                "basis": "def2-svpd",
            },
        ),
    ]

    result = assess_path_ensemble(members)

    assert result["accepted"] is False
    assert "path_methods_or_numerical_settings_differ" in result["reasons"]


def test_matching_class_but_unconverged_profiles_are_rejected() -> None:
    members = [
        _member(
            "shallow",
            [-10.0, -10.1, -10.2],
            input_hash="input3",
            output_hash="output3",
        ),
        _member(
            "curved",
            [-10.0, -10.01, -10.05, -10.15, -10.2],
            input_hash="input5",
            output_hash="output5",
        ),
    ]

    result = assess_path_ensemble(members, profile_energy_tolerance_kcal_mol=0.05)

    assert result["accepted"] is False
    assert "path_energy_profiles_not_converged" in result["reasons"]


def test_profiles_are_compared_by_geometry_arc_length_not_image_index() -> None:
    members = [
        _member(
            "three_images",
            [-10.0, -10.2, -10.4],
            input_hash="input3",
            output_hash="output3",
            coordinate=[0.0, 0.2, 1.0],
        ),
        _member(
            "five_images",
            [-10.0, -10.2, -10.25, -10.3, -10.4],
            input_hash="input5",
            output_hash="output5",
            coordinate=[0.0, 0.2, 0.4, 0.6, 1.0],
        ),
    ]

    result = assess_path_ensemble(
        members, profile_energy_tolerance_kcal_mol=0.001
    )

    assert result["accepted"] is True
    assert result["profile_coordinate"].startswith("normalized_cumulative")


def test_missing_members_report_completion_instead_of_a_method_mismatch() -> None:
    result = assess_path_ensemble([])

    assert result["accepted"] is False
    assert result["next_action"] == "complete_missing_path_ensemble_members"
    assert "eligible_path_count_insufficient" in result["reasons"]
    assert "path_methods_or_numerical_settings_differ" not in result["reasons"]

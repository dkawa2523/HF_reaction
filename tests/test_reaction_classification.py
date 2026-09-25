from pathlib import Path

import numpy as np

from hfauto.chemistry.reaction_classification import classify_reaction
from hfauto.chemistry.reaction_profile import analyze_reaction_path
from hfauto.chemistry.xyz import XYZ, write_xyz
from hfauto.chemistry.xyz_trajectory import write_xyz_trajectory
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import StageContext
from hfauto.stages.path_intermediates import PathIntermediatesStage
from hfauto.workflow.reaction_evidence import (
    latest_minimum_attempts_by_species,
)


def _reaction() -> Artifact:
    return Artifact(
        artifact_id="rxn",
        artifact_type="reaction",
        data={
            "reaction_id": "rxn",
            "reactant_species_id": "reactant",
            "product_species_id": "product",
        },
    )


def _basin(kind: str) -> dict:
    return {
        "accepted": True,
        "status": kind,
        "same_basin": kind == "same_basin",
        "distinct_basin": kind == "distinct_basin",
    }


def _attempt(identifier: str, energies: list[float]) -> Artifact:
    path = analyze_reaction_path(
        reaction_id="rxn",
        engine="test_path",
        comments=[f"energy_hartree={energy}" for energy in energies],
        converged=True,
        barrier_threshold_kcal_mol=0.05,
    )
    return Artifact(
        artifact_id=identifier,
        artifact_type="path_attempt",
        data={"reaction_id": "rxn", "path": path.model_dump()},
    )


def _minimum(identifier: str, path: Path, state: str) -> Artifact:
    return Artifact(
        artifact_id=identifier,
        artifact_type="species_optimized",
        paths={"xyz": str(path)},
        data={
            "species_id": identifier,
            "state": state,
            "reaction_id": "rxn",
            "charge": 0,
            "multiplicity": 1,
            "xyz_path": str(path),
        },
        qc={
            "minimum_accepted": True,
            "is_minimum": True,
            "n_imag": 0,
            "fallback_dummy": False,
        },
    )


def test_rejected_latest_minimum_attempt_is_not_hidden_by_an_old_minimum() -> None:
    old = Artifact(
        artifact_id="old_minimum",
        artifact_type="calculation",
        method={"stage": "dft-minima"},
        data={"species_id": "reactant", "task": "opt_freq"},
        qc={"minimum_accepted": True},
    )
    rejected = Artifact(
        artifact_id="tight_refinement",
        artifact_type="calculation",
        method={"stage": "dft-minima"},
        data={"species_id": "reactant", "task": "opt_freq"},
        qc={"minimum_accepted": False, "n_imag": 1},
    )
    manifest = Manifest.new(run_id="test", stage="dft-minima")
    manifest.extend([old, rejected])

    selected = latest_minimum_attempts_by_species(manifest)

    assert selected["reactant"] is rejected


def test_same_basin_class_has_priority() -> None:
    result = classify_reaction(_reaction(), _basin("same_basin"), require_real_qm=False)

    assert result["classification"] == "same_basin_relaxation"
    assert result["scientific_conclusion_supported"] is True


def test_elementary_step_requires_frequency_ts_and_irc() -> None:
    ts = Artifact(
        artifact_id="ts",
        artifact_type="reaction_validated",
        data={"reaction_id": "rxn"},
        qc={"ts_validated_by_frequency": True, "n_imag": 1},
    )
    irc = Artifact(
        artifact_id="irc",
        artifact_type="reaction_path_validated",
        parents=[ts.artifact_id],
        data={"reaction_id": "rxn"},
        qc={"irc_validated": True},
    )

    result = classify_reaction(
        _reaction(),
        _basin("distinct_basin"),
        ts_artifacts=[ts],
        irc_artifacts=[irc],
        require_real_qm=False,
    )

    assert result["classification"] == "elementary_first_order_saddle"
    assert result["scientific_conclusion_supported"] is True


def test_unrelated_ts_and_irc_cannot_be_combined() -> None:
    ts = Artifact(
        artifact_id="ts",
        artifact_type="reaction_validated",
        data={"reaction_id": "rxn"},
        qc={"ts_validated_by_frequency": True, "n_imag": 1},
    )
    irc = Artifact(
        artifact_id="irc",
        artifact_type="reaction_path_validated",
        parents=["different_ts"],
        data={"reaction_id": "rxn"},
        qc={"irc_validated": True},
    )

    result = classify_reaction(
        _reaction(),
        _basin("distinct_basin"),
        ts_artifacts=[ts],
        irc_artifacts=[irc],
        require_real_qm=False,
    )

    assert result["classification"] == "unresolved"
    assert "irc_not_bound_to_frequency_validated_ts" in result["reasons"]


def test_two_unaudited_paths_do_not_support_barrierless_classification() -> None:
    result = classify_reaction(
        _reaction(),
        _basin("distinct_basin"),
        path_attempts=[
            _attempt("path_5", [-10.0, -10.1, -10.2, -10.3, -10.4]),
            _attempt(
                "path_9",
                [-10.0, -10.05, -10.1, -10.15, -10.2, -10.25, -10.3, -10.35, -10.4],
            ),
        ],
        require_real_qm=False,
    )

    assert result["classification"] == "unresolved"
    assert result["evidence"]["barrierless"]["scope"] == (
        "electronic_energy_path_at_declared_resolution"
    )


def test_one_monotonic_path_stays_unresolved() -> None:
    result = classify_reaction(
        _reaction(),
        _basin("distinct_basin"),
        path_attempts=[_attempt("path_5", [-10.0, -10.1, -10.2, -10.3, -10.4])],
        require_real_qm=False,
    )

    assert result["classification"] == "unresolved"
    assert result["candidate_classification"] == "effectively_barrierless"
    assert result["next_action"] == "repeat_path_at_independent_tighter_resolution"


def test_rejected_endpoint_refinement_precedes_further_path_work() -> None:
    basin = {
        "accepted": False,
        "status": "invalid_endpoints",
        "same_basin": False,
        "distinct_basin": False,
        "reasons": ["reactant_latest_minimum_refinement_rejected"],
    }

    result = classify_reaction(
        _reaction(),
        basin,
        path_attempts=[
            _attempt("path_5", [-10.0, -10.1, -10.2, -10.3, -10.4])
        ],
        require_real_qm=False,
    )

    assert result["classification"] == "unresolved"
    assert result["candidate_classification"] == "effectively_barrierless"
    assert result["next_action"] == (
        "follow_endpoint_imaginary_mode_in_both_directions"
    )


def test_subresolution_path_ripple_is_not_published_as_activation_barrier() -> None:
    path = analyze_reaction_path(
        reaction_id="rxn",
        engine="test_path",
        comments=[
            "energy_hartree=-10.0",
            "energy_hartree=-9.99999",
            "energy_hartree=-10.1",
        ],
        converged=True,
        barrier_threshold_kcal_mol=0.05,
    )

    assert path.classification == "monotonic_no_internal_maximum"
    assert path.barrier_from_reactant_kcal_mol is None
    assert path.maximum_internal_rise_from_reactant_kcal_mol is not None


def test_conflicting_converged_profile_blocks_barrierless_claim() -> None:
    result = classify_reaction(
        _reaction(),
        _basin("distinct_basin"),
        path_attempts=[
            _attempt("path_5", [-10.0, -10.1, -10.2, -10.3, -10.4]),
            _attempt(
                "path_9",
                [-10.0, -10.05, -10.1, -10.15, -10.2, -10.25, -10.3, -10.35, -10.4],
            ),
            _attempt("path_with_peak", [-10.0, -9.9, -10.2]),
        ],
        require_real_qm=False,
    )

    assert result["classification"] == "unresolved"
    assert "converged_path_profiles_disagree" in result["reasons"]


def test_accepted_ensemble_supersedes_its_coarse_seed() -> None:
    coarse = _attempt("coarse_peak", [-10.0, -9.9, -10.2])
    ensemble = Artifact(
        artifact_id="ensemble",
        artifact_type="path_ensemble_assessment",
        data={
            "reaction_id": "rxn",
            "accepted": True,
            "barrierless_at_resolution": True,
            "path_barrier_resolution_kcal_mol": 0.05,
            "profile_energy_tolerance_kcal_mol": 0.05,
            "supersedes_path_attempt_ids": [coarse.artifact_id],
        },
        qc={"path_discretization_converged": True},
    )

    result = classify_reaction(
        _reaction(),
        _basin("distinct_basin"),
        path_attempts=[coarse],
        path_ensemble_assessments=[ensemble],
        require_real_qm=False,
    )

    assert result["classification"] == "effectively_barrierless"
    assert coarse.artifact_id in result["evidence"]["superseded_path_attempt_ids"]


def test_coarse_ensemble_cannot_satisfy_a_tighter_classification_policy() -> None:
    ensemble = Artifact(
        artifact_id="coarse_ensemble",
        artifact_type="path_ensemble_assessment",
        data={
            "reaction_id": "rxn",
            "accepted": True,
            "barrierless_at_resolution": True,
            "path_barrier_resolution_kcal_mol": 0.5,
            "profile_energy_tolerance_kcal_mol": 0.5,
        },
    )

    result = classify_reaction(
        _reaction(),
        _basin("distinct_basin"),
        path_ensemble_assessments=[ensemble],
        require_real_qm=False,
        barrier_resolution_kcal_mol=0.05,
    )

    assert result["classification"] == "unresolved"
    assert (
        "path_ensemble_resolution_coarser_than_classification_policy"
        in result["reasons"]
    )


def test_newer_rejected_ensemble_blocks_an_outdated_barrierless_result() -> None:
    accepted = Artifact(
        artifact_id="accepted_ensemble",
        artifact_type="path_ensemble_assessment",
        data={
            "reaction_id": "rxn",
            "accepted": True,
            "barrierless_at_resolution": True,
        },
    )
    rejected = Artifact(
        artifact_id="rejected_ensemble",
        artifact_type="path_ensemble_assessment",
        data={
            "reaction_id": "rxn",
            "accepted": False,
            "barrierless_at_resolution": False,
            "reasons": ["path_energy_profiles_not_converged"],
        },
        status={"status": "partial", "category": "path_ensemble_unresolved"},
    )

    result = classify_reaction(
        _reaction(),
        _basin("distinct_basin"),
        path_ensemble_assessments=[accepted, rejected],
        require_real_qm=False,
    )

    assert result["classification"] == "unresolved"
    assert result["evidence"]["barrierless"]["path_ensemble_assessment_id"] == (
        rejected.artifact_id
    )


def test_multistep_requires_a_distinct_frequency_minimum(tmp_path: Path) -> None:
    reactant_xyz = write_xyz(
        XYZ(
            ["N", "H", "F"],
            np.array([[0.0, 0.0, 0.0], [1.6, 0.0, 0.0], [2.53, 0.0, 0.0]]),
        ),
        tmp_path / "reactant.xyz",
    )
    product_xyz = write_xyz(
        XYZ(
            ["N", "H", "F"],
            np.array([[0.0, 0.0, 0.0], [1.05, 0.0, 0.0], [2.45, 0.0, 0.0]]),
        ),
        tmp_path / "product.xyz",
    )
    intermediate_xyz = write_xyz(
        XYZ(
            ["N", "H", "F"],
            np.array([[0.0, 0.0, 0.0], [1.2, 1.2, 0.0], [2.8, 0.0, 0.0]]),
        ),
        tmp_path / "intermediate.xyz",
    )
    reactant = _minimum("reactant", reactant_xyz, "reactant")
    product = _minimum("product", product_xyz, "product")
    attempt = _attempt("multi_path", [-10.0, -9.8, -10.2, -9.7, -10.3])
    intermediate = _minimum(
        "intermediate", intermediate_xyz, "reaction_intermediate_candidate"
    )
    intermediate.data.update(
        {"source_path_attempt_id": attempt.artifact_id, "source_image_index": 2}
    )

    result = classify_reaction(
        _reaction(),
        _basin("distinct_basin"),
        path_attempts=[attempt],
        intermediate_minima=[intermediate],
        reactant=reactant,
        product=product,
        require_real_qm=False,
    )

    assert result["classification"] == "multistep_with_intermediate"
    assert result["next_action"] == "split_path_into_elementary_segments"


def test_two_path_peaks_without_minimum_are_only_a_multistep_candidate() -> None:
    attempt = _attempt("multi_path", [-10.0, -9.8, -10.2, -9.7, -10.3])

    result = classify_reaction(
        _reaction(),
        _basin("distinct_basin"),
        path_attempts=[attempt],
        require_real_qm=False,
    )

    assert result["classification"] == "unresolved"
    assert result["candidate_classification"] == "multistep_with_intermediate"
    assert result["next_action"] == "optimize_and_frequency_check_path_well"


def test_production_classification_rejects_unverified_endpoint_evidence() -> None:
    result = classify_reaction(_reaction(), _basin("same_basin"), require_real_qm=True)

    assert result["classification"] == "unresolved"
    assert result["candidate_classification"] == "same_basin_relaxation"


def test_off_coordinate_mode_following_blocks_barrierless_candidate() -> None:
    assessment = Artifact(
        artifact_id="mode_assessment",
        artifact_type="minimum_mode_following_assessment",
        data={
            "reaction_ids": ["rxn_test"],
            "source_species_id": "reactant",
            "outcome": "distinct_basin_descents_off_declared_coordinate",
            "local_basin_conclusion_supported": True,
            "reaction_gate_reasons": [
                "imaginary_mode_does_not_match_declared_reaction_coordinate"
            ],
        },
    )
    reactant = Artifact(
        artifact_id="reactant",
        artifact_type="species_optimized",
        data={"species_id": "reactant"},
    )

    result = classify_reaction(
        _reaction(),
        {"accepted": False, "reasons": ["reactant_refinement_rejected"]},
        minimum_mode_following_assessments=[assessment],
        reactant=reactant,
        require_real_qm=False,
    )

    assert result["classification"] == "unresolved"
    assert result["candidate_classification"] == "unresolved"
    assert result["next_action"] == (
        "revise_declared_reaction_coordinate_or_endpoint_hypothesis"
    )
    assert "stationary_point_mode_does_not_match_declared_reaction_coordinate" in (
        result["reasons"]
    )


def test_mode_projection_is_scoped_to_the_classified_reaction() -> None:
    assessment = Artifact(
        artifact_id="mode_assessment",
        artifact_type="minimum_mode_following_assessment",
        data={
            "reaction_ids": ["rxn", "rxn_other"],
            "source_species_id": "reactant",
            "outcome": "distinct_basin_descents",
            "local_basin_conclusion_supported": True,
            "declared_reaction_mode_projections": [
                {
                    "reaction_id": "rxn",
                    "accepted": False,
                    "reasons": [
                        "imaginary_mode_does_not_match_declared_reaction_coordinate"
                    ],
                },
                {
                    "reaction_id": "rxn_other",
                    "accepted": True,
                    "reasons": [],
                },
            ],
        },
    )
    reactant = Artifact(
        artifact_id="reactant",
        artifact_type="species_optimized",
        data={"species_id": "reactant"},
    )

    result = classify_reaction(
        _reaction(),
        {"accepted": False, "reasons": ["reactant_refinement_rejected"]},
        minimum_mode_following_assessments=[assessment],
        reactant=reactant,
        require_real_qm=False,
    )

    assert result["classification"] == "unresolved"
    assert result["candidate_classification"] == "unresolved"
    assert result["next_action"] == (
        "revise_declared_reaction_coordinate_or_endpoint_hypothesis"
    )
    selected = result["evidence"]["minimum_mode_following"][
        "selected_reaction_mode_projections"
    ]
    assert [item["reaction_id"] for item in selected] == ["rxn"]


def test_path_intermediate_stage_only_extracts_resolved_wells(tmp_path: Path) -> None:
    energies = [-10.0, -9.8, -10.2, -9.7, -10.3]
    images = [
        XYZ(
            ["N", "H", "F"],
            np.array([[0.0, 0.0, 0.0], [1.1 + 0.1 * index, 0.0, 0.0], [2.5, 0.0, 0.0]]),
            f"energy_hartree={energy}",
        )
        for index, energy in enumerate(energies)
    ]
    trajectory = write_xyz_trajectory(images, tmp_path / "path.xyz")
    attempt = _attempt("multi_path", energies)
    attempt.paths["path_xyz"] = str(trajectory)
    manifest = Manifest.new(run_id="test", stage="ts-search")
    manifest.extend([_reaction(), attempt])

    output = PathIntermediatesStage().run(
        manifest,
        {},
        StageContext(
            out_dir=tmp_path / "path-intermediates",
            run_id="test",
            global_config={},
        ),
    )
    seeds = [
        item
        for item in output.latest_artifacts("species")
        if item.data.get("state") == "reaction_intermediate_candidate"
    ]

    assert len(seeds) == 1
    assert seeds[0].data["source_image_index"] == 2
    assert seeds[0].qc["minimum_optimization_required"] is True
    assert seeds[0].qc["scientific_conclusion_supported"] is False


def test_path_intermediate_stage_ignores_superseded_coarse_path(
    tmp_path: Path,
) -> None:
    energies = [-10.0, -9.8, -10.2, -9.7, -10.3]
    images = [
        XYZ(
            ["N", "H", "F"],
            np.array(
                [[0.0, 0.0, 0.0], [1.1 + 0.1 * index, 0.0, 0.0], [2.5, 0.0, 0.0]]
            ),
            f"energy_hartree={energy}",
        )
        for index, energy in enumerate(energies)
    ]
    attempt = _attempt("coarse_path", energies)
    attempt.paths["path_xyz"] = str(
        write_xyz_trajectory(images, tmp_path / "coarse_path.xyz")
    )
    assessment = Artifact(
        artifact_id="ensemble",
        artifact_type="path_ensemble_assessment",
        data={
            "reaction_id": "rxn",
            "accepted": True,
            "supersedes_path_attempt_ids": [attempt.artifact_id],
        },
    )
    manifest = Manifest.new(run_id="test", stage="path-ensemble-assess")
    manifest.extend([_reaction(), attempt, assessment])

    output = PathIntermediatesStage().run(
        manifest,
        {},
        StageContext(
            out_dir=tmp_path / "path-intermediates",
            run_id="test",
            global_config={},
        ),
    )

    assert output.metadata["path_intermediate_seed_count"] == 0

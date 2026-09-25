from pathlib import Path

import numpy as np
import pytest

from hfauto.backends.ts.base import (
    plan_saddle_recovery,
    validate_saddle_seed_hessian,
    validate_ts_frequency_calculation,
)
from hfauto.backends.ts.nwchem_saddle import NWChemSaddleEngine
from hfauto.chemistry.saddle_seeds import (
    aligned_interpolation_xyz,
    bracket_energy_maximum,
    endpoint_biased_fractions,
)
from hfauto.chemistry.xyz import read_xyz
from hfauto.chemistry.xyz_trajectory import endpoint_mode_reference
from hfauto.core.schemas.artifact import Artifact

ROOT = Path(__file__).parents[1]
NEUTRAL = ROOT / "examples/m3_trimethylamine_hf2/neutral.xyz"
SHARED_PROTON = ROOT / "examples/m3_trimethylamine_hf2/shared_proton.xyz"


def test_endpoint_biased_fractions_are_bounded_and_include_endpoints() -> None:
    fractions = endpoint_biased_fractions()

    assert fractions == sorted(fractions)
    assert fractions[0] == 0.0
    assert fractions[-1] == 1.0
    assert fractions[1] < 0.01


def test_nonconverged_saddle_plans_one_hessian_refresh(tmp_path: Path) -> None:
    final_xyz = tmp_path / "final.xyz"
    final_xyz.write_text("1\nfinal\nH 0 0 0\n", encoding="utf-8")
    calculation = Artifact.failure(
        "calc_failed",
        "calculation",
        "optimizer did not converge",
        category="nwchem_failed",
        data={"geometry_converged": False},
    )
    calculation.paths["final_xyz"] = str(final_xyz)

    recovery = plan_saddle_recovery(
        calculation,
        {},
        {
            "driver_trust": 0.15,
            "driver_saddle_step": 0.05,
            "saddle_max_hessian_restarts": 1,
        },
    )

    assert recovery["diagnosis"] == "saddle_optimizer_failed"
    assert recovery["next_action"] == "refresh_hessian_and_restart"
    assert recovery["restart_allowed"] is True
    assert recovery["method_updates"]["driver_trust"] == pytest.approx(0.075)
    assert recovery["method_updates"]["driver_saddle_step"] == pytest.approx(0.02)


def test_saddle_hessian_restart_is_bounded(tmp_path: Path) -> None:
    final_xyz = tmp_path / "final.xyz"
    final_xyz.write_text("1\nfinal\nH 0 0 0\n", encoding="utf-8")
    calculation = Artifact.failure(
        "calc_failed",
        "calculation",
        "optimizer did not converge",
        category="nwchem_failed",
        data={"geometry_converged": False},
    )
    calculation.paths["final_xyz"] = str(final_xyz)

    recovery = plan_saddle_recovery(
        calculation,
        {},
        {"saddle_hessian_restart_count": 1, "saddle_max_hessian_restarts": 1},
    )

    assert recovery["restart_allowed"] is False
    assert recovery["next_action"] == "adaptive_double_ended_path"


def test_first_order_saddle_below_an_endpoint_is_not_validated() -> None:
    species = Artifact(
        artifact_id="ts",
        artifact_type="species",
        paths={"xyz": str(NEUTRAL)},
        data={"species_id": "ts", "xyz_path": str(NEUTRAL)},
    )
    calculation = Artifact(
        artifact_id="calc_ts",
        artifact_type="calculation",
        paths={"final_xyz": str(NEUTRAL)},
        data={
            "electronic_energy_hartree": -100.1,
            "n_imag": 1,
            "imag_freq_cm1": -250.0,
        },
        qc={
            "real_qm_executed": True,
            "fallback_dummy": False,
        },
    )
    method = {
        "basin_assessment": {
            "reactant_evidence": {"electronic_energy_hartree": -100.0},
            "product_evidence": {"electronic_energy_hartree": -101.0},
        }
    }

    assessment = validate_ts_frequency_calculation(
        species,
        calculation,
        method,
        backend="test",
        search_method_evidence_validated=True,
        overlap_evaluator=lambda *_args, **_kwargs: (0.9, "test_mode"),
    )

    assert assessment["ts_energy_above_endpoints"] is False
    assert assessment["validated"] is False
    assert calculation.qc["ts_energy_above_endpoints"] is False


def test_subresolution_barrier_is_separate_from_ts_topology() -> None:
    species = Artifact(
        artifact_id="ts",
        artifact_type="species",
        paths={"xyz": str(NEUTRAL)},
        data={"species_id": "ts", "xyz_path": str(NEUTRAL)},
    )
    calculation = Artifact(
        artifact_id="calc_ts",
        artifact_type="calculation",
        paths={"final_xyz": str(NEUTRAL)},
        data={
            "electronic_energy_hartree": -99.9999985,
            "n_imag": 1,
            "imag_freq_cm1": -250.0,
        },
        qc={
            "real_qm_executed": True,
            "fallback_dummy": False,
        },
    )
    method = {
        "ts_endpoint_energy_tolerance_hartree": 1.0e-5,
        "basin_assessment": {
            "reactant_evidence": {"electronic_energy_hartree": -100.0},
            "product_evidence": {"electronic_energy_hartree": -101.0},
        },
    }

    assessment = validate_ts_frequency_calculation(
        species,
        calculation,
        method,
        backend="test",
        search_method_evidence_validated=True,
        overlap_evaluator=lambda *_args, **_kwargs: (0.9, "test_mode"),
    )

    assert assessment["validated"] is True
    assert assessment["ts_energy_above_endpoints"] is True
    assert assessment["activation_energy_resolved"] is False
    assert calculation.data["electronic_activation_energy_kcal_mol"] == (
        pytest.approx(0.0009412642)
    )
    assert calculation.qc["activation_energy_resolved"] is False


def test_aligned_interpolation_removes_rigid_rotation(
    tmp_path: Path,
) -> None:
    start = tmp_path / "start.xyz"
    target = tmp_path / "target.xyz"
    start.write_text(
        "3\nstart\nO 0 0 0\nH 1 0 0\nH 0 1 0\n", encoding="utf-8"
    )
    target.write_text(
        "3\ntarget\nO 4 2 0\nH 4 3 0\nH 3 2 0\n", encoding="utf-8"
    )

    midpoint = read_xyz(
        aligned_interpolation_xyz(start, target, tmp_path / "mid.xyz", 0.5)
    )

    assert np.allclose(midpoint.coords, read_xyz(start).coords, atol=1.0e-10)


def test_bracket_energy_maximum_requires_an_internal_negative_curvature() -> None:
    profile = bracket_energy_maximum(
        [
            {"fraction": 0.0, "energy_hartree": -10.0, "xyz_path": "r.xyz"},
            {"fraction": 0.1, "energy_hartree": -9.9, "xyz_path": "ts.xyz"},
            {"fraction": 0.4, "energy_hartree": -10.2, "xyz_path": "p.xyz"},
        ]
    )

    assert profile["classification"] == "resolved_internal_maximum"
    assert profile["candidate"]["fraction"] == 0.1
    assert profile["candidate"]["directional_curvature_hartree"] < 0.0

    monotonic = bracket_energy_maximum(
        [
            {"fraction": 0.0, "energy_hartree": -10.0},
            {"fraction": 0.1, "energy_hartree": -10.1},
            {"fraction": 1.0, "energy_hartree": -10.2},
        ]
    )
    assert monotonic["classification"] == "no_bracketed_maximum"
    assert monotonic["candidate"] is None


def test_endpoint_energy_is_never_promoted_to_saddle() -> None:
    profile = bracket_energy_maximum(
        [
            {"fraction": 0.0, "energy_hartree": -9.0},
            {"fraction": 0.5, "energy_hartree": -10.0},
            {"fraction": 1.0, "energy_hartree": -11.0},
        ]
    )

    assert profile["candidate"] is None


def test_multiple_maxima_are_not_collapsed_to_one_saddle() -> None:
    profile = bracket_energy_maximum(
        [
            {"fraction": 0.0, "energy_hartree": -10.0},
            {"fraction": 0.2, "energy_hartree": -9.8},
            {"fraction": 0.4, "energy_hartree": -10.1},
            {"fraction": 0.6, "energy_hartree": -9.7},
            {"fraction": 1.0, "energy_hartree": -10.2},
        ]
    )

    assert profile["classification"] == "multiple_internal_maxima"
    assert len(profile["internal_maxima"]) == 2
    assert profile["candidate"] is None


def test_fraction_grid_rejects_missing_endpoints() -> None:
    with pytest.raises(ValueError, match="include 0 and 1"):
        endpoint_biased_fractions([0.1, 0.5, 0.9])


def test_nwchem_saddle_backend_brackets_before_refinement(
    tmp_path: Path, monkeypatch
) -> None:
    reaction = Artifact(
        artifact_id="rxn_test",
        artifact_type="reaction",
        data={
            "reaction_id": "rxn_test",
            "reaction_coordinate": {
                "terms": [
                    {"kind": "distance", "atoms": [13, 14], "coefficient": 1.0},
                    {"kind": "distance", "atoms": [1, 13], "coefficient": -1.0},
                ]
            },
        },
    )
    reactant = Artifact(
        artifact_id="opt_reactant",
        artifact_type="species_optimized",
        paths={"xyz": str(NEUTRAL)},
        data={
            "species_id": "reactant",
            "xyz_path": str(NEUTRAL),
            "charge": 0,
            "multiplicity": 1,
        },
    )
    product = Artifact(
        artifact_id="opt_product",
        artifact_type="species_optimized",
        paths={"xyz": str(SHARED_PROTON)},
        data={
            "species_id": "product",
            "xyz_path": str(SHARED_PROTON),
            "charge": 0,
            "multiplicity": 1,
        },
    )

    def fake_single_point(self, species, method, workdir):
        fraction = float(species.data["interpolation_fraction"])
        energy = -10.0 + 0.003 - 2.0 * (fraction - 0.02) ** 2
        return Artifact(
            artifact_id=f"calc_{fraction}",
            artifact_type="calculation",
            parents=[species.artifact_id],
            data={
                "electronic_energy_hartree": energy,
                "program_version": "7.2.3",
            },
            qc={
                "scf_converged": True,
                "normal_termination": True,
                "real_qm_executed": True,
                "fallback_dummy": False,
                "dispersion_applied": True,
            },
        )

    monkeypatch.setattr(
        "hfauto.backends.ts.nwchem_saddle.NWChemEngine.single_point",
        fake_single_point,
    )
    result = NWChemSaddleEngine().search_ts(
        reaction,
        reactant,
        product,
        {
            "path_strategy": "bracketed_saddle_search",
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
            "required_program_version": "7.2.3",
            "saddle_bracket_fractions": [0.0, 0.01, 0.02, 0.04, 1.0],
            "saddle_bracket_only": True,
            "basin_assessment": {
                "reactant_evidence": {"electronic_energy_hartree": -10.0},
                "product_evidence": {"electronic_energy_hartree": -10.2},
            },
        },
        tmp_path / "saddle",
    )

    scan = next(
        artifact
        for artifact in result.artifacts
        if artifact.artifact_type == "saddle_seed_scan"
    )
    attempt = next(
        artifact
        for artifact in result.artifacts
        if artifact.artifact_type == "saddle_attempt"
    )
    assert result.success is False
    assert scan.data["classification"] == "resolved_internal_maximum"
    assert scan.data["candidate"]["fraction"] == 0.02
    assert attempt.data["diagnosis"] == "resolved_saddle_candidate"


def test_default_bracket_does_not_promote_subresolution_ripple(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reaction = Artifact(
        artifact_id="rxn_test",
        artifact_type="reaction",
        data={"reaction_id": "rxn_test"},
    )
    endpoints = [
        Artifact(
            artifact_id=name,
            artifact_type="species_optimized",
            paths={"xyz": str(path)},
            data={
                "species_id": name,
                "xyz_path": str(path),
                "charge": 0,
                "multiplicity": 1,
            },
        )
        for name, path in (("reactant", NEUTRAL), ("product", SHARED_PROTON))
    ]

    def fake_single_point(self, species, method, workdir):
        fraction = float(species.data["interpolation_fraction"])
        energy = -10.0 + 5.0e-6 - 2.0e-5 * (fraction - 0.5) ** 2
        return Artifact(
            artifact_id=f"calc_{fraction}",
            artifact_type="calculation",
            data={
                "electronic_energy_hartree": energy,
                "program_version": "7.2.3",
            },
            qc={
                "scf_converged": True,
                "normal_termination": True,
                "real_qm_executed": True,
                "fallback_dummy": False,
                "dispersion_applied": True,
            },
        )

    monkeypatch.setattr(
        "hfauto.backends.ts.nwchem_saddle.NWChemEngine.single_point",
        fake_single_point,
    )
    result = NWChemSaddleEngine().search_ts(
        reaction,
        endpoints[0],
        endpoints[1],
        {
            "path_strategy": "bracketed_saddle_search",
            "disp_vdw": 3,
            "required_program_version": "7.2.3",
            "saddle_bracket_fractions": [0.0, 0.25, 0.5, 0.75, 1.0],
            "saddle_bracket_only": True,
            "basin_assessment": {
                "reactant_evidence": {"electronic_energy_hartree": -10.0},
                "product_evidence": {"electronic_energy_hartree": -10.00001},
            },
        },
        tmp_path / "subresolution",
    )

    scan = next(
        artifact
        for artifact in result.artifacts
        if artifact.artifact_type == "saddle_seed_scan"
    )
    assert scan.data["classification"] == "no_bracketed_maximum"
    assert scan.data["candidate"] is None


def test_nwchem_segment_bracket_uses_local_path_tangent(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    segment_start = tmp_path / "segment_start.xyz"
    segment_end = tmp_path / "segment_end.xyz"
    segment_start.write_text(
        "3\nsegment start\nO 0 0 0\nH 1 0 0\nH 0 1 0\n",
        encoding="utf-8",
    )
    segment_end.write_text(
        "3\nsegment end\nO 0 0 0\nH 1.1 0 0\nH 0 1 0\n",
        encoding="utf-8",
    )
    reaction = Artifact(
        artifact_id="rxn_segment",
        artifact_type="reaction",
        data={"reaction_id": "rxn_segment"},
    )
    endpoints = [
        Artifact(
            artifact_id=name,
            artifact_type="species_optimized",
            paths={"xyz": str(path)},
            data={
                "species_id": name,
                "xyz_path": str(path),
                "charge": 0,
                "multiplicity": 1,
            },
        )
        for name, path in (("reactant", NEUTRAL), ("product", SHARED_PROTON))
    ]

    def fake_single_point(self, species, method, workdir):
        fraction = float(species.data["interpolation_fraction"])
        energy = -10.0 + 0.001 - (fraction - 0.5) ** 2
        return Artifact(
            artifact_id=f"calc_{fraction}",
            artifact_type="calculation",
            data={
                "electronic_energy_hartree": energy,
                "program_version": "7.2.3",
            },
            qc={
                "scf_converged": True,
                "normal_termination": True,
                "real_qm_executed": True,
                "fallback_dummy": False,
                "dispersion_applied": True,
            },
        )

    monkeypatch.setattr(
        "hfauto.backends.ts.nwchem_saddle.NWChemEngine.single_point",
        fake_single_point,
    )
    result = NWChemSaddleEngine().search_ts(
        reaction,
        endpoints[0],
        endpoints[1],
        {
            "path_strategy": "bracketed_saddle_search",
            "disp_vdw": 3,
            "required_program_version": "7.2.3",
            "saddle_bracket_fractions": [0.0, 0.25, 0.5, 0.75, 1.0],
            "saddle_bracket_only": True,
            "saddle_bracket_segment": {
                "accepted": True,
                "start_index": 3,
                "end_index": 4,
                "start_xyz": str(segment_start),
                "target_xyz": str(segment_end),
                "start_energy_hartree": -10.25,
                "target_energy_hartree": -10.25,
            },
        },
        tmp_path / "segment_scan",
    )

    scan = next(
        artifact
        for artifact in result.artifacts
        if artifact.artifact_type == "saddle_seed_scan"
    )
    attempt = next(
        artifact
        for artifact in result.artifacts
        if artifact.artifact_type == "saddle_attempt"
    )
    candidate = attempt.data["candidate"]
    reference = candidate["reaction_mode_reference"]
    displacements = np.asarray(reference["cartesian_displacements"])
    expected = np.asarray(
        endpoint_mode_reference(
            read_xyz(segment_start), read_xyz(segment_end)
        )["cartesian_displacements"]
    )
    assert reference["source"] == "validated_path_segment:3-4"
    assert scan.data["points"][0]["source"] == "validated_reaction_path_image"
    assert np.allclose(displacements, expected)


def test_standalone_refinement_rejects_unvalidated_seed(tmp_path: Path) -> None:
    reaction = Artifact(
        artifact_id="rxn_test",
        artifact_type="reaction",
        data={"reaction_id": "rxn_test"},
    )
    reactant = Artifact(
        artifact_id="reactant",
        artifact_type="species_optimized",
        paths={"xyz": str(NEUTRAL)},
        data={
            "species_id": "reactant",
            "xyz_path": str(NEUTRAL),
            "charge": 0,
            "multiplicity": 1,
        },
    )
    product = Artifact(
        artifact_id="product",
        artifact_type="species_optimized",
        paths={"xyz": str(SHARED_PROTON)},
        data={
            "species_id": "product",
            "xyz_path": str(SHARED_PROTON),
            "charge": 0,
            "multiplicity": 1,
        },
    )

    result = NWChemSaddleEngine().search_ts(
        reaction,
        reactant,
        product,
        {
            "path_strategy": "saddle_refinement",
            "saddle_seed_candidate": {"xyz_path": str(NEUTRAL)},
            "saddle_seed_evidence_validated": False,
        },
        tmp_path,
    )

    failure = result.artifacts[-1]
    assert result.success is False
    assert failure.status.category == "saddle_refinement_seed_invalid"


def test_seed_hessian_gate_does_not_promote_nonstationary_seed() -> None:
    seed = Artifact(
        artifact_id="seed",
        artifact_type="species",
        paths={"xyz": str(NEUTRAL)},
        data={"species_id": "seed", "xyz_path": str(NEUTRAL)},
    )
    calculation = Artifact(
        artifact_id="seed_hessian",
        artifact_type="calculation",
        paths={"final_xyz": str(NEUTRAL)},
        data={
            "n_imag": 1,
            "imag_freq_cm1": -500.0,
            "frequency_count_complete": True,
            "program_version": "7.2.3",
        },
        qc={
            "real_qm_executed": True,
            "fallback_dummy": False,
            "dispersion_applied": True,
        },
    )

    assessment = validate_saddle_seed_hessian(
        seed,
        calculation,
        {"disp_vdw": 3, "required_program_version": "7.2.3"},
        backend="nwchem_saddle",
        overlap_evaluator=lambda *args, **kwargs: (0.9, "test_projection"),
    )

    assert assessment["accepted"] is True
    assert assessment["stationary_point_validated"] is False
    assert calculation.qc["saddle_seed_hessian_accepted"] is True
    assert calculation.qc["ts_validated_by_frequency"] is False


def test_seed_hessian_selects_target_from_multiple_negative_modes() -> None:
    structure = read_xyz(NEUTRAL)
    reference = np.zeros_like(structure.coords)
    reference[13, 2] = 1.0
    transverse = np.zeros_like(structure.coords)
    transverse[15, 0] = 1.0
    seed = Artifact(
        artifact_id="seed",
        artifact_type="species",
        paths={"xyz": str(NEUTRAL)},
        data={
            "species_id": "seed",
            "xyz_path": str(NEUTRAL),
            "reaction_mode_reference": {
                "kind": "local_path_tangent",
                "component_units": "angstrom",
                "cartesian_displacements": reference.tolist(),
            },
        },
    )
    calculation = Artifact(
        artifact_id="seed_hessian",
        artifact_type="calculation",
        paths={"final_xyz": str(NEUTRAL)},
        data={
            "n_imag": 2,
            "imag_freq_cm1": -500.0,
            "frequency_count_complete": True,
            "program_version": "7.2.3",
            "projected_imaginary_modes": [
                {
                    "mode_number": 1,
                    "frequency_cm1": -500.0,
                    "component_units": "angstrom",
                    "cartesian_displacements": transverse.tolist(),
                },
                {
                    "mode_number": 2,
                    "frequency_cm1": -200.0,
                    "component_units": "angstrom",
                    "cartesian_displacements": reference.tolist(),
                },
            ],
        },
        qc={
            "real_qm_executed": True,
            "fallback_dummy": False,
            "dispersion_applied": True,
        },
    )

    assessment = validate_saddle_seed_hessian(
        seed,
        calculation,
        {"disp_vdw": 3, "required_program_version": "7.2.3"},
        backend="nwchem_saddle",
    )

    assert assessment["accepted"] is True
    assert assessment["spectral_qc"]["selected_mode"]["mode_number"] == 2
    assert calculation.qc["seed_has_one_target_imaginary_mode"] is False
    assert calculation.qc["seed_has_target_imaginary_mode"] is True


def test_nwchem_refinement_rejects_ambiguous_multiple_mode_mapping(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    structure = read_xyz(NEUTRAL)
    reference = np.zeros_like(structure.coords)
    reference[13, 2] = 1.0
    transverse = np.zeros_like(structure.coords)
    transverse[15, 0] = 1.0
    reaction = Artifact(
        artifact_id="rxn_test",
        artifact_type="reaction",
        data={"reaction_id": "rxn_test"},
    )
    endpoints = [
        Artifact(
            artifact_id=name,
            artifact_type="species_optimized",
            paths={"xyz": str(path)},
            data={
                "species_id": name,
                "xyz_path": str(path),
                "charge": 0,
                "multiplicity": 1,
            },
        )
        for name, path in (("reactant", NEUTRAL), ("product", SHARED_PROTON))
    ]
    seed_hessian = Artifact(
        artifact_id="seed_hessian",
        artifact_type="calculation",
        paths={"final_xyz": str(NEUTRAL)},
        data={
            "n_imag": 2,
            "frequency_count_complete": True,
            "program_version": "7.2.3",
            "projected_imaginary_modes": [
                {
                    "mode_number": 1,
                    "frequency_cm1": -500.0,
                    "component_units": "angstrom",
                    "cartesian_displacements": transverse.tolist(),
                },
                {
                    "mode_number": 2,
                    "frequency_cm1": -200.0,
                    "component_units": "angstrom",
                    "cartesian_displacements": reference.tolist(),
                },
            ],
        },
        qc={
            "real_qm_executed": True,
            "fallback_dummy": False,
            "dispersion_applied": True,
        },
    )

    monkeypatch.setattr(
        "hfauto.backends.ts.nwchem_saddle.NWChemEngine.saddle_frequency",
        lambda *_args, **_kwargs: pytest.fail("ambiguous MODDIR must not run"),
    )
    result = NWChemSaddleEngine().search_ts(
        reaction,
        endpoints[0],
        endpoints[1],
        {
            "path_strategy": "saddle_refinement",
            "saddle_seed_candidate": {
                "xyz_path": str(NEUTRAL),
                "reaction_mode_reference": {
                    "kind": "local_path_tangent",
                    "component_units": "angstrom",
                    "cartesian_displacements": reference.tolist(),
                },
            },
            "saddle_seed_evidence_validated": True,
            "saddle_seed_hessian_evidence": seed_hessian.model_dump(),
            "disp_vdw": 3,
            "required_program_version": "7.2.3",
        },
        tmp_path,
    )

    attempt = next(
        artifact
        for artifact in result.artifacts
        if artifact.artifact_type == "saddle_attempt"
    )
    assert result.success is False
    assert result.artifacts[-1].status.category == (
        "saddle_seed_mode_mapping_ambiguous"
    )
    assert attempt.qc["seed_driver_mode_mapping_validated"] is False


def test_refinement_reuses_seed_hessian_and_stops_without_target_mode(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reaction = Artifact(
        artifact_id="rxn_test",
        artifact_type="reaction",
        data={"reaction_id": "rxn_test"},
    )
    reactant = Artifact(
        artifact_id="reactant",
        artifact_type="species_optimized",
        paths={"xyz": str(NEUTRAL)},
        data={
            "species_id": "reactant",
            "xyz_path": str(NEUTRAL),
            "charge": 0,
            "multiplicity": 1,
        },
    )
    product = Artifact(
        artifact_id="product",
        artifact_type="species_optimized",
        paths={"xyz": str(SHARED_PROTON)},
        data={
            "species_id": "product",
            "xyz_path": str(SHARED_PROTON),
            "charge": 0,
            "multiplicity": 1,
        },
    )
    saddle_called = False

    seed_hessian = Artifact(
        artifact_id="seed_hessian",
        artifact_type="calculation",
        paths={"final_xyz": str(NEUTRAL)},
        data={
            "n_imag": 0,
            "imag_freq_cm1": None,
            "frequency_count_complete": True,
            "program_version": "7.2.3",
        },
        qc={
            "real_qm_executed": True,
            "fallback_dummy": False,
            "dispersion_applied": True,
        },
    )

    def fake_frequency(self, species, method, workdir):
        raise AssertionError("a supplied seed Hessian must not be recomputed")

    def fake_saddle_frequency(self, species, method, workdir):
        nonlocal saddle_called
        saddle_called = True
        raise AssertionError("rejected seed must not start saddle optimization")

    monkeypatch.setattr(
        "hfauto.backends.ts.nwchem_saddle.NWChemEngine.frequency",
        fake_frequency,
    )
    monkeypatch.setattr(
        "hfauto.backends.ts.nwchem_saddle.NWChemEngine.saddle_frequency",
        fake_saddle_frequency,
    )

    result = NWChemSaddleEngine().search_ts(
        reaction,
        reactant,
        product,
        {
            "path_strategy": "saddle_refinement",
            "saddle_seed_candidate": {"xyz_path": str(NEUTRAL)},
            "saddle_seed_evidence_validated": True,
            "saddle_seed_hessian_evidence": seed_hessian.model_dump(),
            "disp_vdw": 3,
            "required_program_version": "7.2.3",
        },
        tmp_path,
    )

    attempt = next(
        artifact
        for artifact in result.artifacts
        if artifact.artifact_type == "saddle_attempt"
    )
    assert result.success is False
    assert saddle_called is False
    assert result.artifacts[-1].status.category == "saddle_seed_hessian_rejected"
    assert attempt.data["diagnosis"] == "saddle_seed_hessian_rejected"
    assert attempt.qc["seed_hessian_reused"] is True

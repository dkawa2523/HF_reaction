from pathlib import Path

import numpy as np
import pytest
import yaml

from hfauto.backends.ts import pysisyphus_saddle as saddle_module
from hfauto.backends.ts.pysisyphus import (
    PysisyphusEngine,
    _dimer_convergence_evidence,
    _endpoint_optimization_evidence,
)
from hfauto.backends.ts.pysisyphus_saddle import PysisyphusSaddleEngine
from hfauto.chemistry.basin_identity import (
    assess_basin_pair,
    assess_minimum_pair_basin_identity,
    match_geometry_to_basin,
)
from hfauto.chemistry.path_diagnostics import (
    assess_optimizer_stagnation,
    diagnose_path_attempt,
    make_path_attempt_record,
    parse_neb_optimization_history,
)
from hfauto.chemistry.reaction_path_qc import (
    make_displaced_reaction_coordinate_seed,
)
from hfauto.chemistry.reaction_profile import analyze_reaction_path
from hfauto.chemistry.xyz import XYZ, read_xyz, write_xyz
from hfauto.core.executables import CommandResult
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import StageContext
from hfauto.stages.irc import IRCStage
from hfauto.stages.ts_search import TSSearchStage
from hfauto.workflow.reaction_state import decide_reaction_case

ROOT = Path(__file__).parents[1]
NEUTRAL = ROOT / "examples/m3_trimethylamine_hf2/neutral.xyz"
SHARED_PROTON = ROOT / "examples/m3_trimethylamine_hf2/shared_proton.xyz"


def _species(artifact_id: str, xyz_path: Path) -> Artifact:
    return Artifact(
        artifact_id=artifact_id,
        artifact_type="species_optimized",
        paths={"xyz": str(xyz_path)},
        data={
            "species_id": artifact_id.removeprefix("opt_"),
            "xyz_path": str(xyz_path),
            "charge": 0,
            "multiplicity": 1,
            "n_imag": 0,
        },
        qc={
            "minimum_accepted": True,
            "is_minimum": True,
            "scf_converged": True,
            "geometry_converged": True,
            "geometry_sane": True,
            "fallback_dummy": False,
        },
    )


def _calculation(species_id: str, energy: float) -> Artifact:
    return Artifact(
        artifact_id=f"calc_{species_id}",
        artifact_type="calculation",
        data={
            "species_id": species_id,
            "n_imag": 0,
            "frequency_count_complete": True,
            "electronic_energy_hartree": energy,
        },
        qc={"is_minimum": True},
    )


def _reaction() -> Artifact:
    return Artifact(
        artifact_id="rxn_tma_hf2_proton_reorganization",
        artifact_type="reaction",
        data={
            "reaction_id": "rxn_tma_hf2_proton_reorganization",
            "reactant_species_id": "tma_hf2_neutral",
            "product_species_id": "tma_hf2_shared_proton",
            "mechanism_family": "hf_cluster_proton_reorganization",
            "reaction_coordinate": {
                "min_change": 0.05,
                "terms": [
                    {"kind": "distance", "atoms": [13, 14], "coefficient": 1.0},
                    {"kind": "distance", "atoms": [1, 13], "coefficient": -1.0},
                ],
            },
        },
    )


def test_tma_hf2_frequency_minima_are_distinct_basins() -> None:
    reactant = _species("opt_tma_hf2_neutral", NEUTRAL)
    product = _species("opt_tma_hf2_shared_proton", SHARED_PROTON)

    assessment = assess_basin_pair(
        _reaction(),
        reactant,
        product,
        _calculation("tma_hf2_neutral", -374.75344230428698),
        _calculation("tma_hf2_shared_proton", -374.75493512258799),
    )

    assert assessment["status"] == "distinct_basin"
    assert assessment["reactant_evidence"]["frequency_validated"] is True
    assert assessment["product_evidence"]["frequency_validated"] is True
    assert assessment["endpoint_chemistry_qc"]["accepted"] is True


def test_basin_identity_requires_explicit_complete_frequency_evidence() -> None:
    reactant = _species("opt_tma_hf2_neutral", NEUTRAL)
    product = _species("opt_tma_hf2_shared_proton", SHARED_PROTON)
    incomplete = _calculation("tma_hf2_neutral", -374.75344230428698)
    incomplete.data.pop("frequency_count_complete")

    assessment = assess_basin_pair(
        _reaction(),
        reactant,
        product,
        incomplete,
        _calculation("tma_hf2_shared_proton", -374.75493512258799),
    )

    assert assessment["accepted"] is False
    assert assessment["status"] == "invalid_endpoints"
    assert assessment["reactant_evidence"]["frequency_validated"] is False
    assert "frequency_validated_minimum_pair_missing" in assessment["reasons"]


def test_atom_reordering_cannot_create_a_distinct_basin(tmp_path: Path) -> None:
    first_xyz = read_xyz(NEUTRAL)
    permutation = list(reversed(range(len(first_xyz.symbols))))
    reordered_path = write_xyz(
        XYZ(
            symbols=[first_xyz.symbols[index] for index in permutation],
            coords=first_xyz.coords[permutation],
            comment="same geometry with a different atom order",
        ),
        tmp_path / "reordered.xyz",
    )
    first = _species("opt_first", NEUTRAL)
    second = _species("opt_second", Path(reordered_path))

    assessment = assess_minimum_pair_basin_identity(
        first,
        second,
        _calculation("first", -374.75344230428698),
        _calculation("second", -374.75344230428698),
    )

    assert assessment["same_basin"] is True
    assert assessment["distinct_basin"] is False
    assert assessment["distance_spectrum_rmsd_A"] == pytest.approx(0.0)


def test_same_minimum_is_not_promoted_to_an_elementary_step() -> None:
    reactant = _species("opt_tma_hf2_neutral", NEUTRAL)
    duplicate = _species("opt_duplicate", NEUTRAL)
    assessment = assess_basin_pair(
        _reaction(),
        reactant,
        duplicate,
        _calculation("tma_hf2_neutral", -374.75344230428698),
        _calculation("duplicate", -374.75344230428698),
    )

    case = decide_reaction_case(
        reaction_id="rxn",
        basin_assessment=assessment,
    )

    assert assessment["status"] == "same_basin"
    assert case["state"] == "same_basin"
    assert case["workflow_complete"] is True
    assert case["ts_search_allowed"] is False


def test_real_tma_neb_signature_is_stagnant_and_endpoint_inconsistent() -> None:
    history = """
 # NEB Path iteration = 64
 # Gmax = 4.6252009397410900E-004
 # Grms = 9.6494026269543833E-005
 # Xmax = 9.2266285306807916E-004
 # Xrms = 3.0729060543820065E-004
 # NEB Path iteration = 65
 # Gmax = 4.8600314520801314E-004
 # Grms = 9.4660546580671164E-005
 # Xmax = 8.3261668354017715E-004
 # Xrms = 3.1964615845197337E-004
 # NEB Path iteration = 66
 # Gmax = 4.8608222920231987E-004
 # Grms = 9.4720025094706791E-005
 # Xmax = 2.0019541580040823E-011
 # Xrms = 2.2483216679388994E-010
"""
    optimization = parse_neb_optimization_history(history, converged=False)
    path = analyze_reaction_path(
        reaction_id="rxn_tma_hf2_proton_reorganization",
        engine="nwchem_neb",
        comments=[
            "energy_hartree=-374.75344230428698",
            "energy_hartree=-374.75353285203198",
            "energy_hartree=-374.75389507311286",
            "energy_hartree=-374.75449350782264",
            "energy_hartree=-374.75493512258799",
        ],
        converged=False,
        barrier_threshold_kcal_mol=0.05,
    )
    decision = diagnose_path_attempt(path, "distinct_basin", optimization)

    assert optimization.iterations == 66
    assert optimization.stagnant is True
    assert path.classification == "monotonic_no_internal_maximum"
    assert decision["diagnosis"] == "endpoint_path_inconsistency"
    assert decision["next_action"] == "bracket_narrow_saddle"


def test_optimizer_stagnation_requires_distance_and_lack_of_improvement() -> None:
    stalled = assess_optimizer_stagnation(
        [0.0045] * 20 + [0.0042] * 20,
        convergence_target=0.00002,
    )
    improving = assess_optimizer_stagnation(
        [0.0045] * 20 + [0.0020] * 20,
        convergence_target=0.00002,
    )
    near_target = assess_optimizer_stagnation(
        [0.0045] * 20 + [0.0005] * 19 + [0.0009],
        convergence_target=0.00002,
    )

    assert stalled["stagnant"] is True
    assert improving["stagnant"] is False
    assert near_target["stagnant"] is False


def test_monotonic_path_replans_to_bracketed_saddle_search() -> None:
    path = analyze_reaction_path(
        reaction_id="rxn",
        engine="nwchem_neb",
        comments=["energy=-10.0", "energy=-10.1", "energy=-10.2"],
        converged=False,
    )
    attempt = make_path_attempt_record(
        attempt_id="attempt_1",
        reaction_id="rxn",
        strategy="double_ended_path",
        engine="nwchem_neb",
        endpoint_basin_status="distinct_basin",
        path=path,
    )
    case = decide_reaction_case(
        reaction_id="rxn",
        basin_assessment={"status": "distinct_basin", "accepted": True},
        path_attempts=[attempt],
    )

    assert case["state"] == "path_replan_ready"
    assert case["next_strategy"] == "bracketed_saddle_search"
    assert case["authorized_path_engines"] == ["nwchem_saddle"]
    assert case["ts_search_allowed"] is True
    assert "double_ended_path" in case["attempted_strategies"]

    failed_bracket = {
        "reaction_id": "rxn",
        "strategy": "bracketed_saddle_search",
        "diagnosis": "saddle_seed_bracket_failed",
    }
    blocked = decide_reaction_case(
        reaction_id="rxn",
        basin_assessment={"status": "distinct_basin", "accepted": True},
        path_attempts=[attempt, failed_bracket],
    )
    assert blocked["state"] == "path_replan_ready"
    assert blocked["next_strategy"] == "adaptive_double_ended_path"
    assert blocked["authorized_path_engines"] == ["nwchem_string"]
    assert blocked["path_attempt_count"] == 2
    assert blocked["ts_search_allowed"] is True


def test_unconverged_internal_maximum_is_only_a_saddle_seed() -> None:
    path = analyze_reaction_path(
        reaction_id="rxn",
        engine="nwchem_neb",
        comments=["energy=-10.0", "energy=-9.9", "energy=-10.1"],
        converged=False,
        barrier_threshold_kcal_mol=0.05,
    )

    decision = diagnose_path_attempt(path, "distinct_basin")

    assert decision["diagnosis"] == "resolved_saddle_candidate"
    assert decision["next_action"] == "refine_saddle"
    assert "path_not_converged" in decision["reasons"]
    assert "unconverged_path_energy_is_seed_only_not_a_barrier" in decision["reasons"]


def test_artifact_saddle_candidate_routes_to_standalone_refinement() -> None:
    candidate = {
        "xyz_path": "candidate.xyz",
        "directional_curvature_hartree": -0.1,
    }
    attempt = Artifact(
        artifact_id="saddle_attempt",
        artifact_type="saddle_attempt",
        data={
            "reaction_id": "rxn",
            "strategy": "bracketed_saddle_search",
            "diagnosis": "resolved_saddle_candidate",
            "candidate": candidate,
            "all_scan_points_valid": True,
        },
    )

    case = decide_reaction_case(
        reaction_id="rxn",
        basin_assessment={"status": "distinct_basin", "accepted": True},
        path_attempts=[attempt],
    )

    assert case["state"] == "saddle_refinement_ready"
    assert case["next_strategy"] == "saddle_refinement"
    assert case["authorized_path_engines"] == [
        "pysisyphus_saddle",
        "nwchem_saddle",
    ]
    assert case["saddle_seed_candidate"] == candidate
    assert case["saddle_seed_evidence_validated"] is True


def test_resolved_profile_without_geometry_seed_replans_instead_of_looping() -> None:
    attempt = Artifact(
        artifact_id="incomplete_saddle_attempt",
        artifact_type="saddle_attempt",
        data={
            "reaction_id": "rxn",
            "strategy": "saddle_refinement",
            "diagnosis": "resolved_saddle_candidate",
            "candidate": {},
            "all_scan_points_valid": False,
        },
    )

    case = decide_reaction_case(
        reaction_id="rxn",
        basin_assessment={"status": "distinct_basin", "accepted": True},
        path_attempts=[attempt],
    )

    assert case["state"] == "path_replan_ready"
    assert case["next_strategy"] == "adaptive_double_ended_path"
    assert case["authorized_path_engines"] == ["nwchem_string"]
    assert "resolved_maximum_missing_validated_geometry_seed" in case["reasons"]


def test_rejected_seed_routes_to_string_with_seed_lineage() -> None:
    candidate = {"xyz_path": "rejected_seed.xyz"}
    attempt = Artifact(
        artifact_id="saddle_attempt",
        artifact_type="saddle_attempt",
        data={
            "reaction_id": "rxn",
            "strategy": "bracketed_saddle_search",
            "diagnosis": "saddle_seed_hessian_rejected",
            "candidate": candidate,
        },
    )

    case = decide_reaction_case(
        reaction_id="rxn",
        basin_assessment={"status": "distinct_basin", "accepted": True},
        path_attempts=[attempt],
    )

    assert case["state"] == "path_replan_ready"
    assert case["next_strategy"] == "adaptive_double_ended_path"
    assert case["authorized_path_engines"] == ["nwchem_string"]
    assert case["adaptive_path_seed_candidate"] == candidate
    assert case["ts_search_allowed"] is True


def test_nonconverged_saddle_routes_to_one_hessian_refreshed_restart() -> None:
    attempt = {
        "reaction_id": "rxn",
        "strategy": "saddle_refinement",
        "diagnosis": "saddle_optimizer_failed",
        "saddle_recovery": {
            "next_action": "refresh_hessian_and_restart",
            "restart_allowed": True,
            "restart_seed_xyz": "final-080.xyz",
            "method_updates": {
                "saddle_hessian_restart_count": 1,
                "driver_trust": 0.075,
                "driver_saddle_step": 0.02,
            },
        },
    }

    case = decide_reaction_case(
        reaction_id="rxn",
        basin_assessment={"status": "distinct_basin", "accepted": True},
        path_attempts=[attempt],
    )

    assert case["state"] == "saddle_hessian_refresh_ready"
    assert case["next_strategy"] == "saddle_refinement"
    assert case["saddle_seed_candidate"] == {
        "xyz_path": "final-080.xyz",
        "source": "failed_saddle_final_geometry",
    }
    assert case["method_updates"]["driver_trust"] == 0.075


def test_converged_saddle_without_final_frequency_routes_to_frequency_only() -> None:
    attempt = {
        "reaction_id": "rxn",
        "strategy": "saddle_refinement",
        "diagnosis": "final_frequency_missing",
        "saddle_recovery": {"restart_seed_xyz": "stationary.xyz"},
    }

    case = decide_reaction_case(
        reaction_id="rxn",
        basin_assessment={"status": "distinct_basin", "accepted": True},
        path_attempts=[attempt],
    )

    assert case["state"] == "saddle_frequency_completion_ready"
    assert case["next_action"] == "run_fixed_geometry_frequency"
    assert case["fixed_geometry_xyz"] == "stationary.xyz"
    assert case["ts_search_allowed"] is False


def test_endpoint_optimization_requires_both_optimizer_conclusions(
    tmp_path: Path,
) -> None:
    for direction in ("forward", "backward"):
        (tmp_path / f"{direction}_end_opt.xyz").write_text(
            "1\nend\nH 0 0 0\n", encoding="utf-8"
        )
        (tmp_path / f"{direction}_end_optimizer.log").write_text(
            "Converged!\n", encoding="utf-8"
        )
    assert _endpoint_optimization_evidence(tmp_path)["accepted"] is True

    (tmp_path / "backward_end_optimizer.log").write_text(
        "Converged!\nNumber of cycles exceeded!\n", encoding="utf-8"
    )
    evidence = _endpoint_optimization_evidence(tmp_path)
    assert evidence["accepted"] is False
    assert evidence["branches"]["backward"]["converged"] is False


def test_endpoint_dimer_input_uses_declared_coordinate_and_real_qm(
    tmp_path: Path,
) -> None:
    endpoint = _species("opt_tma_hf2_neutral", NEUTRAL)
    input_path = PysisyphusEngine().render_dimer_input(
        _reaction(),
        endpoint,
        {
            "program": "nwchem",
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
            "required_program_version": "7.2.3",
            "ncores": 4,
        },
        tmp_path / "pysis_dimer.yaml",
    )
    config = yaml.safe_load(input_path.read_text(encoding="utf-8"))
    orientation = Path(config["calc"]["N_raw"])

    assert config["geom"]["type"] == "cart"
    assert config["calc"]["type"] == "dimer"
    assert config["calc"]["calc"]["type"] == "qcengine"
    assert config["calc"]["calc"]["program"] == "nwchem"
    assert config["calc"]["calc"]["keywords"]["dft__disp"] == "vdw 3"
    assert config["calc"]["calc"]["keywords"]["basis__spherical"] is True
    assert config["tsopt"]["type"] == "plbfgs"
    assert len(orientation.read_text(encoding="utf-8").splitlines()) == 51


def test_path_seed_rsirfo_uses_exact_hessian_and_bounded_trust(
    tmp_path: Path,
) -> None:
    seed = _species("string_saddle_seed", NEUTRAL)
    input_path = PysisyphusSaddleEngine().render_input(
        seed,
        {
            "program": "nwchem",
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
            "required_program_version": "7.2.3",
            "ncores": 4,
            "ts_coordinate_type": "tric",
            "ts_trust_radius": 0.08,
            "ts_trust_min": 0.004,
            "ts_trust_max": 0.16,
        },
        tmp_path / "pysis_saddle.yaml",
    )
    config = yaml.safe_load(input_path.read_text(encoding="utf-8"))

    assert config["geom"]["type"] == "tric"
    assert config["calc"]["type"] == "qcengine"
    assert config["calc"]["keywords"]["dft__disp"] == "vdw 3"
    assert config["calc"]["keywords"]["basis__spherical"] is True
    assert config["tsopt"] == {
        "type": "rsirfo",
        "hessian_init": "calc",
        "hessian_update": "bofill",
        "root": 0,
        "assert_neg_eigval": True,
        "trust_radius": 0.08,
        "trust_min": 0.004,
        "trust_max": 0.16,
        "max_cycles": 100,
        "thresh": "gau_tight",
        "dump": True,
    }


def test_path_seed_rsirfo_fails_closed_without_seed_hessian(
    tmp_path: Path,
) -> None:
    result = PysisyphusSaddleEngine().search_ts(
        _reaction(),
        _species("opt_tma_hf2_neutral", NEUTRAL),
        _species("opt_tma_hf2_shared_proton", SHARED_PROTON),
        {
            "path_strategy": "saddle_refinement",
            "saddle_seed_evidence_validated": True,
            "saddle_seed_candidate": {"xyz_path": str(NEUTRAL)},
            "require_supplied_seed_hessian": True,
            "program": "nwchem",
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
            "required_program_version": "7.2.3",
        },
        tmp_path / "rsirfo",
    )

    assert result.success is False
    assert result.artifacts[-1].status.category == (
        "pysisyphus_saddle_seed_hessian_missing"
    )


def test_rsirfo_hessian_cache_requires_matching_seed_and_pes(
    tmp_path: Path, monkeypatch
) -> None:
    h5py = pytest.importorskip("h5py")
    seed = read_xyz(NEUTRAL)
    cache = tmp_path / "seed_hessian.h5"
    with h5py.File(cache, "w") as handle:
        handle.create_dataset(
            "coords3d", data=seed.coords * saddle_module._BOHR_PER_ANGSTROM
        )
        handle.create_dataset(
            "hessian", data=np.eye(3 * len(seed.symbols), dtype=float)
        )
    evidence_dir = tmp_path / "raw"
    evidence_dir.mkdir()
    monkeypatch.setattr(
        saddle_module,
        "_nwchem_qcengine_evidence",
        lambda *args, **kwargs: {"accepted": True, "raw_outputs": []},
    )
    method = {
        "ts_hessian_init_h5": str(cache),
        "ts_hessian_evidence_dir": str(evidence_dir),
        "required_program_version": "7.2.3",
        "basis": "def2-svpd",
    }

    hessian_init, evidence = saddle_module._validated_hessian_initialization(
        NEUTRAL, method
    )

    assert hessian_init == str(cache.resolve())
    assert evidence["accepted"] is True
    assert evidence["coordinate_max_abs_delta_bohr"] == pytest.approx(0.0)

    with h5py.File(cache, "r+") as handle:
        handle["coords3d"][0, 0] += 0.01
    with pytest.raises(ValueError, match="coordinate, shape, or PES"):
        saddle_module._validated_hessian_initialization(NEUTRAL, method)


def test_endpoint_dimer_seed_moves_off_minimum_toward_target(
    tmp_path: Path,
) -> None:
    seed = make_displaced_reaction_coordinate_seed(
        _reaction().data,
        NEUTRAL,
        SHARED_PROTON,
        tmp_path / "seed.xyz",
        displacement_A=0.15,
    )

    assert seed["geometry_qc"]["geometry_sane"] is True
    assert seed["direction_sign"] == 1.0
    assert seed["endpoint_coordinate"] < seed["seed_coordinate"]
    assert seed["seed_coordinate"] < seed["target_coordinate"]

    reverse = make_displaced_reaction_coordinate_seed(
        _reaction().data,
        SHARED_PROTON,
        NEUTRAL,
        tmp_path / "reverse_seed.xyz",
        displacement_A=0.15,
    )
    assert reverse["direction_sign"] == -1.0
    assert reverse["target_coordinate"] < reverse["seed_coordinate"]
    assert reverse["seed_coordinate"] < reverse["endpoint_coordinate"]


def test_dimer_candidate_requires_rotation_convergence_and_negative_curvature(
    tmp_path: Path,
) -> None:
    (tmp_path / "ts_opt.xyz").write_text(
        "1\ncandidate\nH 0 0 0\n", encoding="utf-8"
    )
    (tmp_path / "ts_optimizer.log").write_text(
        "Force convergence overachieved\nConverged!\n", encoding="utf-8"
    )
    (tmp_path / "dimer.log").write_text(
        "Doing dimer rotations\n"
        "00: rms(rot_force)=0.001003 C= 0.017020\n"
        "Rotation did not yield a negative curvature. "
        "Restoring previous unrotated N.\n",
        encoding="utf-8",
    )

    evidence = _dimer_convergence_evidence(tmp_path)

    assert evidence["optimizer_converged"] is True
    assert evidence["rotation_converged"] is False
    assert evidence["negative_curvature"] is False
    assert evidence["accepted"] is False

    (tmp_path / "dimer.log").write_text(
        "Doing dimer rotations\n"
        "Dimer rotation converged in 3 cycles.\n"
        "02: rms(rot_force)=0.000020 C=-0.001500\n",
        encoding="utf-8",
    )
    accepted = _dimer_convergence_evidence(tmp_path)
    assert accepted["final_curvature_hartree_per_bohr2"] == -0.0015
    assert accepted["accepted"] is True


def test_endpoint_dimer_failure_is_auditable_and_bounded(
    tmp_path: Path, monkeypatch
) -> None:
    def fake_command(command, *, cwd, **kwargs):
        workdir = Path(cwd)
        (workdir / "ts_opt.xyz").write_text(
            NEUTRAL.read_text(encoding="utf-8"), encoding="utf-8"
        )
        (workdir / "ts_optimizer.log").write_text(
            "Number of cycles exceeded!\n", encoding="utf-8"
        )
        stdout = workdir / kwargs["stdout_name"]
        stderr = workdir / kwargs["stderr_name"]
        stdout.write_text("dimer stopped\n", encoding="utf-8")
        stderr.write_text("", encoding="utf-8")
        return CommandResult(
            command=list(command),
            cwd=str(workdir),
            returncode=0,
            stdout_path=str(stdout),
            stderr_path=str(stderr),
            duration_s=0.1,
            executable=str(command[0]),
        )

    monkeypatch.setattr(
        "hfauto.backends.ts.pysisyphus.resolve_executable",
        lambda *args, **kwargs: "/usr/bin/true",
    )
    monkeypatch.setattr(
        "hfauto.backends.ts.pysisyphus.nwchem_environment",
        lambda *args, **kwargs: {},
    )
    monkeypatch.setattr(
        "hfauto.backends.ts.pysisyphus.run_command", fake_command
    )
    monkeypatch.setattr(
        "hfauto.backends.ts.pysisyphus._nwchem_qcengine_evidence",
        lambda *args, **kwargs: {
            "accepted": True,
            "raw_output_count": 1,
            "raw_outputs": [],
        },
    )
    result = PysisyphusEngine().search_ts(
        _reaction(),
        _species("opt_tma_hf2_neutral", NEUTRAL),
        _species("opt_tma_hf2_shared_proton", SHARED_PROTON),
        {
            "path_strategy": "endpoint_local_saddle_search",
            "allow_subprocess": True,
            "program": "nwchem",
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
            "required_program_version": "7.2.3",
        },
        tmp_path / "dimer",
    )

    attempts = [
        artifact
        for artifact in result.artifacts
        if artifact.artifact_type == "saddle_attempt"
    ]
    assert result.success is False
    assert len(attempts) == 2
    assert {attempt.data["endpoint_label"] for attempt in attempts} == {
        "reactant",
        "product",
    }
    assert all(
        attempt.data["diagnosis"] == "saddle_optimizer_failed"
        for attempt in attempts
    )


def test_basin_match_allows_symmetry_but_rejects_mirror_images(
    tmp_path: Path,
) -> None:
    reference = tmp_path / "reference.xyz"
    permuted = tmp_path / "permuted.xyz"
    reference.write_text(
        "5\nreference\nC 0 0 0\nH 1.1 0 0\nH -0.3 1.0 0\n"
        "H -0.2 -0.4 0.9\nH -0.5 -0.6 -0.8\n",
        encoding="utf-8",
    )
    permuted.write_text(
        "5\npermuted\nC 0 0 0\nH -0.2 -0.4 0.9\nH -0.5 -0.6 -0.8\n"
        "H 1.1 0 0\nH -0.3 1.0 0\n",
        encoding="utf-8",
    )
    symmetric_match = match_geometry_to_basin(
        permuted,
        reference,
        ordered_rmsd_threshold_A=1.0e-6,
        permutation_rmsd_threshold_A=1.0e-6,
    )
    assert symmetric_match["accepted"] is True
    assert symmetric_match["match_method"] == "covalent_graph_isomorphism_rmsd"

    chiral = tmp_path / "chiral.xyz"
    mirror = tmp_path / "mirror.xyz"
    chiral.write_text(
        "5\nchiral\nC 0 0 0\nH 0.64 0.64 0.64\nF -0.75 -0.75 0.75\n"
        "Cl -0.98 0.98 -0.98\nBr 1.04 -1.04 -1.04\n",
        encoding="utf-8",
    )
    mirror.write_text(
        "5\nmirror\nC 0 0 0\nH -0.64 0.64 0.64\nF 0.75 -0.75 0.75\n"
        "Cl 0.98 0.98 -0.98\nBr -1.04 -1.04 -1.04\n",
        encoding="utf-8",
    )
    mirror_match = match_geometry_to_basin(
        mirror,
        chiral,
        ordered_rmsd_threshold_A=0.10,
        permutation_rmsd_threshold_A=0.10,
    )
    assert mirror_match["distance_spectrum_rmsd_A"] == 0.0
    assert mirror_match["accepted"] is False


def test_ts_stage_obeys_reaction_case_authorization(tmp_path: Path) -> None:
    reaction = _reaction()
    manifest = Manifest.new(run_id="test", stage="reaction-plan")
    manifest.extend(
        [
            reaction,
            _species("opt_tma_hf2_neutral", NEUTRAL),
            _species("opt_tma_hf2_shared_proton", SHARED_PROTON),
            Artifact(
                artifact_id="reaction_case_rxn",
                artifact_type="reaction_case",
                data={
                    "reaction_id": reaction.data["reaction_id"],
                    "ts_search_allowed": False,
                    "next_action": "implement_path_strategy",
                },
            ),
        ]
    )

    result = TSSearchStage().run(
        manifest,
        {
            "engine": "dummy",
            "require_reaction_case": True,
            "require_validated_dft_minima": True,
        },
        StageContext(
            out_dir=tmp_path / "ts",
            run_id="test",
            global_config={"mode": "development"},
        ),
    )

    failure = result.latest_artifacts("ts_result")[-1]
    assert failure.status.category == "reaction_case_not_ts_ready"


def test_irc_stage_obeys_reaction_case_authorization(tmp_path: Path) -> None:
    manifest = Manifest.new(run_id="test", stage="ts-search")
    manifest.extend(
        [
            Artifact(
                artifact_id="rxn_with_ts",
                artifact_type="reaction_validated",
                data={
                    "reaction_id": "rxn",
                    "ts_species_id": "ts",
                    "reactant_species_id": "reactant",
                    "product_species_id": "product",
                    "mechanism_family": "generic_rearrangement",
                },
                qc={"ts_validated_by_frequency": True},
            ),
            Artifact(
                artifact_id="reaction_case_rxn",
                artifact_type="reaction_case",
                data={
                    "reaction_id": "rxn",
                    "irc_allowed": False,
                    "next_action": "implement_saddle_refiner",
                },
            ),
        ]
    )

    result = IRCStage().run(
        manifest,
        {"engine": "dummy", "require_reaction_case": True},
        StageContext(
            out_dir=tmp_path / "irc",
            run_id="test",
            global_config={"mode": "development"},
        ),
    )

    failure = result.latest_artifacts("irc")[-1]
    assert failure.status.category == "reaction_case_not_irc_ready"

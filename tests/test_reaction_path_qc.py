from pathlib import Path

import numpy as np

from hfauto.chemistry.reaction_path_qc import (
    endpoint_pair_match_qc,
    estimate_reaction_mode_overlap,
    proton_transfer_coordinate_vector,
    reaction_coordinate_vector,
)
from hfauto.chemistry.xyz import XYZ, write_xyz
from hfauto.core.frequency_qc import DEFAULT_IMAGINARY_FREQUENCY_CUTOFF_CM1
from hfauto.core.qc import ts_qc


def _write(path: Path, h_x: float, f_x: float) -> Path:
    return write_xyz(
        XYZ(
            ["N", "H", "F"],
            np.array([[0.0, 0.0, 0.0], [h_x, 0.0, 0.0], [f_x, 0.0, 0.0]]),
            path.stem,
        ),
        path,
    )


def _coordinate_data() -> dict:
    return {
        "reaction_coordinate": {
            "type": "distance_difference",
            "atoms": {"base_atom": 0, "transfer_h": 1, "leaving_f": 2},
        }
    }


def test_mode_overlap_fails_closed_without_literal_displacement(tmp_path):
    ts_xyz = _write(tmp_path / "ts.xyz", 1.20, 2.40)

    overlap, method = estimate_reaction_mode_overlap(
        _coordinate_data(),
        ts_xyz,
        n_imag=1,
        imag_freq_cm1=-420.0,
    )

    assert overlap == 0.0
    assert method == "normal_mode_displacement_unavailable"


def test_mode_overlap_projects_literal_cartesian_eigenvector(tmp_path):
    ts_xyz = _write(tmp_path / "ts.xyz", 1.20, 2.40)
    coordinate_vector = proton_transfer_coordinate_vector(_coordinate_data(), ts_xyz)
    assert coordinate_vector is not None

    overlap, method = estimate_reaction_mode_overlap(
        _coordinate_data(),
        ts_xyz,
        n_imag=1,
        imag_freq_cm1=-420.0,
        mode_displacements=coordinate_vector,
    )

    assert overlap == 1.0
    assert method == "cartesian_normal_mode_projection"


def test_nwchem_mode_projection_uses_common_mass_weighted_space(tmp_path):
    xyz = write_xyz(
        XYZ(
            ["N", "H", "H", "H"],
            np.array(
                [
                    [0.0, 0.0, 0.35],
                    [0.94, 0.0, 0.0],
                    [-0.47, 0.814, 0.0],
                    [-0.47, -0.814, 0.0],
                ]
            ),
            "ammonia",
        ),
        tmp_path / "ammonia.xyz",
    )
    data = {
        "reaction_coordinate": {
            "terms": [
                {"kind": "dihedral", "atoms": [0, 1, 2, 3], "coefficient": 1.0}
            ]
        }
    }
    coordinate_vector = reaction_coordinate_vector(data, xyz)
    assert coordinate_vector is not None
    nwchem_displacement = coordinate_vector.copy()

    overlap, method = estimate_reaction_mode_overlap(
        data,
        xyz,
        n_imag=1,
        imag_freq_cm1=-750.0,
        mode_displacements=nwchem_displacement,
        mode_component_units="amu^-1/2",
    )

    assert overlap > 0.999
    assert method == "mass_weighted_normal_mode_projection"


def test_shallow_but_significant_imaginary_mode_can_validate_ts(tmp_path):
    """A real -5 cm-1 mode is not rejected by the former -100 cm-1 rule."""

    ts_xyz = _write(tmp_path / "shallow_ts.xyz", 1.20, 2.40)
    coordinate_vector = proton_transfer_coordinate_vector(_coordinate_data(), ts_xyz)
    assert coordinate_vector is not None

    overlap, method = estimate_reaction_mode_overlap(
        _coordinate_data(),
        ts_xyz,
        n_imag=1,
        imag_freq_cm1=-5.0,
        mode_displacements=coordinate_vector,
    )
    qc = ts_qc(1, -5.0, overlap)

    assert DEFAULT_IMAGINARY_FREQUENCY_CUTOFF_CM1 == -1.0
    assert method == "cartesian_normal_mode_projection"
    assert overlap == 1.0
    assert qc["ts_imaginary_frequency_gate_passed"] is True
    assert qc["ts_target_mode_projection_gate_passed"] is True
    assert qc["ts_validated_by_frequency"] is True
    assert (
        ts_qc(
            1,
            -5.0,
            overlap,
            imaginary_frequency_cutoff_cm1=-10.0,
        )["ts_validated_by_frequency"]
        is False
    )


def test_numerical_noise_only_cannot_validate_ts(tmp_path):
    ts_xyz = _write(tmp_path / "noise_ts.xyz", 1.20, 2.40)
    coordinate_vector = proton_transfer_coordinate_vector(_coordinate_data(), ts_xyz)
    assert coordinate_vector is not None

    overlap, method = estimate_reaction_mode_overlap(
        _coordinate_data(),
        ts_xyz,
        n_imag=1,
        imag_freq_cm1=-0.8,
        mode_displacements=coordinate_vector,
    )
    qc = ts_qc(1, -0.8, 1.0)

    assert overlap == 0.0
    assert method == "frequency_gate_failed"
    assert qc["ts_imaginary_frequency_gate_passed"] is False
    assert qc["ts_validated_by_frequency"] is False


def test_legacy_mode_heuristic_never_invents_fixed_overlap_without_coordinate(
    tmp_path,
):
    ts_xyz = _write(tmp_path / "ts.xyz", 1.20, 2.40)

    overlap, method = estimate_reaction_mode_overlap(
        {},
        ts_xyz,
        n_imag=1,
        imag_freq_cm1=-420.0,
        allow_legacy_geometry_heuristic=True,
    )

    assert overlap == 0.0
    assert method == "legacy_geometry_coordinate_unavailable"


def test_irc_endpoint_pair_requires_identity_invariant_proton_basin(tmp_path):
    reactant = _write(tmp_path / "reactant.xyz", 1.60, 2.53)
    product = _write(tmp_path / "product.xyz", 1.05, 2.45)
    wrong_product = _write(tmp_path / "wrong_product.xyz", 1.55, 2.48)

    strict = endpoint_pair_match_qc(
        wrong_product,
        reactant,
        reactant,
        product,
        _coordinate_data(),
    )
    legacy = endpoint_pair_match_qc(
        wrong_product,
        reactant,
        reactant,
        product,
        _coordinate_data(),
        require_identity_invariant_geometry=False,
    )

    assert strict["irc_validated"] is False
    assert (
        "observed_endpoint_is_neutral_complex"
        in strict["forward_endpoint_qc"]["endpoint_match_reasons"]
    )
    assert legacy["irc_validated"] is True


def test_irc_endpoint_pair_requires_identity_invariant_q_agreement(tmp_path):
    reactant = _write(tmp_path / "reactant.xyz", 1.60, 2.53)
    product = _write(tmp_path / "product.xyz", 1.05, 2.45)
    shifted_product = _write(tmp_path / "shifted_product.xyz", 1.05, 2.90)

    qc = endpoint_pair_match_qc(
        shifted_product,
        reactant,
        reactant,
        product,
        _coordinate_data(),
        q_tolerance_A=0.30,
    )

    assert qc["irc_validated"] is False
    assert (
        "identity_invariant_q_mismatch"
        in qc["forward_endpoint_qc"]["endpoint_match_reasons"]
    )


def test_irc_endpoint_pair_accepts_exact_validated_basin_pair(tmp_path):
    reactant = _write(tmp_path / "reactant.xyz", 1.60, 2.53)
    product = _write(tmp_path / "product.xyz", 1.05, 2.45)

    qc = endpoint_pair_match_qc(
        product,
        reactant,
        reactant,
        product,
        _coordinate_data(),
    )

    assert qc["irc_validated"] is True
    assert qc["endpoint_orientation"] == "forward_product_backward_reactant"

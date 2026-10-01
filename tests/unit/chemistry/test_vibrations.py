"""Projected harmonic analysis (design §5.5, chem 6, CH-26)."""

import math

import numpy as np
import pytest

from hfauto.chemistry import vibrations as vib
from hfauto.chemistry.elements import mass
from hfauto.chemistry.identity import BASIN_DE_HARTREE
from hfauto.core.constants import AMU_TO_ME, BOHR_TO_ANGSTROM, CM1_TO_HARTREE

WATER = (
    ["O", "H", "H"], np.array([[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]])
)


def _hcn(angle_deg: float) -> np.ndarray:
    t = math.radians(180.0 - angle_deg)
    h = [-1.066 * math.cos(t), 1.066 * math.sin(t), 0.0]
    return np.array([[0.0, 0.0, 0.0], [1.153, 0.0, 0.0], h])


def test_projection_removes_rotational_contamination():
    symbols, x = WATER
    external = vib.external_basis(symbols, x)
    assert external.shape == (9, 6)
    rng = np.random.default_rng(1)
    q, _ = np.linalg.qr(np.hstack([external, rng.standard_normal((9, 3))]))
    internal = q[:, 6:]
    rotations = external[:, -2:]  # smallest singular values: rotations
    h_mw = internal @ np.diag([-0.05, 0.3, 0.6]) @ internal.T - 0.02 * rotations @ rotations.T
    root = np.sqrt(np.repeat([mass(s) for s in symbols], 3))
    freqs, modes, k = vib.projected_frequencies(h_mw * np.outer(root, root), symbols, x)
    assert np.count_nonzero(np.linalg.eigvalsh(h_mw) < -1e-8) == 3  # unprojected
    assert k == 6 and np.count_nonzero(freqs < 0) == 1
    factor = 1.0 / (math.sqrt(AMU_TO_ME) * CM1_TO_HARTREE)
    expected = np.sign([-0.05, 0.3, 0.6]) * np.sqrt([0.05, 0.3, 0.6]) * factor
    assert freqs == pytest.approx(expected, rel=1e-10)
    assert np.linalg.norm(modes, axis=1) == pytest.approx(1.0)


@pytest.mark.parametrize(("angle", "k"), [(180.0, 5), (179.99, 5), (179.0, 6)])
def test_linearity_from_svd_rank(angle, k):  # without a point group: engines, trials
    assert vib.external_basis(["C", "N", "H"], _hcn(angle)).shape[1] == k


@pytest.mark.parametrize(("angle", "linear"), [(179.99, False), (179.0, True), (180.0, True)])
def test_a_point_groups_linearity_sets_the_external_count(angle, linear):
    """G5-2: with ``linear`` given, k is 5 or 6 by definition, whatever the numerical rank."""
    x = _hcn(angle)
    h = np.diag(np.linspace(0.1, 0.9, 9))
    freqs, modes, k = vib.projected_frequencies(h, ["C", "N", "H"], x, linear=linear)
    assert k == (5 if linear else 6) and len(freqs) == len(modes) == 9 - k
    assert vib.projected_frequencies(np.zeros((3, 3)), ["Ar"], np.zeros((1, 3)), linear=False
                                     )[2] == 3  # an atom: no mode


def test_rotational_constants_of_linear_hcn():
    a, b, c = vib.rotational_constants_ghz(["C", "N", "H"], _hcn(180.0), linear=True)
    assert a == 0.0 and b == pytest.approx(c) and b == pytest.approx(44.5, abs=0.2)
    bent = vib.rotational_constants_ghz(["C", "N", "H"], _hcn(170.0), linear=False)
    assert bent[0] > 1000.0 and bent[1] > bent[2] > 0.0


def _four_atoms() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Coordinates, rigid motions (12, 6, unweighted) and an orthonormal internal basis (12, 6)."""
    x = np.array([[0.0, 0.0, 0.0], [1.45, 0.0, 0.0], [-0.3, 0.9, 0.1], [1.8, 0.2, 0.9]])
    c = x - x.mean(axis=0)
    rigid = [np.tile(a, 4) for a in np.eye(3)] + [np.cross(a, c).ravel() for a in np.eye(3)]
    noise = np.random.default_rng(3).standard_normal((12, 6))
    q, _ = np.linalg.qr(np.hstack([np.array(rigid).T, noise]))  # q[:, :6] spans rigid motion
    return x, q[:, :6], q[:, 6:]


def _model_input() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Two negative modes, a nearly degenerate soft pair (a rotor) and a rigid-motion leak."""
    x, rigid, internal = _four_atoms()
    curvature = np.array([-0.2, -0.1, 2.8e-3, 2.9e-3, 0.4, 0.6])
    h = internal @ np.diag(curvature) @ internal.T - 0.3 * np.outer(rigid[:, 5], rigid[:, 5])
    positive = internal @ np.diag(np.abs(curvature)) @ internal.T  # |λ| above the 1e-3 floor
    return x, rigid, internal, h, positive


@pytest.mark.parametrize("mix", [(1.0, 0, 0, 0, 0.2, 0),  # κ = d̂ᵀH₊d̂ = 0.21
                                 (0, 0, 1.0, 1.0, 0, 0),  # the rotor pair: κ = KAPPA_MIN
                                 (0, 0.5, 0.3, 0, 0, 1.0)])
def test_shape_hessian_puts_the_only_negative_curvature_along_the_direction(mix):
    """G1-P1: exactly one negative eigenvalue, −κ along d̂ (the direction's internal part), even
    for a direction between two nearly degenerate modes; the rest is P·H₊·P."""
    x, rigid, internal, h, positive = _model_input()
    direction = internal @ np.array(mix) + 0.5 * rigid[:, 0]
    d = internal @ np.array(mix) / np.linalg.norm(mix)
    model = vib.shape_hessian(h, x, direction)
    values, vectors = np.linalg.eigh(model)
    kappa = max(d @ positive @ d, vib.KAPPA_MIN)
    assert np.count_nonzero(values < -1e-8) == 1 and kappa >= 0.05
    assert values[0] == pytest.approx(-kappa) and abs(vectors[:, 0] @ d) == pytest.approx(1.0)
    p = np.eye(12) - np.outer(d, d)
    assert model + kappa * np.outer(d, d) == pytest.approx(p @ positive @ p, abs=1e-10)
    assert np.abs(model @ rigid).max() < 1e-10
    with pytest.raises(ValueError, match="no internal component"):
        vib.shape_hessian(h, x, rigid[:, 2])


def test_without_a_direction_the_model_is_positive_definite_on_the_same_modes():
    x, rigid, internal, h, positive = _model_input()
    model = vib.shape_hessian(h, x)
    assert model == pytest.approx(positive, abs=1e-10)  # |λ|, no sign flip
    assert np.abs(model @ rigid).max() < 1e-10  # rigid motions projected out
    assert np.linalg.eigvalsh(internal.T @ model @ internal).min() >= 2.8e-3 - 1e-12
    soft = internal @ np.diag([-0.2, -0.1, 1e-6, 0.4, 0.5, 0.6]) @ internal.T
    assert np.linalg.eigvalsh(internal.T @ vib.shape_hessian(soft, x) @ internal).min() == (
        pytest.approx(1e-3))  # the floor


def _mass_weighted_model(curvatures, gamma, rigid=0.0):
    """(H Eh/bohr², g Eh/bohr) of water with mass-weighted internal curvatures ``curvatures`` and
    mass-weighted internal gradient components ``gamma`` (on one internal basis), plus ``rigid``
    times a mass-weighted rigid motion in g."""
    symbols, x = WATER
    external = vib.external_basis(symbols, x)
    q, _ = np.linalg.qr(np.hstack([external, np.random.default_rng(7).standard_normal((9, 3))]))
    internal = q[:, 6:]
    root = np.sqrt(np.repeat([mass(s) for s in symbols], 3))
    h = (internal @ np.diag(curvatures) @ internal.T) * np.outer(root, root)
    g = (internal @ np.asarray(gamma) + rigid * external[:, 0]) * root
    return symbols, x, h, g


def test_the_stationarity_gap_weighs_a_residual_gradient_by_the_curvature_of_its_mode():
    """G5-P1: ΔE_N = Σ γ_i²/(2|λ_i|) over the modes of projected_frequencies, rigid motions out.
    The same residual gradient is a shoulder (W3 1ee97b40: 2.7e-4 Eh) along a soft mode and a
    stationary point (VAL7 minima: ≤ 1.2e-5 Eh) along stiff ones; |λ| of an imaginary mode."""
    shoulder = _mass_weighted_model([2e-5, 0.2, 0.5], [1e-4, 1e-5, 1e-5], rigid=0.3)
    stiff = _mass_weighted_model([2e-5, 0.2, 0.5], [0.0, 1e-4, 1e-5], rigid=0.3)
    imaginary = _mass_weighted_model([-2e-5, 0.2, 0.5], [1e-4, 1e-5, 1e-5])
    gap = 1e-8 / 4e-5 + 1e-10 / 0.4 + 1e-10 / 1.0
    assert vib.stationarity_gap(shoulder[2], shoulder[3], *shoulder[:2]) == pytest.approx(gap)
    assert vib.stationarity_gap(imaginary[2], imaginary[3], *imaginary[:2]) == pytest.approx(gap)
    assert gap > BASIN_DE_HARTREE > 1e3 * vib.stationarity_gap(stiff[2], stiff[3], *stiff[:2])


def test_the_newton_step_descends_every_internal_mode_of_the_positive_model():
    """−H₊⁺g on shape_hessian's modes (|λ| floored at 1e-3): the quadratic model's minimum
    along a positive mode, downhill along a negative one, rigid motions untouched."""
    x, rigid, internal, h, _ = _model_input()
    gamma = np.array([1e-3, -2e-3, 1e-4, 0.0, 4e-3, -6e-3])
    step = vib.newton_step(h, internal @ gamma + 0.1 * rigid[:, 1], x)
    assert step.shape == (4, 3)
    curvature = np.maximum(np.abs([-0.2, -0.1, 2.8e-3, 2.9e-3, 0.4, 0.6]), 1e-3)
    q = step.ravel() / BOHR_TO_ANGSTROM
    assert internal.T @ q == pytest.approx(-gamma / curvature, abs=1e-12)
    assert np.abs(rigid.T @ q).max() < 1e-12
    soft = internal @ np.diag([1e-6, 0.1, 0.2, 0.3, 0.4, 0.5]) @ internal.T
    assert internal[:, 0] @ vib.newton_step(soft, internal[:, 0] * 1e-3, x).ravel() == (
        pytest.approx(-BOHR_TO_ANGSTROM))  # 1e-3 / max(1e-6, 1e-3) bohr


def test_canonical_npy(tmp_path):
    h = np.diag([0.5] + [0.0] * 5)  # Eh/bohr², atom 0 along x
    out = vib.to_canonical_npy(h, tmp_path / "freq" / "hessian.npy")
    assert np.array_equal(np.load(out), h)

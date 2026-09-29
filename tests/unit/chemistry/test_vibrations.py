"""Projected harmonic analysis (design §5.5, chem 6, CH-26)."""

import math

import numpy as np
import pytest

from hfauto.chemistry import vibrations as vib
from hfauto.chemistry.elements import mass
from hfauto.core.constants import AMU_TO_ME, CM1_TO_HARTREE

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
def test_linearity_from_svd_rank(angle, k):
    assert vib.external_basis(["C", "N", "H"], _hcn(angle)).shape[1] == k


def test_rotational_constants_of_linear_hcn():
    a, b, c = vib.rotational_constants_ghz(["C", "N", "H"], _hcn(180.0))
    assert a == 0.0 and b == pytest.approx(c) and b == pytest.approx(44.5, abs=0.2)


def _four_atoms() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Coordinates, rigid motions (12, 6, unweighted) and an orthonormal internal basis (12, 6)."""
    x = np.array([[0.0, 0.0, 0.0], [1.45, 0.0, 0.0], [-0.3, 0.9, 0.1], [1.8, 0.2, 0.9]])
    c = x - x.mean(axis=0)
    rigid = [np.tile(a, 4) for a in np.eye(3)] + [np.cross(a, c).ravel() for a in np.eye(3)]
    noise = np.random.default_rng(3).standard_normal((12, 6))
    q, _ = np.linalg.qr(np.hstack([np.array(rigid).T, noise]))  # q[:, :6] spans rigid motion
    return x, q[:, :6], q[:, 6:]


@pytest.mark.parametrize("reaction", [0, 1, 4])
def test_shape_hessian_leaves_only_the_reaction_mode_negative(reaction):
    x, rigid, internal = _four_atoms()
    curvature = np.array([-0.2, -0.1, 1e-6, 0.4, 0.5, 0.6])  # two negative modes, one soft
    h = internal @ np.diag(curvature) @ internal.T - 0.3 * np.outer(rigid[:, 5], rigid[:, 5])
    mode = internal[:, reaction] + 0.2 * internal[:, 3] + 0.5 * rigid[:, 0]  # mixed
    values, vectors = np.linalg.eigh(vib.shape_hessian(h, x, mode))
    assert np.count_nonzero(values < -1e-8) == 1  # inertia: one negative eigenvalue
    assert abs(vectors[:, 0] @ internal[:, reaction]) == pytest.approx(1.0)  # mode chosen
    expected = np.maximum(np.abs(curvature), 1e-3)  # |λ| with the floor, rigid motions at 0
    expected[reaction] *= -1.0
    assert values == pytest.approx(sorted([*expected, 0, 0, 0, 0, 0, 0]), abs=1e-10)


def test_without_a_reaction_mode_the_model_is_positive_definite_on_the_same_modes():
    x, rigid, internal = _four_atoms()
    curvature = np.array([-0.2, -0.1, 1e-6, 0.4, 0.5, 0.6])
    h = internal @ np.diag(curvature) @ internal.T - 0.3 * np.outer(rigid[:, 5], rigid[:, 5])
    model = vib.shape_hessian(h, x, None)
    expected = np.maximum(np.abs(curvature), 1e-3)  # |λ| with the floor, no sign flip
    assert model == pytest.approx(internal @ np.diag(expected) @ internal.T, abs=1e-10)
    assert np.abs(model @ rigid).max() < 1e-10  # rigid motions projected out
    assert np.linalg.eigvalsh(internal.T @ model @ internal).min() >= 1e-3 - 1e-12


def test_canonical_npy(tmp_path):
    h = np.diag([0.5] + [0.0] * 5)  # Eh/bohr², atom 0 along x
    out = vib.to_canonical_npy(h, tmp_path / "freq" / "hessian.npy")
    assert np.array_equal(np.load(out), h)

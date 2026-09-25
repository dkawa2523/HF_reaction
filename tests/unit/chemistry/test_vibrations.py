"""Projected harmonic analysis (design §5.5, chem 6, CH-26)."""

import math

import numpy as np
import pytest

from hfauto.chemistry import vibrations as vib
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
    root = np.sqrt(np.repeat([vib.ISOTOPIC_MASSES[s] for s in symbols], 3))
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


def test_curvature_and_canonical_npy(tmp_path):
    h = np.diag([0.5] + [0.0] * 5)  # Eh/bohr², atom 0 along x
    mode = np.array([[2.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    assert vib.curvature_along(h, mode) == pytest.approx(0.5 / 0.529177210903**2)
    out = vib.to_canonical_npy(h, tmp_path / "freq" / "hessian.npy")
    assert np.array_equal(np.load(out), h)
    with pytest.raises(ValueError):
        vib.projected_frequencies(h, ["Xx", "H"], np.zeros((2, 3)))

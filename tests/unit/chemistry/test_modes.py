"""Mode displacements, QRC amplitudes and mode-follow outcomes (design §5.5, chem 7)."""

import numpy as np
import pytest

from hfauto.chemistry import modes
from hfauto.core.constants import BOHR_TO_ANGSTROM


def _quadratic_pes(curvature: float) -> tuple[np.ndarray, np.ndarray]:
    """Hessian (Eh/bohr²) of 3 atoms with one mode of the given eigenvalue."""
    rng = np.random.default_rng(7)
    q, _ = np.linalg.qr(rng.standard_normal((9, 9)))
    eigenvalues = np.array([curvature] + list(np.linspace(0.1, 0.5, 8)))
    return q @ np.diag(eigenvalues) @ q.T, q[:, 0]


def test_qrc_amplitude_hits_target_energy():
    hessian, mode = _quadratic_pes(-0.08)
    s = modes.qrc_amplitude(hessian, mode, target_hartree=3e-4, bounds_A=(0.01, 5.0))
    plus, _ = modes.displace(np.zeros((3, 3)), mode, s)
    step_bohr = plus.ravel() / BOHR_TO_ANGSTROM
    assert 0.5 * step_bohr @ hessian @ step_bohr == pytest.approx(-3e-4)


def test_qrc_amplitude_is_clipped():
    stiff, mode = _quadratic_pes(-50.0)
    assert modes.qrc_amplitude(stiff, mode, target_hartree=3e-4) == 0.05
    flat, mode = _quadratic_pes(-1e-6)
    assert modes.qrc_amplitude(flat, mode, target_hartree=3e-4) == 0.4


def test_displace_and_overlap():
    coords = np.zeros((2, 3))
    mode = np.array([[0.0, 0.0, 2.0], [0.0, 1.0, 0.0]])
    plus, minus = modes.displace(coords, mode, 0.1)
    assert np.linalg.norm(plus, axis=1).max() == pytest.approx(0.1)
    assert np.allclose(plus, -minus)
    assert modes.overlap(mode, -mode) == pytest.approx(1.0)
    assert modes.overlap([1.0, 0.0], [0.0, 1.0]) == 0.0


@pytest.mark.parametrize(
    ("plus", "minus", "expected"),
    [
        ("x", "x", "replace"),
        ("a", "b", "ts_candidate"),
        ("a", None, "one_side"),
        ("a", "src", "one_side"),
        ("src", None, "same_as_source"),
        (None, None, "same_as_source"),
    ],
)
def test_classify_mode_follow(plus, minus, expected):
    assert modes.classify_mode_follow("src", plus, minus) == expected

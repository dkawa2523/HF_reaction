"""Mode displacements, QRC amplitudes and mode-follow outcomes (design §5.5, chem 7)."""

from pathlib import Path

import numpy as np
import pytest

from hfauto.backends.nwchem.output import hess_to_npy
from hfauto.chemistry import modes
from hfauto.chemistry.vibrations import projected_frequencies
from hfauto.chemistry.xyz import read_xyz
from hfauto.core.constants import BOHR_TO_ANGSTROM

G07 = Path(__file__).parents[2] / "golden" / "data" / "nwchem" / "G07"  # the HCN TS, -1131i


def _former_qrc_rule(hessian, mode, target):  # uᵀHu of the 1 Å mode (removed in S-B)
    u = np.reshape(mode, (-1, 3)) / np.linalg.norm(np.reshape(mode, (-1, 3)), axis=1).max()
    curvature = abs(u.ravel() @ hessian @ u.ravel()) / BOHR_TO_ANGSTROM**2
    return float(np.clip(np.sqrt(2.0 * target / curvature), 0.05, 0.4))


@pytest.mark.golden
def test_amplitude_hits_the_target_energy_of_the_former_qrc_rule(tmp_path):
    """κ = ω²·Σm|u|² equals uᵀHu for a normal mode, so QRC amplitudes do not change."""
    xyz = read_xyz(G07 / "final.xyz")
    hessian = np.load(hess_to_npy(G07 / "hfauto_job.hess", 3, tmp_path / "h.npy"))
    freqs, vectors, _ = projected_frequencies(hessian, xyz.symbols, xyz.coords)
    s = modes.amplitude(freqs[0], vectors[0], xyz.symbols)
    assert freqs[0] < -1000 and 0.05 < s < 0.4  # not clipped
    assert s == pytest.approx(_former_qrc_rule(hessian, vectors[0], 3e-4), rel=1e-9)
    plus, _ = modes.displace(np.zeros((3, 3)), vectors[0], s)
    step_bohr = plus.ravel() / BOHR_TO_ANGSTROM
    assert 0.5 * step_bohr @ hessian @ step_bohr == pytest.approx(-3e-4, rel=1e-6)


def test_amplitude_is_clipped_and_mass_weighted():
    mode = np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 0.0]])
    assert modes.amplitude(-3000.0, mode, ["H", "H"]) == 0.05  # stiff
    assert modes.amplitude(-5.0, mode, ["H", "H"]) == 0.4  # flat
    assert modes.amplitude(0.0, mode, ["H", "H"]) == 0.4
    light = modes.amplitude(-300.0, mode, ["H", "H"], bounds_A=(0.0, 9.0))
    heavy = modes.amplitude(-300.0, mode, ["I", "H"], bounds_A=(0.0, 9.0))
    assert light / heavy == pytest.approx(np.sqrt(126.904 / 1.008), rel=1e-3)


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

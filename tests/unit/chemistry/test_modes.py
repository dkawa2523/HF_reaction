"""Mode displacements, QRC amplitudes and mode-follow outcomes (design §5.5, chem 7)."""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from hfauto.backends.nwchem.output import hess_to_npy
from hfauto.chemistry import modes
from hfauto.chemistry.vibrations import projected_frequencies
from hfauto.chemistry.xyz import read_xyz
from hfauto.core.constants import BOHR_TO_ANGSTROM
from hfauto.core.evidence import Level

G07 = Path(__file__).parents[2] / "golden" / "data" / "nwchem" / "G07"  # the HCN TS, -1131i
H4 = ["H"] * 4
AXES = [np.eye(12)[3 * i + i % 3].reshape(4, 3) for i in range(4)]  # atom i along axis i % 3


def _freq(nu, imaginary, scf_tol=1e-7):
    """What the displacements read of a freq job: its level, frequencies and imaginary modes."""
    level = Level(program="fake", version="0", method="pbe0", charge=0, multiplicity=1,
                  scf_tol=scf_tol)
    return SimpleNamespace(level=level, frequencies_cm1=tuple(nu),
                           imaginary_modes=tuple(tuple(np.ravel(m)) for m in imaginary))


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
    light = modes.amplitude(-300.0, mode, ["H", "H"])  # κ = ω²m: one κ, one amplitude
    heavy = modes.amplitude(-300.0 / np.sqrt(126.904 / 1.008), mode, ["I", "H"])
    assert 0.05 < light < 0.4 and light == pytest.approx(heavy, rel=1e-3)


def test_qrc_step_is_one_mode_at_the_qrc_target_capped():
    """The bundled methods' 3 x qrc_drop (3e-5 Eh) lies below TARGET_HARTREE, so their QRC and
    mode-follow steps take TARGET_HARTREE; a looser SCF (qrc_drop 2e-3 Eh) takes 6e-3 Eh."""
    freq = _freq((-900.0, -100.0, 50.0), AXES[:2])
    s = modes.amplitude(-900.0, AXES[0], H4)
    assert np.allclose(modes.qrc_step(freq, 0, H4), s * AXES[0]) and 0.05 < s < 0.2
    loose = modes.qrc_step(_freq((-900.0,), AXES[:1], scf_tol=1e-4), 0, H4)
    assert np.abs(loose).max() == pytest.approx(modes.amplitude(-900.0, AXES[0], H4,
                                                                target_hartree=6e-3))
    assert np.abs(modes.qrc_step(freq, 0, H4, factor=2.0)).max() == pytest.approx(2 * s)
    assert np.abs(modes.qrc_step(freq, 0, H4, factor=8.0)).max() == modes.BOUNDS_A[1]


def test_off_saddle_leaves_a_third_order_saddle_along_both_other_modes():
    """Shown on a fake freq only (no run has met a third-order saddle): keeping the reaction
    mode, the push sums the two other modes below the threshold, each at its qrc_step, where
    the former push took the most negative of them only. A soft mode joins only below a lower
    threshold (a soft point), and the sum is capped at 0.4 A."""
    freq = _freq((-900.0, -600.0, -300.0, -30.0, 100.0), AXES)
    steps = [modes.qrc_step(freq, i, H4) for i in range(4)]
    assert np.allclose(modes.off_saddle(freq, H4, below_cm1=50.0, keep=0), steps[1] + steps[2])
    assert np.allclose(modes.off_saddle(freq, H4, below_cm1=50.0), sum(steps[:3]))
    assert np.allclose(modes.off_saddle(freq, H4, below_cm1=10.0, keep=0), sum(steps[1:]))
    flat = _freq((-60.0, -55.0), [AXES[0], AXES[0]])  # two flat modes moving one atom 0.4 A
    step = modes.off_saddle(flat, H4, below_cm1=50.0)
    assert np.allclose(step, modes.BOUNDS_A[1] * AXES[0])


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

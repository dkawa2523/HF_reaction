"""Path class, peak interpolation and settled string energies (design §5.5, CH-08)."""

import numpy as np
import pytest

from hfauto.chemistry import profile as prof


@pytest.mark.parametrize(
    ("energies", "expected"),
    [
        ([0.0, 0.4, 0.2, 0.9, 2.0, 3.0], "barrierless"),  # wiggle below the resolution
        ([0.0, 1.0, 2.0, 3.0], "barrierless"),
        ([0.0, 2.0, 5.0, 3.0, 1.0], "single"),
        ([0.0, 10.0, 9.99, 10.0, 0.0], "single"),  # shallow dip counts once
        ([0.0, 5.0, 1.0], "single"),  # three points: a low-level TS between two minima
        ([0.0, 0.5, 1.0], "barrierless"),
        ([0.0, 5.0, 2.0, 5.0, 0.0], "intermediate"),  # two peaks always have a well between
        ([0.0, -1.5, -3.0, -1.0, 1.0, 3.0, 5.0], "intermediate"),  # well below both ends
        ([0.0, 2.0, 0.5, 3.0, 6.0, 9.0, 10.0], "intermediate"),  # uphill two-step path
    ],
)
def test_classify_counts_hills_and_wells_alike(energies, expected):
    assert prof.classify(energies, resolution=1.0) == expected


def test_hei_interpolates_parabola_and_coordinates_at_the_given_peak():
    s = np.arange(6.0)
    energies = -((s - 2.3) ** 2)
    frames = [np.array([[x, 0.0, 0.0]]) for x in s]
    index, energy, coords = prof.hei(frames, energies, 2)
    assert index == pytest.approx(2.3) and energy == pytest.approx(0.0, abs=1e-12)
    assert coords == pytest.approx(np.array([[2.3, 0.0, 0.0]]))
    two = [0.0, 3.0, 1.0, 5.0, 0.0]  # the detected peak, not only the global maximum
    assert prof.hei(frames[:5], two, 1)[0] == pytest.approx(1.1)  # toward the higher side
    with pytest.raises(ValueError, match="interior"):
        prof.hei(frames, energies, 0)


def test_energies_settled():
    settled = [(0.0, 5.0, 1.0), (0.0, 3.0, 1.0), (0.0, 3.06, 1.0), (0.0, 3.0, 1.02)]
    assert prof.energies_settled(settled, 0.1)  # the first step lies outside the window
    assert not prof.energies_settled(settled, 0.05)
    assert not prof.energies_settled([*settled, (0.0, 3.2, 1.0)], 0.1)  # still moving
    assert not prof.energies_settled(settled[-2:], 0.1)  # fewer than 3 iterations
    assert not prof.energies_settled([], 0.1)


def test_spacing_and_tangent():
    frames = [np.array([[x, 0.0, 0.0]]) for x in (0.0, 0.1, 0.5)]
    assert prof.max_node_spacing(frames) == pytest.approx(0.4)
    assert prof.tangent(frames, 1) == pytest.approx(np.array([[1.0, 0.0, 0.0]]))
    assert prof.tangent(frames, 0) == pytest.approx(np.array([[1.0, 0.0, 0.0]]))

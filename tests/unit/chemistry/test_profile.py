"""Path shape, HEI, string convergence and stagnation (design §5.5, CH-08)."""

import numpy as np
import pytest

from hfauto.chemistry import profile as prof


@pytest.mark.parametrize(
    ("energies", "expected"),
    [
        ([0.0, 0.4, 0.2, 0.9, 2.0, 3.0], "monotonic"),  # wiggle below the resolution
        ([0.0, 2.0, 5.0, 3.0, 1.0], "single_max"),
        ([0.0, 10.0, 9.99, 10.0, 0.0], "single_max"),  # shallow dip counts once
        ([0.0, 5.0, 2.0, 5.0, 0.0], "multi_max"),
        ([0.0, 1.0, 2.0, 3.0], "monotonic"),
    ],
)
def test_shape(energies, expected):
    assert prof.shape(energies, resolution=1.0) == expected


def test_hei_interpolates_parabola_and_coordinates():
    s = np.arange(6.0)
    energies = -((s - 2.3) ** 2)
    frames = [np.array([[x, 0.0, 0.0]]) for x in s]
    index, energy, coords = prof.hei(frames, energies)
    assert index == pytest.approx(2.3) and energy == pytest.approx(0.0, abs=1e-12)
    assert coords == pytest.approx(np.array([[2.3, 0.0, 0.0]]))


def test_string_converged():
    assert prof.string_converged([5e-3, 2e-3, 9e-4, 8e-4, 7e-4])
    assert not prof.string_converged([5e-3, 7e-4, 8e-4, 7e-4])  # worsened within the last 3
    assert not prof.string_converged([2e-3])
    assert not prof.string_converged([])


def test_stagnated():
    assert prof.stagnated([1.0] + [0.5] * 11)
    assert not prof.stagnated([0.9**i for i in range(20)])
    assert prof.stagnated([], [1e-3] + [1e-8] * 5)
    assert not prof.stagnated([], [1e-8] * 4)


def test_spacing_and_tangent():
    frames = [np.array([[x, 0.0, 0.0]]) for x in (0.0, 0.1, 0.5)]
    assert prof.max_node_spacing(frames) == pytest.approx(0.4)
    assert prof.tangent(frames, 1) == pytest.approx(np.array([[1.0, 0.0, 0.0]]))
    assert prof.tangent(frames, 0) == pytest.approx(np.array([[1.0, 0.0, 0.0]]))

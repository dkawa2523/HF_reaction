"""Permutation-invariant identity with proper rotations only (design §5.5, chem 13)."""

import numpy as np
import pytest

from hfauto.chemistry import identity as idn

NH3_SYMBOLS = ["N", "H", "H", "H"]
NH3 = np.array([[0.0, 0.0, 0.38], [0.94, 0.0, 0.0], [-0.47, 0.814, 0.0], [-0.47, -0.814, 0.0]])
INVERTED = NH3 * [1.0, 1.0, -1.0]
CHFCLBR_SYMBOLS = ["C", "H", "F", "Cl", "Br"]
CHFCLBR = np.array(
    [[0.0, 0.0, 0.0], [0.63, 0.63, 0.63], [-0.8, -0.8, 0.8], [-1.0, 1.0, -1.0], [1.1, -1.1, -1.1]]
)


def _proper_rotation(seed: int) -> np.ndarray:
    q, _ = np.linalg.qr(np.random.default_rng(seed).standard_normal((3, 3)))
    return q * np.sign(np.linalg.det(q))


def test_ammonia_inversion_is_a_degenerate_pair():
    assert idn.compare_minima(NH3_SYMBOLS, NH3, INVERTED, -56.5, -56.5) == "same"
    assert idn.mapped_equivalent(NH3_SYMBOLS, NH3, INVERTED)
    assert not idn.mapped_equivalent(NH3_SYMBOLS, NH3, NH3)


def test_enantiomers_stay_distinct():
    mirror = CHFCLBR * [-1.0, 1.0, 1.0]
    assert idn.compare_minima(CHFCLBR_SYMBOLS, CHFCLBR, mirror, -1.0, -1.0) == "distinct"


def test_rotated_permuted_copy_is_recovered():
    order = [0, 3, 1, 2]
    moved = (NH3 @ _proper_rotation(3).T + [1.0, -2.0, 0.5])[order]
    rmsd, perm = idn.permutation_invariant_rmsd(NH3_SYMBOLS, NH3, moved)
    assert rmsd < 1e-6
    assert [order[i] for i in perm] == [0, 1, 2, 3]


def test_thresholds_and_assignment():
    assert idn.compare_minima(NH3_SYMBOLS, NH3, NH3, 0.0, 3e-5) == "ambiguous"
    squeezed = NH3 * [1.0, 1.0, 0.5]
    candidates = {"up": (NH3, 0.0), "flat": (squeezed, 0.0)}
    assert idn.assign(NH3_SYMBOLS, NH3 + 1e-4, 1e-6, candidates) == "up"
    assert idn.assign(NH3_SYMBOLS, NH3, 1e-3, candidates) is None  # energy off
    twins = {"a": (NH3, 0.0), "b": (NH3 + 0.004, 0.0)}
    assert idn.assign(NH3_SYMBOLS, NH3 + 0.002, 0.0, twins) is None  # not unique


def test_periodic_nearest_wraps():
    assert idn.periodic_nearest(355.0, [0.0, 120.0, 240.0]) == 0
    assert idn.periodic_nearest(170.0, [-170.0, 90.0]) == 0
    with pytest.raises(ValueError):
        idn.periodic_nearest(0.0, [])

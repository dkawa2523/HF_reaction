"""Permutation-invariant identity: mirror images are one basin, labels stay proper (§5.5)."""

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
MIRROR = CHFCLBR * [-1.0, 1.0, 1.0]
CH4 = np.array([[0, 0, 0], [0.63, 0.63, 0.63], [-0.63, -0.63, 0.63], [-0.63, 0.63, -0.63],
                [0.63, -0.63, -0.63]])


def _proper_rotation(seed: int) -> np.ndarray:
    q, _ = np.linalg.qr(np.random.default_rng(seed).standard_normal((3, 3)))
    return q * np.sign(np.linalg.det(q))


def test_ammonia_inversion_is_one_basin_and_a_degenerate_pair():
    assert idn.same_basin(NH3_SYMBOLS, NH3, INVERTED, -56.5, -56.5)
    assert idn.mapped_equivalent(NH3_SYMBOLS, NH3, INVERTED)  # labels differ: proper only
    assert not idn.mapped_equivalent(NH3_SYMBOLS, NH3, NH3)
    assert not idn.is_chiral(NH3_SYMBOLS, NH3)


def test_enantiomers_are_one_chiral_basin():
    assert idn.permutation_invariant_rmsd(CHFCLBR_SYMBOLS, CHFCLBR, MIRROR)[0] > 1.0  # proper
    moved = MIRROR @ _proper_rotation(5).T + [0.3, 0.0, -1.0]
    assert idn.same_basin(CHFCLBR_SYMBOLS, CHFCLBR, moved, -1.0, -1.0)
    assert idn.is_chiral(CHFCLBR_SYMBOLS, CHFCLBR) and idn.is_chiral(CHFCLBR_SYMBOLS, moved)
    assert idn.mapped_equivalent(CHFCLBR_SYMBOLS, CHFCLBR, MIRROR)  # R -> S is degenerate


def test_methane_is_achiral():
    assert not idn.is_chiral(["C", "H", "H", "H", "H"], CH4)


def test_basin_coords_take_the_atom_order_and_handedness_of_the_member():
    member = MIRROR @ _proper_rotation(7).T + 0.02  # the S member of an R basin
    x = idn.basin_coords(CHFCLBR_SYMBOLS, CHFCLBR, member)
    assert idn.mapped_rmsd(x, member) < 0.05 < idn.mapped_rmsd(x, CHFCLBR)
    order = [0, 2, 3, 1]  # an achiral basin is only relabelled
    x = idn.basin_coords(NH3_SYMBOLS, NH3, NH3[order] @ _proper_rotation(3).T + 0.03)
    assert np.allclose(x, NH3[order])


def test_rotated_permuted_copy_is_recovered():
    order = [0, 3, 1, 2]
    moved = (NH3 @ _proper_rotation(3).T + [1.0, -2.0, 0.5])[order]
    rmsd, perm = idn.permutation_invariant_rmsd(NH3_SYMBOLS, NH3, moved)
    assert rmsd < 1e-6
    assert [order[i] for i in perm] == [0, 1, 2, 3]


def test_one_criterion_energy_first_then_a_unique_structure():
    assert idn.same_basin(NH3_SYMBOLS, NH3, NH3, 0.0, 3e-5)  # 5e-5 Eh
    assert not idn.same_basin(NH3_SYMBOLS, NH3, NH3, 0.0, 6e-5)
    squeezed = NH3 * [1.0, 1.0, 0.5]
    assert not idn.same_basin(NH3_SYMBOLS, NH3, squeezed, 0.0, 0.0)
    candidates = {"up": (NH3, 0.0), "flat": (squeezed, 0.0)}
    assert idn.assign(NH3_SYMBOLS, NH3 + 1e-4, 1e-6, candidates) == "up"
    assert idn.assign(NH3_SYMBOLS, NH3, 1e-3, candidates) is None  # energy off
    twins = {"a": (NH3, 0.0), "b": (NH3 + 0.004, 0.0)}
    assert idn.assign(NH3_SYMBOLS, NH3 + 0.002, 0.0, twins) is None  # not unique
    twins["b"] = (NH3 + 0.004, 1e-3)  # another energy never competes as the runner-up
    assert idn.assign(NH3_SYMBOLS, NH3 + 0.002, 0.0, twins) == "a"


def test_periodic_nearest_wraps():
    assert idn.periodic_nearest(355.0, [0.0, 120.0, 240.0]) == 0
    assert idn.periodic_nearest(170.0, [-170.0, 90.0]) == 0
    with pytest.raises(ValueError):
        idn.periodic_nearest(0.0, [])

"""Hysteresis bonds, fragments, proton transfer and state labels (design §5.5, CH-07)."""

import numpy as np
import pytest

from hfauto.chemistry import topology as top

FHF = (["F", "H", "F"], np.array([[-1.14, 0.0, 0.0], [0.0, 0.0, 0.0], [1.14, 0.0, 0.0]]))
NH3_HF = ["N", "H", "H", "H", "H", "F"]
_NH3 = [[0.0, 0.0, 0.0], [0.94, 0.0, -0.38], [-0.47, 0.814, -0.38], [-0.47, -0.814, -0.38]]
NEUTRAL = np.array(_NH3 + [[0.0, 0.0, 1.75], [0.0, 0.0, 2.70]])  # N···H–F
ION_PAIR = np.array(_NH3 + [[0.0, 0.0, 1.05], [0.0, 0.0, 2.70]])  # N–H···F⁻


def test_symmetric_hydrogen_bond_is_one_fragment():
    assert top.fragments(*FHF) == ((0, 1, 2),)
    assert top.state_label(*FHF).startswith("F2H_")


def test_hysteresis_keeps_previous_state_in_the_band():
    x = np.array([[0.0, 0.0, 0.0], [1.3 * (0.31 + 0.57), 0.0, 0.0]])
    assert top.bonds(["H", "F"], x, previous={(0, 1)}) == {(0, 1)}
    assert top.bonds(["H", "F"], x, previous=frozenset()) == frozenset()
    assert top.bonds(["H", "F"], x) == {(0, 1)}


def test_proton_transfer_changes():
    formed, broken = top.bond_changes(NH3_HF, NEUTRAL, ION_PAIR)
    assert formed == {(0, 4)} and broken == {(4, 5)}
    assert top.transferred_hydrogens(NH3_HF, NEUTRAL, ION_PAIR) == 1
    assert top.transferred_hydrogens(NH3_HF, NEUTRAL, NEUTRAL) == 0
    assert top.proton_coordinate(NEUTRAL, 5, 4, 0) == pytest.approx(0.95 - 1.75)


def test_labile_hydrogens_and_acceptors():
    assert top.labile_hydrogens(NH3_HF, NEUTRAL) == (1, 2, 3, 4)
    assert top.acceptor_atoms(NH3_HF, NEUTRAL) == (0, 5)
    assert top.acceptor_atoms(NH3_HF, ION_PAIR) == (5,)


def test_state_label_is_permutation_invariant():
    order = [5, 3, 0, 4, 1, 2]
    permuted = top.state_label([NH3_HF[i] for i in order], NEUTRAL[order])
    assert permuted == top.state_label(NH3_HF, NEUTRAL)
    assert top.state_label(NH3_HF, ION_PAIR) != top.state_label(NH3_HF, NEUTRAL)


def test_vdw_radius_table():
    assert top.vdw_radius("H") == 1.20 and top.vdw_radius("N") == 1.55
    with pytest.raises(ValueError):
        top.vdw_radius("Xx")

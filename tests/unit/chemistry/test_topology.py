"""Hysteresis bonds, fragments, bond changes, WL atom classes and state labels (design §5.5,
CH-07)."""

import numpy as np

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
    assert top.bond_changes(NH3_HF, NEUTRAL, NEUTRAL) == (frozenset(), frozenset())


def test_labile_hydrogens_and_acceptors():
    assert top.labile_hydrogens(NH3_HF, NEUTRAL) == (1, 2, 3, 4)
    assert top.acceptor_atoms(NH3_HF, NEUTRAL) == (0, 5)
    assert top.acceptor_atoms(NH3_HF, ION_PAIR) == (5,)


def test_state_label_is_permutation_invariant():
    order = [5, 3, 0, 4, 1, 2]
    permuted = top.state_label([NH3_HF[i] for i in order], NEUTRAL[order])
    assert permuted == top.state_label(NH3_HF, NEUTRAL)
    assert top.state_label(NH3_HF, ION_PAIR) != top.state_label(NH3_HF, NEUTRAL)


def test_wl_classes_are_atom_equivalence_and_permutation_invariant():
    bonded = top.bonds(NH3_HF, NEUTRAL)
    classes = top.wl_classes(NH3_HF, bonded)
    assert classes[1] == classes[2] == classes[3] and len(set(classes)) == 4  # N, 3 H, H, F
    order = [5, 3, 0, 4, 1, 2]
    permuted = top.wl_classes([NH3_HF[i] for i in order], top.bonds(
        [NH3_HF[i] for i in order], NEUTRAL[order]))
    assert permuted == tuple(classes[i] for i in order)

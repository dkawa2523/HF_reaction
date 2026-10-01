"""Additive bond rule, fragments, bond changes resolved by one band (analysis X3), WL atom
classes and state labels (design §5.5, review §6.3)."""

import numpy as np
import pytest

from hfauto.chemistry import topology as top
from hfauto.chemistry.elements import covalent_radius

NH3_HF = ["N", "H", "H", "H", "H", "F"]
_NH3 = [[0.0, 0.0, 0.0], [0.94, 0.0, -0.38], [-0.47, 0.814, -0.38], [-0.47, -0.814, -0.38]]
NEUTRAL = np.array(_NH3 + [[0.0, 0.0, 1.75], [0.0, 0.0, 2.70]])  # N···H–F
ION_PAIR = np.array(_NH3 + [[0.0, 0.0, 1.05], [0.0, 0.0, 2.70]])  # N–H···F⁻


def _line(symbols: str, r: float) -> tuple[list[str], np.ndarray]:
    atoms = symbols.split()
    return atoms, np.array([[k * r, 0.0, 0.0] for k in range(len(atoms))])


@pytest.mark.parametrize("symbols,r,bonded", [
    ("F H F", 1.14, True),  # FHF⁻: one fragment
    ("I I I", 2.92, True),  # I3⁻
    ("N I", 2.8, False),  # halogen bond
    ("I I", 3.5, False),
])
def test_additive_bond_rule(symbols, r, bonded):
    atoms, x = _line(symbols, r)
    assert (len(top.fragments(atoms, x)) == 1) is bonded


def test_halogen_bonded_complex_is_two_fragments():
    x = np.array(_NH3 + [[0.0, 0.0, 2.8], [0.0, 0.0, 2.8 + 2.32]])  # N···I 2.8 Å, I–Cl 2.32 Å
    assert top.fragments(["N", "H", "H", "H", "I", "Cl"], x) == ((0, 1, 2, 3), (4, 5))


def test_bond_changes_are_symmetric_and_within_the_bond_graph_differences():
    rng = np.random.default_rng(7)
    for _ in range(50):
        a, b = (NEUTRAL + rng.normal(scale=0.4, size=NEUTRAL.shape) for _ in range(2))
        formed, broken = top.bond_changes(NH3_HF, a, b)
        assert top.bond_changes(NH3_HF, b, a) == (broken, formed)
        assert formed <= top.bonds(NH3_HF, b) - top.bonds(NH3_HF, a)
        assert broken <= top.bonds(NH3_HF, a) - top.bonds(NH3_HF, b)


def test_proton_transfer_changes():
    formed, broken = top.bond_changes(NH3_HF, NEUTRAL, ION_PAIR)
    assert formed == {(0, 4)} and broken == {(4, 5)}
    assert top.bond_changes(NH3_HF, NEUTRAL, NEUTRAL) == (frozenset(), frozenset())


def _nh3_hf(nh: float) -> np.ndarray:
    """NH3·HF with N···H at r_thr + nh (Å) and H–F 0.95 Å."""
    r = covalent_radius("N") + covalent_radius("H") + top.BOND_TOLERANCE_A + nh
    return np.array(_NH3 + [[0.0, 0.0, r], [0.0, 0.0, r + 0.95]])


def test_bond_changes_clear_the_band_on_both_sides():
    """One band for every bond change (analysis X3, G8-P5): a threshold crossing inside one basin
    (S19's N···H, ±0.006 Å about r_thr) splits the bond graph and the state label but is no
    change; S6's H–O (+2.00 / −0.35 Å) is one; one side inside ±RESOLVED_A is none."""
    noise = _nh3_hf(+0.006), _nh3_hf(-0.006)
    assert top.bonds(NH3_HF, noise[1]) - top.bonds(NH3_HF, noise[0]) == {(0, 4)}
    assert top.state_label(NH3_HF, noise[0]) != top.state_label(NH3_HF, noise[1])
    assert top.bond_changes(NH3_HF, *noise) == (set(), set())
    assert top.bond_changes(NH3_HF, _nh3_hf(+2.0), _nh3_hf(-0.35)) == ({(0, 4)}, set())
    assert top.bond_changes(NH3_HF, _nh3_hf(-0.35), _nh3_hf(+2.0)) == (set(), {(0, 4)})
    for inside in (top.RESOLVED_A - 0.01, -top.RESOLVED_A + 0.01):
        assert top.bond_changes(NH3_HF, _nh3_hf(inside), _nh3_hf(-0.35)) == (set(), set())
        assert top.bond_changes(NH3_HF, _nh3_hf(+2.0), _nh3_hf(inside)) == (set(), set())
    outside = _nh3_hf(top.RESOLVED_A + 0.01), _nh3_hf(-top.RESOLVED_A - 0.01)
    assert top.bond_changes(NH3_HF, *outside) == ({(0, 4)}, set())


def test_state_label_is_permutation_invariant():
    order = [5, 3, 0, 4, 1, 2]
    permuted = top.state_label([NH3_HF[i] for i in order], NEUTRAL[order])
    assert permuted == top.state_label(NH3_HF, NEUTRAL)
    assert top.state_label(NH3_HF, ION_PAIR) != top.state_label(NH3_HF, NEUTRAL)
    assert top.state_label(*_line("F H F", 1.14)).startswith("F2H_")


def test_wl_classes_are_atom_equivalence_and_permutation_invariant():
    bonded = top.bonds(NH3_HF, NEUTRAL)
    classes = top.wl_classes(NH3_HF, bonded)
    assert classes[1] == classes[2] == classes[3] and len(set(classes)) == 4  # N, 3 H, H, F
    order = [5, 3, 0, 4, 1, 2]
    permuted = top.wl_classes([NH3_HF[i] for i in order], top.bonds(
        [NH3_HF[i] for i in order], NEUTRAL[order]))
    assert permuted == tuple(classes[i] for i in order)

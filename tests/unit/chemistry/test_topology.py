"""Additive bond rule, fragments, bond changes resolved by one band (analysis X3), canonical atom
classes and state labels with stereo (design §5.5, review U1-P2/U1-P3/U1-P4). RDKit: WSL only."""

import hashlib
import re

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from hfauto.chemistry import topology as top
from hfauto.chemistry.elements import covalent_radius
from hfauto.chemistry.smiles import smiles_to_xyz

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
    assert top.same_bonding(NH3_HF, *noise)
    assert top.bond_changes(NH3_HF, _nh3_hf(+2.0), _nh3_hf(-0.35)) == ({(0, 4)}, set())
    assert not top.same_bonding(NH3_HF, _nh3_hf(+2.0), _nh3_hf(-0.35))
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
    assert re.fullmatch(r"F2H_[0-9a-f]{16}", top.state_label(*_line("F H F", 1.14)))


@pytest.mark.parametrize("fused,linked", [
    ("C1CCC2CCCCC2C1", "C1CCC(C1)C1CCCC1"),  # decalin / bicyclopentyl
    ("C1CC2CCCC2C1", "C1CC(C1)C1CCC1"),  # bicyclo[3.3.0]octane / bicyclobutyl
    ("C1CC2CCC12", "C1CC1C1CC1"),  # bicyclo[2.2.0]hexane / bicyclopropyl
])
def test_state_labels_separate_isomers_that_colour_refinement_cannot(fused, linked):
    """U1-I3: two CH bonded to each other and n CH2, two rings either way; the 3-round WL hash
    gave both one label."""
    labels = set()
    for smiles in (fused, linked):
        x = smiles_to_xyz(smiles, 0)
        assert len(top.bonds(x.symbols, x.coords)) == len(x.symbols) + 1  # two rings
        labels.add(top.state_label(x.symbols, x.coords))
    assert len(labels) == 2 and len({label.split("_")[0] for label in labels}) == 1


def test_labels_and_classes_are_invariant_under_100_atom_permutations():
    x = smiles_to_xyz("C1CCC2CCCCC2C1", 0)  # decalin: classes with ties
    label, classes = top.state_label(x.symbols, x.coords), top.atom_classes(
        x.symbols, top.bonds(x.symbols, x.coords))
    rng = np.random.default_rng(11)
    for _ in range(100):
        order = rng.permutation(len(x.symbols))
        symbols, coords = [x.symbols[i] for i in order], x.coords[order]
        assert top.state_label(symbols, coords) == label
        assert top.atom_classes(symbols, top.bonds(symbols, coords)) == tuple(
            classes[i] for i in order)


def _turned(x: np.ndarray, symbols: list[str], i: int, j: int) -> np.ndarray:
    """x with the side of bond i-j that holds j turned by 180° about it."""
    bonded = top.bonds(symbols, x) - {(min(i, j), max(i, j))}
    side = next(g for g in top._components(len(symbols), bonded) if j in g)
    u = (x[j] - x[i]) / np.linalg.norm(x[j] - x[i])
    y = x.copy()
    y[list(side)] = Rotation.from_rotvec(np.pi * u).apply(x[list(side)] - x[j]) + x[j]
    return y


def _label(smiles: str) -> str:
    x = smiles_to_xyz(smiles, 0)
    return top.state_label(x.symbols, x.coords)


def test_double_bond_stereo_separates_states_and_single_bonds_carry_none():
    """U1-P4: E- and Z-2-butene are two states; s-cis/s-trans butadiene and the two amide
    rotamers of N-methylformamide are one (their C-C and C-N bonds are single in the Lewis
    structure); a charge-separated double bond (C-=O+ of hydroxymethylene) carries none."""
    assert _label("C/C=C/C") != _label("C/C=C\\C")
    for smiles, (i, j) in (("C=CC=C", (1, 2)), ("CNC=O", (1, 2)), ("[CH-]=[OH+]", (0, 1))):
        x = smiles_to_xyz(smiles, 0)
        turned = _turned(x.coords, x.symbols, i, j)
        assert top.state_label(x.symbols, turned) == top.state_label(x.symbols, x.coords)


def test_mirror_images_are_one_state_and_diastereomers_two():
    """The lower stereo text of a structure and its mirror image counts; one centre alone tells
    nothing and keeps the bond graph's label (a flattened centre splits nothing); a radical
    keeps its tetrahedral centres."""
    rr, ss, meso = (_label(s) for s in ("C[C@@H](F)[C@@H](F)C", "C[C@H](F)[C@H](F)C",
                                        "C[C@@H](F)[C@H](F)C"))
    assert rr == ss != meso
    x = smiles_to_xyz("C[C@H](F)Cl", 0)
    graph = top._canonical(x.symbols, top.bonds(x.symbols, x.coords))[0]
    plain = f"C2H4ClF_{hashlib.sha256(graph.encode()).hexdigest()[:16]}"
    assert top.state_label(x.symbols, x.coords) == top.state_label(x.symbols, -x.coords) == plain
    assert _label("C[C@H](F)[C@H](F)[CH2]") != _label("C[C@H](F)[C@@H](F)[CH2]")


@pytest.mark.parametrize("smiles", ["C/C=C/C", "C[C@H](F)[C@H](F)C", "C[CH-][NH+]=C",
                                    "C1CCC2CCCCC2C1"])
def test_stereo_labels_are_invariant_under_atom_permutations(smiles):
    """The Lewis structure is perceived in canonical atom order: a resonance form (the ylide's)
    does not depend on the input order."""
    x = smiles_to_xyz(smiles, 0)
    label, rng = top.state_label(x.symbols, x.coords), np.random.default_rng(13)
    for _ in range(30):
        order = rng.permutation(len(x.symbols))
        assert top.state_label([x.symbols[i] for i in order], x.coords[order]) == label


def test_atom_classes_are_the_refined_equivalence():
    classes = top.atom_classes(NH3_HF, top.bonds(NH3_HF, NEUTRAL))
    assert classes[1] == classes[2] == classes[3] and len(set(classes)) == 4  # N, 3 H, H, F
    x = smiles_to_xyz("CCCCCCCCCl", 0)  # 3 WL rounds left two of the 9 heavy atoms in one class
    classes = top.atom_classes(x.symbols, top.bonds(x.symbols, x.coords))
    assert len({c for c, s in zip(classes, x.symbols, strict=True) if s != "H"}) == 9

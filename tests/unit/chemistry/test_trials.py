"""Reaction trials (design §8.2 explore): one bond-graph edit enumerator. Mechanisms come out of
it as results, its counts equal the P0c/P0f probe's, its classes do not depend on numbering or
orientation, and each class is realised once per state (conformer, placement, linear bend)."""

from collections import Counter

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from hfauto.chemistry import topology, trials
from hfauto.chemistry.smiles import smiles_to_xyz
from hfauto.chemistry.vibrations import external_basis

Chem = pytest.importorskip("rdkit.Chem")

SN2 = (["C", "H", "H", "H", "Cl", "Cl"],  # Cl⁻ ··· CH3Cl ion-dipole complex (backside)
       np.array([[0, 0, 0], [1.034, 0, -0.346], [-0.517, 0.8955, -0.346],
                 [-0.517, -0.8955, -0.346], [0, 0, 1.8], [0, 0, -3.2]]), -1, 1)
CH3O2 = (["C", "H", "H", "H", "O", "O"],  # VAL9 R2 s5 source
         np.array([[0.9538, 0.4127, 0.2711], [0.8582, 0.2654, -0.8029], [0.6098, -0.4718, 0.8165],
                   [0.3573, 1.2698, 0.5996], [2.6451, 0.8372, 1.7518], [2.3262, 0.6409, 0.5042]]),
         0, 2)
CH3_H2O = (["C", "H", "H", "H", "H", "O", "H"],  # VAL9 R2 s6 source: CH3 ··· H2–O5–H6
           np.array([[0.1859, 0.0949, -0.0402], [0.7868, -0.4360, 0.6731],
                     [1.3051, 0.8140, -0.8816], [-0.3630, 0.9447, 0.3178],
                     [-0.2675, -0.4837, -0.8211], [2.0548, 1.3095, -1.3642],
                     [2.5480, 1.7931, -0.6988]]), 0, 2)
NH3_HF = (["N", "H", "H", "H", "H", "F"],  # VAL9 R2 s10 source: H3N ··· H4–F5
          np.array([[-1.3139, 0, 0], [-1.6621, 0.1954, -0.9302], [-1.6617, -0.9035, 0.2958],
                    [-1.6622, 0.7077, 0.6345], [0.2409, 0.0007, 0.0004], [1.2204, 0, 0]]), 0, 1)
BH3_NH3 = (["B", "H", "H", "H", "N", "H", "H", "H"],  # configs bh3_nh3 complex: B···N 3.2 Å
           np.array([[0, 0, 0], [1.19, 0, 0], [-0.595, 1.03057, 0], [-0.595, -1.03057, 0],
                     [0, 0, 3.2], [0.470154, 0.814331, 3.582157], [-0.940309, 0, 3.582157],
                     [0.470154, -0.814331, 3.582157]]), 0, 1)
H3N_H_F = (["N", "H", "H", "H", "H", "F"],  # H4 bonded to N0 (1.30 Å) and F5 (1.00 Å)
           np.vstack([NH3_HF[1][:4], [[-0.0139, 0, 0], [0.9861, 0, 0]]]), 0, 1)
HCN = (["H", "C", "N"], np.array([[0, 0, -1.066], [0, 0, 0], [0, 0, 1.156]]), 0, 1)


def _smiles(*smiles, gap=6.0):
    """ETKDG structures, each next one 6 Å further along z: separate fragments."""
    parts = [smiles_to_xyz(s, 0) for s in smiles]
    return ([e for p in parts for e in p.symbols],
            np.vstack([p.coords - p.coords.mean(axis=0) + [0, 0, gap * k]
                       for k, p in enumerate(parts)]), 0, 1)


def _state(symbols, bonds):
    """Canonical SMILES of a bond graph: a state, whatever the atom numbering."""
    mol = Chem.RWMol(Chem.MolFromSmiles(".".join(f"[{s}]" for s in symbols), sanitize=False))
    for i, j in bonds:
        mol.AddBond(int(i), int(j), Chem.BondType.SINGLE)
    return Chem.MolToSmiles(mol)


def _smiles_state(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    return _state([a.GetSymbol() for a in mol.GetAtoms()],
                  [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()])


def _edit(formed, broken):
    return (frozenset(tuple(sorted(p)) for p in formed),
            frozenset(tuple(sorted(p)) for p in broken))


def _product(system, formed, broken):
    symbols, x, *_ = system
    f, b = _edit(formed, broken)
    return _state(symbols, (topology.bonds(symbols, x) - b) | f)


def _found(system):
    """(class -> members, the product states of the classes, every member edit)."""
    symbols, x, charge, multiplicity = system
    classes, bonds = trials.edits(symbols, x, charge, multiplicity), topology.bonds(symbols, x)
    products = {_state(symbols, (bonds - b) | f) for f, b in (m[0] for m in classes.values())}
    return classes, products, {e for m in classes.values() for e in m}


FLUOROETHANE, HEXADIENE, ETHYL_FORMATE = _smiles("CCF"), _smiles("C=CCCC=C"), _smiles("O=COCC")
BUTADIENE_ETHENE, NME3_SO2 = _smiles("C=CC=C", "C=C"), _smiles("CN(C)C", "O=S=O")


def test_mechanisms_are_results_of_one_enumeration():
    """P0c (i) with no per-mechanism code: the 4-centre 1,2-HF elimination, the [3,3] shift,
    the [4+2] cycloaddition, the retro-ene and the plain 1,5-H shift (a 1,5-zwitterion, kept
    because P0f left charge separation unlimited)."""
    assert _smiles_state("C=C.F") in _found(FLUOROETHANE)[1]
    assert _edit([(0, 5)], [(2, 3)]) in _found(HEXADIENE)[2]  # degenerate: the class itself
    assert _smiles_state("C1=CCCCC1") in _found(BUTADIENE_ETHENE)[1]
    _, products, members = _found(ETHYL_FORMATE)  # O0=C1–O2–C3–C4, H8 on C4
    assert _smiles_state("OC=O.C=C") in products
    assert _edit([(0, 8)], [(4, 8)]) in members


def test_ion_pairs_adducts_and_sn2():
    assert _product(NH3_HF, [(0, 4)], [(4, 5)]) in _found(NH3_HF)[1]  # NH4+ F−
    assert _smiles_state("[NH3+][BH3-]") in _found(BH3_NH3)[1]
    assert _smiles_state("C[N+](C)(C)[S-](=O)=O") in _found(NME3_SO2)[1]
    assert _edit([(0, 5)], [(0, 4)]) in _found(SN2)[2]  # backside Cl' in, Cl out


def test_every_r2_product_state_is_a_class_product():
    """VAL9 R2 products at state level: s5 CH2OOH, s6 CH3OH + H (eddde71's atom-mapped f2b2
    edit swaps a spectator H and is not local; its state comes from the f1b1 class C0–O5 formed,
    H2–O5 broken) and the s10 degenerate double H exchange (a ring)."""
    assert _product(CH3O2, [(2, 4)], [(0, 2)]) in _found(CH3O2)[1]
    _, products, members = _found(CH3_H2O)
    eddde71 = _product(CH3_H2O, [(0, 2), (0, 5)], [(0, 4), (2, 5)])
    assert eddde71 == _product(CH3_H2O, [(0, 5)], [(5, 6)]) in products  # 316385
    assert _edit([(0, 2), (0, 5)], [(0, 4), (2, 5)]) not in members
    assert _edit([(0, 5)], [(2, 5)]) in members
    assert _edit([(0, 4), (2, 5)], [(0, 2), (4, 5)]) in _found(NH3_HF)[2]


@pytest.mark.parametrize(("system", "by_type"), [
    (FLUOROETHANE, {"f1b2": 8, "f2b2": 2}),
    (HEXADIENE, {"f1b0": 3, "f1b1": 25, "f1b2": 71, "f2b0": 2, "f2b1": 41, "f2b2": 297}),
    (BUTADIENE_ETHENE, {"f1b0": 4, "f1b1": 19, "f1b2": 33, "f2b0": 5, "f2b1": 42, "f2b2": 199}),
    (ETHYL_FORMATE, {"f1b0": 1, "f1b1": 21, "f1b2": 53, "f2b1": 28, "f2b2": 158}),
    (NME3_SO2, {"f1b0": 3, "f1b1": 11, "f1b2": 13, "f2b1": 25, "f2b2": 68}),
    (SN2, {"f1b1": 4, "f1b2": 2}), (CH3O2, {"f1b1": 4, "f1b2": 8, "f2b2": 6}),
    (CH3_H2O, {"f1b1": 4, "f1b2": 6, "f2b2": 2}), (NH3_HF, {"f1b1": 2, "f1b2": 3, "f2b2": 2}),
    (BH3_NH3, {"f1b0": 1, "f1b1": 2, "f1b2": 4, "f2b2": 2})])
def test_class_counts_equal_the_probe(system, by_type):
    """The P0c counts after the Lewis filter (P0f kept P0c's charge rule), with no cap."""
    classes = _found(system)[0]
    assert Counter(f"f{len(m[0][0])}b{len(m[0][1])}" for m in classes.values()) == by_type


@pytest.mark.parametrize("system", [ETHYL_FORMATE, SN2, CH3_H2O, NH3_HF])
def test_classes_do_not_depend_on_numbering_or_orientation(system):
    symbols, x, charge, multiplicity = system
    keys = set(trials.edits(symbols, x, charge, multiplicity))
    rng = np.random.default_rng(20261005)
    for _ in range(20):
        order = rng.permutation(len(symbols))
        assert set(trials.edits([symbols[i] for i in order], x[order], charge,
                                multiplicity)) == keys
    turned = Rotation.random(random_state=rng).apply(x) + 3.0
    assert set(trials.edits(symbols, turned, charge, multiplicity)) == keys


@pytest.mark.parametrize(("symbols", "bonds", "charge", "multiplicity", "expected"), [
    (["C", "N", "O", "O", "H", "H", "H"], [(0, 1), (1, 2), (1, 3), (0, 4), (0, 5), (0, 6)],
     0, 1, True),  # nitromethane: N+ and O−
    (["O", "S", "O"], [(0, 1), (1, 2)], 0, 1, True),  # S valence 4
    (["N", "C", "B", *"HHHHHHHH"], [(0, 1), (1, 2), (0, 3), (0, 4), (0, 5), (1, 6), (1, 7),
                                    (2, 8), (2, 9), (2, 10)], 0, 1, True),  # 1,3-zwitterion
    (["O", "H"], [(0, 1)], 0, 2, True), (["O", "O"], [(0, 1)], 0, 3, True),
    (["O", "H"], [(0, 1)], 0, 1, False),  # a singlet OH has no Lewis structure
    (["Fe", "Cl", "Cl", "Cl"], [(0, 1), (0, 2), (0, 3)], 0, 6, False),  # d block: out of scope
    (H3N_H_F[0], [(0, 1), (0, 2), (0, 3), (0, 4), (4, 5)], 0, 1, False)])  # H with two bonds
def test_lewis_structures(symbols, bonds, charge, multiplicity, expected):
    assert trials.lewis(symbols, frozenset(bonds), charge, multiplicity) is expected


def test_a_source_without_a_lewis_structure_is_still_enumerated():
    """Only products are filtered (VAL9 R2 s19: an H bonded to N and F)."""
    symbols, x, *_ = H3N_H_F
    assert not trials.lewis(symbols, topology.bonds(symbols, x), 0, 1)
    assert _product(H3N_H_F, [(1, 5)], [(0, 1), (4, 5)]) in _found(H3N_H_F)[1]


def test_each_class_runs_once_on_its_closest_conformer_in_budget_order():
    symbols, x, charge, multiplicity = NH3_HF
    far = x.copy()
    far[4:] += [0.3, 0, 0]  # the same state, H4···N0 0.3 Å longer
    ts = trials.trials(symbols, [far, x], charge, multiplicity)
    assert sorted(t.key for t in ts) == sorted(trials.edits(*NH3_HF))
    assert [t.order for t in ts] == sorted(t.order for t in ts)
    assert [len(t.formed) + len(t.broken) for t in ts] == [2, 2, 3, 3, 3, 4, 4]
    [ion_pair] = [t for t in ts if (t.formed, t.broken) == (((0, 4),), ((4, 5),))]
    assert ion_pair.conformer == 1 and np.array_equal(ion_pair.start, x)
    assert ion_pair.drive == pytest.approx(np.linalg.norm(x[0] - x[4]) / (0.71 + 0.31))


def test_an_intermolecular_class_without_contact_starts_placed():
    """B···N 6 Å: N's fragment turns rigidly to face B and moves to 1.5 Σr_cov; at 3.2 Å the
    pair is in contact (Σr_vdW 3.47 Å) and the conformer is the start."""
    symbols, x, charge, multiplicity = BH3_NH3
    apart = x.copy()
    apart[4:] += [0, 0, 2.8]

    def start(conformer):
        [adduct] = [t for t in trials.trials(symbols, [conformer], charge, multiplicity)
                    if (t.formed, t.broken) == (((0, 4),), ())]
        return adduct.start

    assert np.array_equal(start(x), x)
    s = start(apart)
    assert np.array_equal(s[:4], x[:4])  # BH3 stays: equal sizes, the pair's second atom moves
    assert np.linalg.norm(s[4] - s[0]) == pytest.approx(1.5 * (0.84 + 0.71))
    outward = s[4] - s[4:].mean(axis=0)  # N's lone-pair side faces B
    assert np.cross(outward, s[0] - s[4]) == pytest.approx(np.zeros(3), abs=1e-9)
    assert outward @ (s[0] - s[4]) > 0
    assert np.allclose(np.linalg.norm(s[4:, None] - s[None, 4:], axis=-1),
                       np.linalg.norm(x[4:, None] - x[None, 4:], axis=-1))  # rigid


def test_a_linear_start_is_bent_by_10_degrees_deterministically():
    symbols, x, charge, multiplicity = HCN
    [shift] = trials.trials(symbols, [x], charge, multiplicity)  # H 1,2-shift C -> N
    assert (shift.formed, shift.broken) == (((0, 2),), ((0, 1),))
    v1, v2 = shift.start[0] - shift.start[1], shift.start[2] - shift.start[1]
    bend = 180 - np.degrees(np.arccos(v1 @ v2 / np.linalg.norm(v1) / np.linalg.norm(v2)))
    assert bend == pytest.approx(10.0) and external_basis(symbols, shift.start).shape[1] == 6
    assert np.array_equal(shift.start, trials.trials(symbols, [x], charge, multiplicity)[0].start)

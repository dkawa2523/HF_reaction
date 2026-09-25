"""Geometry fingerprints, Molecule and composition keys (design §5.5)."""

import numpy as np

from hfauto.chemistry.xyz import XYZ, Molecule, composition_key, geometry_fingerprint, read_xyz

WATER = XYZ(
    ["O", "H", "H"], np.array([[0.0, 0.0, 0.1234561], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]])
)


def test_fingerprint_rounds_to_micro_angstrom():
    base = geometry_fingerprint(WATER.symbols, WATER.coords)
    assert geometry_fingerprint(WATER.symbols, WATER.coords + 1e-9) == base
    assert geometry_fingerprint(WATER.symbols, WATER.coords + 1e-5) != base
    assert geometry_fingerprint(["S", "H", "H"], WATER.coords) != base
    origin = geometry_fingerprint(["H"], np.zeros((1, 3)))
    assert geometry_fingerprint(["H"], np.array([[-0.0, 0.0, 0.0]])) == origin


def test_molecule_fingerprint_includes_charge_and_multiplicity(tmp_path):
    neutral = Molecule(WATER, 0, 1)
    assert neutral.fingerprint() != Molecule(WATER, 1, 2).fingerprint()
    assert neutral.fingerprint() != Molecule(WATER, 0, 3).fingerprint()
    back = read_xyz(neutral.write(tmp_path / "water.xyz"))
    assert back.symbols == WATER.symbols and np.allclose(back.coords, WATER.coords)


def test_composition_key_hill_order():
    assert composition_key(["H", "C", "N"], 0, 1) == "CHN_q0_m1"
    assert composition_key(["N", "H", "H", "H", "H", "F"], 0, 1) == "FH4N_q0_m1"
    assert composition_key(["F", "H", "F"], -1, 1) == "F2H_q-1_m1"

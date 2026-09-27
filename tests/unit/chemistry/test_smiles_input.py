"""SMILES input: one ETKDG structure; multi-fragment SMILES are compositions (CH-23)."""

import pytest

from hfauto.chemistry.smiles import smiles_to_molecule


def test_smiles_embeds_one_structure_and_rejects_fragments():
    with pytest.raises(ValueError, match="composition"):
        smiles_to_molecule("N.F", 0, 1)  # rejected before RDKit is needed
    pytest.importorskip("rdkit")
    mol = smiles_to_molecule("C[NH3+]", 1, 1)
    assert sorted(mol.xyz.symbols) == ["C", "H", "H", "H", "H", "H", "H", "N"]
    assert (mol.charge, mol.multiplicity, mol.xyz.coords.shape) == (1, 1, (8, 3))
    with pytest.raises(ValueError, match="formal charge"):
        smiles_to_molecule("C[NH3+]", 0, 1)
    with pytest.raises(ValueError, match="isotopes"):
        smiles_to_molecule("[2H]O[2H]", 0, 1)  # D2O would silently get 1H masses


@pytest.mark.parametrize(("smiles", "charge", "declared", "expected"), [
    ("[O][O]", 0, None, 3), ("[CH3]", 0, None, 2), ("C", 0, None, 1),
    ("[Fe+2]", 2, None, 1),  # no radicals: structures then asks for a declaration
    ("[O][O]", 0, 1, 1),  # a declared value is kept (parity is checked by structures)
])
def test_undeclared_multiplicity_is_radical_electrons_plus_one(smiles, charge, declared,
                                                               expected):
    pytest.importorskip("rdkit")
    assert smiles_to_molecule(smiles, charge, declared).multiplicity == expected

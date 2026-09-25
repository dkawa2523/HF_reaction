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

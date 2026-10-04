"""SMILES input: one ETKDG structure; multi-fragment SMILES are compositions (CH-23)."""

import pytest

from hfauto.chemistry.smiles import smiles_to_xyz


def test_smiles_embeds_one_structure_and_rejects_fragments():
    with pytest.raises(ValueError, match="composition"):
        smiles_to_xyz("N.F", 0)  # rejected before RDKit is needed
    pytest.importorskip("rdkit")
    xyz = smiles_to_xyz("C[NH3+]", 1)
    assert sorted(xyz.symbols) == ["C", "H", "H", "H", "H", "H", "H", "N"]
    assert xyz.coords.shape == (8, 3)
    with pytest.raises(ValueError, match="formal charge"):
        smiles_to_xyz("C[NH3+]", 0)
    with pytest.raises(ValueError, match="isotopes"):
        smiles_to_xyz("[2H]O[2H]", 0)  # D2O would silently get 1H masses


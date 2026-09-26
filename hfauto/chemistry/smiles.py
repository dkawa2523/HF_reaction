"""SMILES to one 3D structure with RDKit ETKDG (design §8.2 structures; RDKit is optional)."""

from __future__ import annotations

import numpy as np

from hfauto.chemistry.xyz import XYZ, Molecule


def smiles_to_molecule(smiles: str, charge: int, multiplicity: int, *,
                       random_seed: int = 20260925) -> Molecule:
    """One ETKDG embedding with explicit hydrogens.

    A multi-fragment SMILES ('.') is rejected: complexes are compositions (CH-23). The SMILES
    formal charge must equal ``charge`` and isotopes are rejected. ImportError without RDKit.
    """
    if "." in smiles:
        raise ValueError(f"multi-fragment SMILES {smiles!r}: declare a composition instead")
    from rdkit import Chem
    from rdkit.Chem import rdDistGeom  # not AllChem: it loads far more native modules

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES {smiles!r}")
    formal = Chem.GetFormalCharge(mol)
    if formal != charge:
        raise ValueError(f"SMILES {smiles!r} has formal charge {formal}, declared {charge}")
    if any(a.GetIsotope() for a in mol.GetAtoms()):
        raise ValueError(f"SMILES {smiles!r} specifies isotopes; "
                         "hfauto uses most-abundant-isotope masses")
    mol = Chem.AddHs(mol)
    params = rdDistGeom.ETKDGv3()
    params.randomSeed = random_seed
    if rdDistGeom.EmbedMolecule(mol, params) != 0:
        raise ValueError(f"ETKDG could not embed {smiles!r}")
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    coords = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    return Molecule(XYZ(symbols, coords, comment=f"smiles={smiles}"), charge, multiplicity)

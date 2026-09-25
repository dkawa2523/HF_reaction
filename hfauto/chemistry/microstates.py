"""Bounded molecular-state enumeration for reaction discovery.

The enumerator deliberately makes only one local change at a time.  It is not a
pH/speciation predictor: it supplies chemically valid structures that can be
screened by the downstream energy workflow.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class MolecularState:
    smiles: str
    kind: str
    source_atom_index: int | None = None


def _canonical(mol) -> str:
    from rdkit import Chem

    return Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True, isomericSmiles=True)


def _feature_atoms(mol, family: str) -> set[int]:
    """Return RDKit chemical-feature atom ids without exposing toolkit objects."""

    from pathlib import Path

    from rdkit import Chem, RDConfig
    from rdkit.Chem import ChemicalFeatures

    factory = ChemicalFeatures.BuildFeatureFactory(
        str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef")
    )
    explicit = Chem.AddHs(mol)
    atoms = {
        int(atom_index)
        for feature in factory.GetFeaturesForMol(explicit)
        if feature.GetFamily() == family
        for atom_index in feature.GetAtomIds()
        if atom_index < mol.GetNumAtoms()
    }
    if family == "Acceptor":
        # BaseFeatures intentionally omits small inorganic bases such as NH3
        # and H2O.  The fallback is conservative: neutral/negative p-block
        # heteroatoms are admitted, while cations and aromatic [nH] are not.
        atoms.update(
            atom.GetIdx()
            for atom in explicit.GetAtoms()
            if atom.GetIdx() < mol.GetNumAtoms()
            and atom.GetAtomicNum() in {7, 8, 9, 15, 16, 17, 35, 53}
            and atom.GetFormalCharge() <= 0
            and not (
                atom.GetIsAromatic()
                and any(neighbor.GetAtomicNum() == 1 for neighbor in atom.GetNeighbors())
            )
        )
    return atoms


def _protonate(mol, atom_index: int):
    from rdkit import Chem

    editable = Chem.RWMol(Chem.AddHs(mol))
    atom = editable.GetAtomWithIdx(int(atom_index))
    hydrogen_index = editable.AddAtom(Chem.Atom(1))
    editable.AddBond(int(atom_index), hydrogen_index, Chem.BondType.SINGLE)
    atom.SetFormalCharge(atom.GetFormalCharge() + 1)
    candidate = editable.GetMol()
    Chem.SanitizeMol(candidate)
    return Chem.RemoveHs(candidate)


def _deprotonate(mol, atom_index: int):
    from rdkit import Chem

    editable = Chem.RWMol(Chem.AddHs(mol))
    atom = editable.GetAtomWithIdx(int(atom_index))
    hydrogen = next(
        (neighbor.GetIdx() for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 1),
        None,
    )
    if hydrogen is None:
        raise ValueError("selected donor has no hydrogen")
    editable.RemoveAtom(int(hydrogen))
    # Hydrogens are appended by AddHs, so removing one does not renumber the
    # selected heavy atom.
    atom = editable.GetAtomWithIdx(int(atom_index))
    atom.SetFormalCharge(atom.GetFormalCharge() - 1)
    candidate = editable.GetMol()
    Chem.SanitizeMol(candidate)
    return Chem.RemoveHs(candidate)


def enumerate_molecular_states(
    smiles: str,
    *,
    include_tautomers: bool = True,
    include_protomers: bool = True,
    max_tautomers: int = 8,
    max_protomers: int = 8,
) -> list[MolecularState]:
    """Enumerate the input state plus bounded one-step tautomers/protomers."""

    from rdkit import Chem

    source = Chem.MolFromSmiles(smiles)
    if source is None:
        raise ValueError(f"invalid SMILES: {smiles!r}")
    states: list[MolecularState] = []
    seen: set[str] = set()

    def add(mol, kind: str, atom_index: int | None = None) -> None:
        value = _canonical(mol)
        if value in seen:
            return
        seen.add(value)
        states.append(MolecularState(value, kind, atom_index))

    add(source, "input")
    if include_tautomers:
        from rdkit.Chem.MolStandardize import rdMolStandardize

        enumerator = rdMolStandardize.TautomerEnumerator()
        for tautomer in list(enumerator.Enumerate(source))[: max(0, int(max_tautomers))]:
            add(tautomer, "tautomer")

    if include_protomers:
        acceptors = sorted(_feature_atoms(source, "Acceptor"))
        donors = sorted(_feature_atoms(source, "Donor"))
        # RDKit's generic donor rules do not include hydrogen halides.  An
        # electronegative atom bearing H is nevertheless a valid proton donor.
        explicit = Chem.AddHs(source)
        donors.extend(
            atom.GetIdx()
            for atom in explicit.GetAtoms()
            if atom.GetIdx() < source.GetNumAtoms()
            and atom.GetAtomicNum() in {7, 8, 9, 15, 16, 17, 35, 53}
            and any(neighbor.GetAtomicNum() == 1 for neighbor in atom.GetNeighbors())
        )
        for atom_index in acceptors[: max(0, int(max_protomers))]:
            try:
                add(_protonate(source, atom_index), "protonated", atom_index)
            except (RuntimeError, ValueError):
                continue
        remaining = max(0, int(max_protomers) - len(acceptors))
        for atom_index in sorted(set(donors))[:remaining]:
            try:
                add(_deprotonate(source, atom_index), "deprotonated", atom_index)
            except (RuntimeError, ValueError):
                continue
    return states

"""Covalent bond graph, fragments, bond changes, canonical atom classes and state labels (§5.5).

A pair is bonded when r < r_thr = r_cov,i + r_cov,j + 0.4 Å (additive tolerance, Meng & Lewis
1991; SCINE BondDetector, OpenBabel): FHF⁻ (1.14 Å) and I3⁻ (2.92 Å) are bonded, a halogen bond
I···N at 2.8 Å is not. A bond change between two structures clears the band ±RESOLVED_A around
r_thr on both sides, so ``bond_changes`` is symmetric and a contact the threshold cuts within one
basin is no change. Two structures can then differ in ``bonds`` and state label with no bond
change: a split label only makes another state, never a merge. Graph identity is RDKit's
canonicalisation (Schneider, Sayle & Landrum 2015), whose version the environment fixes.

A state label also carries stereo (review U1-P4): each fragment's bond orders come from
``DetermineBondOrders`` at charge 0 in canonical atom order, so they depend on the graph alone;
where it finds a Lewis structure, E/Z of its double bonds between uncharged atoms and
tetrahedral stereo are read from 3D, elsewhere (radicals, ions) tetrahedral only. A structure
and its mirror image are one state: the lower of the two stereo texts counts, and one
tetrahedral centre alone counts as none. A structure without stereo keeps the label of its bond
graph.
"""

from __future__ import annotations

import hashlib
from collections.abc import Collection, Sequence
from typing import TYPE_CHECKING

import numpy as np

from hfauto.chemistry.elements import covalent_radius
from hfauto.chemistry.xyz import hill_formula

if TYPE_CHECKING:
    from rdkit.Chem import Mol

Bond = tuple[int, int]  # (i, j) with i < j
BOND_TOLERANCE_A = 0.4
# Half-width of the band around r_thr that resolves a bond change: GFN2 and PBE0 N···H of the
# TMA·(HF)2 amine·HF basin differ by 0.15 Å, so a crossing within 0.1 Å of r_thr is noise.
RESOLVED_A = 0.1


def _excess(symbols: Sequence[str], coords: np.ndarray) -> np.ndarray:
    """r − r_thr (Å) per atom pair: negative where bonded."""

    x = np.asarray(coords, dtype=float).reshape(-1, 3)
    if len(x) != len(symbols):
        raise ValueError("symbols and coordinates differ in atom count")
    radii = np.array([covalent_radius(s) for s in symbols])
    r = np.linalg.norm(x[:, None] - x[None], axis=-1)
    return r - (np.add.outer(radii, radii) + BOND_TOLERANCE_A)


def _pairs(mask: np.ndarray) -> frozenset[Bond]:
    i, j = np.nonzero(np.triu(mask, 1))
    return frozenset(zip(i.tolist(), j.tolist(), strict=True))


def bonds(symbols: Sequence[str], coords: np.ndarray) -> frozenset[Bond]:
    return _pairs(_excess(symbols, coords) < 0)


def _components(n_atoms: int, bonded: Collection[Bond]) -> tuple[tuple[int, ...], ...]:
    parent = list(range(n_atoms))

    def root(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j in bonded:
        parent[root(i)] = root(j)
    groups: dict[int, list[int]] = {}
    for i in range(n_atoms):
        groups.setdefault(root(i), []).append(i)
    return tuple(sorted((tuple(g) for g in groups.values()), key=lambda g: (-len(g), g)))


def fragments(symbols: Sequence[str], coords: np.ndarray) -> tuple[tuple[int, ...], ...]:
    """Connected components, largest first."""

    return _components(len(symbols), bonds(symbols, coords))


def bond_changes(
    symbols: Sequence[str], a: np.ndarray, b: np.ndarray
) -> tuple[frozenset[Bond], frozenset[Bond]]:
    """(formed, broken) going from a to b (one atom order): the pairs whose r − r_thr is
    ≥ +RESOLVED_A in one structure and ≤ −RESOLVED_A in the other."""

    ea, eb = _excess(symbols, a), _excess(symbols, b)
    return (_pairs((ea >= RESOLVED_A) & (eb <= -RESOLVED_A)),
            _pairs((ea <= -RESOLVED_A) & (eb >= RESOLVED_A)))


def same_bonding(symbols: Sequence[str], a: np.ndarray, b: np.ndarray) -> bool:
    """No resolved bond change between a and b (one atom order); not transitive: pairs only."""

    return not any(bond_changes(symbols, a, b))


def _molecule(symbols: Sequence[str], bonded: Collection[Bond]) -> Mol:
    """The connectivity molecule: one neutral atom per atom, single bonds, no implicit H, not
    sanitised."""

    from rdkit import Chem
    editable = Chem.RWMol()
    for symbol in symbols:
        atom = Chem.Atom(symbol)
        atom.SetNoImplicit(True)
        editable.AddAtom(atom)
    for i, j in sorted(bonded):
        editable.AddBond(i, j, Chem.BondType.SINGLE)
    return editable.GetMol()


def _canonical(symbols: Sequence[str], bonded: Collection[Bond]) -> tuple[str, tuple[int, ...]]:
    """(canonical SMILES, canonical ranks without tie breaking) of the connectivity molecule."""

    from rdkit import Chem
    mol = _molecule(symbols, bonded)
    return Chem.MolToSmiles(mol), tuple(Chem.CanonicalRankAtoms(mol, breakTies=False))


def _lewis(symbols: Sequence[str], bonded: Collection[Bond], group: Sequence[int]
           ) -> tuple[list[int], Mol]:
    """A fragment's atoms in canonical order (ties broken) and its molecule in that order with
    the bond orders of its Lewis structure at charge 0, single bonds where there is none. In
    canonical order the Lewis structure depends on the bond graph alone, not on the atom order
    (a resonance form, e.g. of an ylide, is picked by order). A double bond to a formally
    charged atom is one resonance form's (C=N+ of an ylide, O+=S, C-=O+): it stays single, so
    it carries no E/Z."""

    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds

    def build(atoms: Sequence[int]) -> Mol:
        index = {atom: k for k, atom in enumerate(atoms)}
        return _molecule([symbols[a] for a in atoms],
                         [(index[i], index[j]) for i, j in bonded if i in index])

    order = [group[k] for k in np.argsort(list(Chem.CanonicalRankAtoms(build(group))))]
    mol = build(order)
    perceived = Chem.Mol(mol)
    try:
        rdDetermineBonds.DetermineBondOrders(perceived, charge=0, embedChiral=False)
    except ValueError:  # no Lewis structure at charge 0 (a radical, an ion)
        return order, mol
    for bond in perceived.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE and (
                bond.GetBeginAtom().GetFormalCharge() or bond.GetEndAtom().GetFormalCharge()):
            bond.SetBondType(Chem.BondType.SINGLE)
    return order, perceived


def _placed(mol: Mol, order: Sequence[int], coords: np.ndarray) -> Mol:
    """mol with the stereo of coords (its atoms in ``order``): E/Z and tetrahedral; a
    hypervalent centre's trigonal-bipyramidal or octahedral tag is fluxional and dropped."""

    from rdkit import Chem
    from rdkit.Geometry import Point3D
    tetrahedral = (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
    placed = Chem.Mol(mol)
    conformer = Chem.Conformer(len(order))
    for k, atom in enumerate(order):
        conformer.SetAtomPosition(k, Point3D(*coords[atom]))
    placed.AddConformer(conformer, assignId=True)
    Chem.AssignStereochemistryFrom3D(placed)
    for atom in placed.GetAtoms():
        if atom.GetChiralTag() not in tetrahedral:
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    return placed


def _stereo(symbols: Sequence[str], coords: np.ndarray, bonded: Collection[Bond],
            groups: Sequence[Sequence[int]]) -> tuple[str, str]:
    """The stereo texts of the structure and of its mirror image: the sorted isomeric SMILES of
    the fragments that have stereo. '' without stereo, and with one tetrahedral centre and no
    E/Z: the mirror image is the same state, so it tells nothing (and a flattened centre that
    goes unassigned splits nothing)."""

    from rdkit import Chem
    texts: tuple[list[str], list[str]] = ([], [])
    centres = double = 0
    for group in groups:
        order, mol = _lewis(symbols, bonded, group)
        here, image = (_placed(mol, order, sign * coords) for sign in (1.0, -1.0))
        for placed, out in zip((here, image), texts, strict=True):
            isomeric = Chem.MolToSmiles(placed)
            if isomeric != Chem.MolToSmiles(placed, isomericSmiles=False):
                out.append(isomeric)
        centres += sum(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
                       for a in here.GetAtoms())
        double += sum(b.GetStereo() != Chem.BondStereo.STEREONONE for b in here.GetBonds())
    if centres <= 1 and not double:
        return "", ""
    return ".".join(sorted(texts[0])), ".".join(sorted(texts[1]))


def atom_classes(symbols: Sequence[str], bonded: Collection[Bond]) -> tuple[int, ...]:
    """Class per atom, numbered invariantly under atom permutations: an equitable partition, never
    finer than the automorphism orbits, so safe as a constraint on permutations and no more."""

    return _canonical(symbols, bonded)[1]


def state_label(symbols: Sequence[str], coords: np.ndarray) -> str:
    """Sorted fragment formulas plus the first 16 hex digits of the sha256 of the canonical
    SMILES of the bond graph, followed by the lower stereo text of the structure and its mirror
    image when there is one."""

    x = np.asarray(coords, dtype=float).reshape(-1, 3)
    bonded = bonds(symbols, x)
    groups = _components(len(symbols), bonded)
    formulas = sorted(hill_formula([symbols[i] for i in g]) for g in groups)
    stereo = min(_stereo(symbols, x, bonded, groups))
    text = _canonical(symbols, bonded)[0] + (f" {stereo}" if stereo else "")
    digest = hashlib.sha256(text.encode()).hexdigest()
    return f"{'+'.join(formulas)}_{digest[:16]}"

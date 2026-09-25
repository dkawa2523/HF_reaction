from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass
class RDKitConformerResult:
    mol: Any
    conf_id: int
    energy_kcal_mol: float
    relative_energy_kcal_mol: float
    boltzmann_weight_298K: float


def require_rdkit():
    try:
        from rdkit import Chem  # noqa: F401
        return True
    except Exception as exc:  # pragma: no cover
        raise RuntimeError("RDKit is required for this operation. Install via conda-forge.") from exc


def mol_to_xyz_block(mol, conf_id: int = 0, comment: str = "generated_by=hfauto") -> str:
    conf = mol.GetConformer(conf_id)
    lines = [str(mol.GetNumAtoms()), comment]
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol():2s} {pos.x: .8f} {pos.y: .8f} {pos.z: .8f}")
    return "\n".join(lines) + "\n"


def write_mol_conformer_xyz(mol, conf_id: int, path: str | Path, comment: str = "generated_by=hfauto") -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(mol_to_xyz_block(mol, conf_id=conf_id, comment=comment), encoding="utf-8")
    return p


def has_3d_coordinates(mol) -> bool:
    if mol.GetNumConformers() == 0:
        return False
    conf = mol.GetConformer()
    z_values = []
    max_abs = 0.0
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        z_values.append(abs(pos.z))
        max_abs = max(max_abs, abs(pos.x), abs(pos.y), abs(pos.z))
    return max_abs > 1e-6 and (max(z_values) > 1e-4 or mol.GetNumAtoms() <= 2)


def molecule_identity(mol) -> dict[str, Any]:
    from rdkit import Chem

    mol_no_h = Chem.RemoveHs(mol)
    smiles = Chem.MolToSmiles(mol_no_h, canonical=True, isomericSmiles=True)
    try:
        inchi = Chem.MolToInchi(mol_no_h)
        inchikey = Chem.InchiToInchiKey(inchi)
    except Exception:
        inchi = None
        inchikey = None
    frags = Chem.GetMolFrags(mol_no_h)
    return {
        "canonical_smiles": smiles,
        "isomeric_smiles": smiles,
        "inchi": inchi,
        "inchikey": inchikey,
        "formal_charge": int(Chem.GetFormalCharge(mol_no_h)),
        "num_atoms": int(mol.GetNumAtoms()),
        "num_heavy_atoms": int(mol_no_h.GetNumAtoms()),
        "fragment_count": len(frags),
        "has_multiple_fragments": bool(len(frags) > 1),
    }


def embed_conformers_from_smiles(
    smiles: str,
    max_conformers: int = 20,
    seed: int = 20260708,
    prune_rms_thresh: float = 0.5,
    optimize: bool = True,
    forcefield: str = "mmff",
    energy_window_kcal_mol: float | None = None,
) -> list[RDKitConformerResult]:
    from math import exp

    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Cannot parse SMILES: {smiles}")
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.pruneRmsThresh = float(prune_rms_thresh)
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=int(max_conformers), params=params))
    if not conf_ids:
        raise RuntimeError("RDKit EmbedMultipleConfs failed")

    energies: list[tuple[int, float]] = []
    ff_name = str(forcefield or "mmff").lower()
    for cid in conf_ids:
        energy = None
        if optimize:
            if ff_name in {"mmff", "mmff94", "mmff94s"}:
                try:
                    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
                    if props is not None:
                        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
                        if ff is not None:
                            ff.Minimize(maxIts=500)
                            energy = float(ff.CalcEnergy())
                except Exception:
                    energy = None
            if energy is None:
                try:
                    ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
                    ff.Minimize(maxIts=500)
                    energy = float(ff.CalcEnergy())
                except Exception:
                    energy = 0.0
        else:
            energy = 0.0
        energies.append((int(cid), float(energy)))

    min_e = min(e for _, e in energies)
    rels = [(cid, e, max(0.0, e - min_e)) for cid, e in energies]
    if energy_window_kcal_mol is not None:
        rels = [x for x in rels if x[2] <= float(energy_window_kcal_mol)] or rels[:1]
    kbt = 0.00198720425864083 * 298.15
    weights_raw = [exp(-rel / kbt) for _, _, rel in rels]
    denom = sum(weights_raw) or 1.0

    results: list[RDKitConformerResult] = []
    for (cid, e, rel), w in zip(rels, weights_raw):
        results.append(
            RDKitConformerResult(
                mol=mol,
                conf_id=cid,
                energy_kcal_mol=e,
                relative_energy_kcal_mol=rel,
                boltzmann_weight_298K=w / denom,
            )
        )
    results.sort(key=lambda r: r.relative_energy_kcal_mol)
    return results


def smiles_to_3d_mol(smiles: str, seed: int = 20260708):
    results = embed_conformers_from_smiles(smiles, max_conformers=1, seed=seed)
    return results[0].mol

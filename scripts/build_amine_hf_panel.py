"""Generate deterministic, hydrogen-complete 3D SDFs for the amine-HF panel."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


def generate(panel_path: Path, output_path: Path, subset: str, seed: int = 20260714) -> int:
    panel = json.loads(panel_path.read_text(encoding="utf-8"))
    selected = [entry for entry in panel["entries"] if subset in entry.get("subsets", [])]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(output_path))
    for index, entry in enumerate(selected):
        mol = Chem.AddHs(Chem.MolFromSmiles(entry["smiles"]))
        params = AllChem.ETKDGv3()
        params.randomSeed = int(seed + index)
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise RuntimeError(f"3D embedding failed for {entry['candidate_id']}")
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94s", maxIters=1000)
        else:
            AllChem.UFFOptimizeMolecule(mol, maxIters=1000)
        mol.SetProp("_Name", entry["name"])
        for key in [
            "candidate_id",
            "class",
            "problem_axis",
            "hypothesis",
            "reference_source",
        ]:
            mol.SetProp(key, str(entry.get(key, "")))
        mol.SetProp("canonical_smiles", Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True))
        mol.SetProp("panel_id", panel["panel_id"])
        mol.SetProp("panel_subset", subset)
        for key in ["proton_affinity_kj_mol", "gas_basicity_kj_mol"]:
            value = entry.get(key)
            if value is not None:
                mol.SetProp(key, str(value))
        writer.write(mol)
    writer.close()
    return len(selected)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel", type=Path, default=Path("examples/production/amine_hf_panel.json"))
    parser.add_argument("--subset", choices=["pilot", "full"], required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    count = generate(args.panel, args.output, args.subset)
    print(f"wrote {count} molecules to {args.output}")


if __name__ == "__main__":
    main()

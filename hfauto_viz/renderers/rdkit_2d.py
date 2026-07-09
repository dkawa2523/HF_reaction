from __future__ import annotations
from html import escape
from pathlib import Path

def fallback_svg(label: str, width: int = 300, height: int = 220) -> str:
    safe=escape(label or 'molecule')
    return f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'><rect x='1' y='1' width='{width-2}' height='{height-2}' fill='#fff' stroke='#999'/><text x='{width/2}' y='{height/2}' text-anchor='middle' font-family='Arial' font-size='14'>{safe}</text><text x='{width/2}' y='{height/2+22}' text-anchor='middle' font-family='Arial' font-size='10' fill='#777'>RDKit depiction unavailable</text></svg>"

def mol_svg_from_smiles(smiles: str | None, legend: str = '', highlight_atoms: list[int] | None = None, width: int = 300, height: int = 220) -> str:
    if not smiles: return fallback_svg(legend or 'no SMILES', width, height)
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D
        mol=Chem.MolFromSmiles(smiles)
        if mol is None: return fallback_svg(smiles, width, height)
        Chem.rdDepictor.Compute2DCoords(mol)
        drawer=rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.DrawMolecule(mol, legend=legend, highlightAtoms=highlight_atoms or [])
        drawer.FinishDrawing(); return drawer.GetDrawingText()
    except Exception:
        return fallback_svg(legend or smiles, width, height)

def write_molecule_svg(path: str | Path, smiles: str | None, legend: str = '', highlight_atoms: list[int] | None = None, width: int = 300, height: int = 220) -> Path:
    p=Path(path); p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(mol_svg_from_smiles(smiles, legend, highlight_atoms, width, height), encoding='utf-8')
    return p

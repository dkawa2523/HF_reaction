from __future__ import annotations

from pathlib import Path

from hfauto_viz.core.html import esc, write_html
from hfauto_viz.data.loaders import RunData, dataframe_to_html
from hfauto_viz.renderers.mol3d_3dmol import html_for_structure
from hfauto_viz.renderers.rdkit_2d import write_molecule_svg
from hfauto_viz.structure.xyz import read_xyz_frames


def render_molecule_dossier(run: RunData, mol_id: str, out_dir: str | Path, render_reactions: bool = True, max_reactions: int = 3, library_mode: str = 'cdn') -> dict[str, Path]:
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    mol = run.molecule_index().get(mol_id); data = mol.data if mol else {'mol_id': mol_id}
    svg = write_molecule_svg(out / 'structure_2d.svg', data.get('canonical_smiles') or data.get('isomeric_smiles'), legend=data.get('name') or mol_id)
    cand = run.read_table('candidate_summary')
    cand_sub = cand[cand['mol_id'].astype(str) == str(mol_id)] if not cand.empty and 'mol_id' in cand.columns else cand.head(0)

    frame = None
    for sp in run.species_index().values():
        if sp.data.get('mol_id') == mol_id and sp.data.get('state') == 'reactant_complex':
            frames = read_xyz_frames(run.artifact_xyz_path(sp)); frame = frames[0] if frames else None; break
    best_html = write_html(out / 'best_hf_complex.html', 'Best HF complex', html_for_structure(frame, 'mol_best', 'Best HF complex', library_mode, asset_prefix='../../assets'))

    reaction_links = ''
    if render_reactions:
        rr = run.read_table('reaction_results')
        if not rr.empty and 'mol_id' in rr.columns:
            sub = rr[rr['mol_id'].astype(str) == str(mol_id)].head(max_reactions)
            items=[]
            for _, row in sub.iterrows():
                rid=str(row.get('reaction_id'))
                items.append(f"<li><a href='../../reactions/{esc(rid)}/dossier.html'>{esc(rid)}</a> ΔG‡={esc(row.get('delta_G_act_kcal_mol',''))}</li>")
            if items:
                reaction_links = "<div class='card'><h2>Linked reactions</h2><ul>" + ''.join(items) + "</ul></div>"

    body = (
        f"<div class='grid'><div class='card'><h2>{esc(data.get('name', mol_id))}</h2><div class='svg-box'>{svg.read_text(encoding='utf-8')}</div></div>"
        f"<div class='card'><h2>Identity</h2>{dataframe_to_html(__import__('pandas').DataFrame([data]))}</div></div>"
        f"<div class='card'><h2>Candidate summary</h2>{dataframe_to_html(cand_sub)}</div>"
        f"<div class='card'><a href='best_hf_complex.html'>Best available HF complex 3D viewer</a></div>{reaction_links}"
    )
    return {'dossier': write_html(out / 'dossier.html', f'Molecule dossier: {mol_id}', body), 'structure_2d': svg, 'best_complex': best_html}

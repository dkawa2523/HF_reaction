from __future__ import annotations

from pathlib import Path

import pandas as pd

from hfauto_viz.core.html import esc, write_html
from hfauto_viz.renderers.rdkit_2d import write_molecule_svg


def write_candidate_thumbnail_table(run, candidate_df: pd.DataFrame, out_path: str | Path, top_n: int = 30) -> Path:
    out = Path(out_path); out.parent.mkdir(parents=True, exist_ok=True)
    thumb_dir = out.parent / 'thumbnails'; thumb_dir.mkdir(exist_ok=True)
    if candidate_df is None or candidate_df.empty:
        return write_html(out, 'Candidate table with thumbnails', '<p>No candidate summary available.</p>')
    mol_index = run.molecule_index()
    rows=[]
    for _, row in candidate_df.head(top_n).iterrows():
        mid = str(row.get('mol_id') or '')
        mol = mol_index.get(mid)
        data = mol.data if mol else {}
        svg = thumb_dir / f'{mid}.svg'
        write_molecule_svg(svg, data.get('canonical_smiles') or row.get('canonical_smiles'), legend=mid, width=220, height=150)
        rel = svg.relative_to(out.parent).as_posix()
        dossier = f"../molecules/{mid}/dossier.html"
        rows.append({
            'structure': f"<a href='{esc(dossier)}'><img class='thumbnail' src='{esc(rel)}' alt='{esc(mid)}'></a>",
            'mol_id': f"<a href='{esc(dossier)}'>{esc(mid)}</a>",
            'site_id': esc(row.get('site_id','')),
            'hf_n': esc(row.get('hf_n','')),
            'quality_tier': esc(row.get('quality_tier','')),
            'confidence_score': esc(row.get('confidence_score','')),
            'scavenger_score': esc(row.get('scavenger_score','')),
            'activation_score': esc(row.get('activation_score','')),
            'recommended_next_action': esc(row.get('recommended_next_action','')),
        })
    header=''.join(f'<th>{esc(k)}</th>' for k in rows[0]) if rows else ''
    body=''.join('<tr>'+''.join(f'<td>{v}</td>' for v in r.values())+'</tr>' for r in rows)
    html=f"<div class='card'><h2>Candidate ranking with thumbnails</h2><table class='data-table'><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table></div>"
    return write_html(out,'Candidate ranking with thumbnails',html)

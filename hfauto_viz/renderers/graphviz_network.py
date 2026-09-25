from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Any

import pandas as pd

from hfauto_viz.core.html import esc, write_html
from hfauto_viz.renderers.rdkit_2d import write_molecule_svg


def _safe_id(value: str) -> str:
    return 'n_' + ''.join(ch if ch.isalnum() else '_' for ch in str(value))[:120]


def _dot_quote(value: Any) -> str:
    return '"' + str(value).replace('\\','/').replace('"','\\"') + '"'


def _fallback_svg(dot: str, out: Path) -> Path:
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        f"<svg xmlns='http://www.w3.org/2000/svg' width='1200' height='720'><rect width='100%' height='100%' fill='white'/><foreignObject x='20' y='20' width='1160' height='680'><pre xmlns='http://www.w3.org/1999/xhtml' style='font-size:12px;white-space:pre-wrap'>{esc(dot)}</pre></foreignObject></svg>",
        encoding='utf-8',
    )
    return out


def build_reaction_network(reaction_df: pd.DataFrame, molecules: dict[str, dict[str, Any]], out_dir: str | Path, top_n: int = 25) -> dict[str, Path]:
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    node_dir = out / 'node_images'; node_dir.mkdir(exist_ok=True)
    rows = reaction_df.head(top_n).to_dict('records') if reaction_df is not None and not reaction_df.empty else []
    lines = [
        'digraph hfauto_reaction_network {',
        'rankdir=LR;',
        'graph [bgcolor=white, pad=0.2, nodesep=0.55, ranksep=0.75];',
        'node [shape=box, style="rounded,filled", fillcolor="#f7fbff", fontsize=10, fontname="Helvetica", labelloc=b];',
        'edge [fontsize=9, fontname="Helvetica", color="#777777", arrowsize=0.7];',
    ]
    seen_molecules: set[str] = set()
    for row in rows:
        rid = str(row.get('reaction_id') or '')
        mid = str(row.get('mol_id') or 'candidate')
        hf = str(row.get('hf_n',''))
        if mid not in seen_molecules:
            img = node_dir / f'{_safe_id(mid)}.svg'
            write_molecule_svg(img, molecules.get(mid, {}).get('canonical_smiles') or row.get('canonical_smiles'), legend=mid, width=220, height=150)
            rel_img = img.resolve().as_posix()
            lines.append(f'{_safe_id(mid)} [label={_dot_quote(mid)}, image={_dot_quote(rel_img)}, URL={_dot_quote("../molecules/"+mid+"/dossier.html")}, tooltip={_dot_quote(mid)}];')
            seen_molecules.add(mid)
        rc = _safe_id(rid + '_rc'); ts = _safe_id(rid + '_ts'); ip = _safe_id(rid + '_ip')
        dossier = f'../reactions/{rid}/dossier.html'
        assoc = row.get('delta_G_assoc_pressure_corrected_kcal_mol', row.get('delta_G_assoc_kcal_mol',''))
        act = row.get('delta_G_act_kcal_mol','')
        qtier = row.get('quality_tier','')
        rc_label = _dot_quote(f"B···(HF){hf}\n{qtier}")
        ts_label = _dot_quote(f"TS\n{qtier}")
        ip_label = _dot_quote(f"ion pair\n{qtier}")
        lines += [
            f'{rc} [label={rc_label}, fillcolor="#e7f1ff", URL={_dot_quote(dossier)}];',
            f'{ts} [label={ts_label}, shape=octagon, fillcolor="#fff3cd", URL={_dot_quote(dossier)}];',
            f'{ip} [label={ip_label}, fillcolor="#e9f7ef", URL={_dot_quote(dossier)}];',
            f'{_safe_id(mid)} -> {rc} [label={_dot_quote("assoc "+str(assoc))}, URL={_dot_quote(dossier)}];',
            f'{rc} -> {ts} [label={_dot_quote("ΔG‡="+str(act))}, color="#d08b00", URL={_dot_quote(dossier)}];',
            f'{ts} -> {ip} [label="PT", URL={_dot_quote(dossier)}];',
        ]
    lines.append('}')
    dot = '\n'.join(lines)
    dot_path = out / 'reaction_network.dot'; dot_path.write_text(dot, encoding='utf-8')
    svg_path = out / 'reaction_network.svg'
    if shutil.which('dot'):
        try:
            subprocess.run(['dot', '-Tsvg', str(dot_path), '-o', str(svg_path)], check=True, timeout=30)
        except Exception:
            _fallback_svg(dot, svg_path)
    else:
        _fallback_svg(dot, svg_path)
    html_path = write_html(out / 'reaction_network.html', 'Reaction network', f"<div class='card'><p><a href='reaction_network.dot'>DOT</a> | <a href='reaction_network.svg'>SVG</a> | <a href='reaction_network_interactive.html'>Interactive Cytoscape view</a></p><div class='svg-box'>{svg_path.read_text(encoding='utf-8')}</div></div>")
    return {'dot': dot_path, 'svg': svg_path, 'html': html_path}


def render_reaction_network(candidate_df, reaction_df, out_dir, max_reactions: int = 25):
    molecules = {}
    if candidate_df is not None and not candidate_df.empty and 'mol_id' in candidate_df.columns:
        for _, row in candidate_df.drop_duplicates('mol_id').iterrows():
            molecules[str(row.get('mol_id'))] = row.to_dict()
    return build_reaction_network(reaction_df, molecules, out_dir, top_n=max_reactions)

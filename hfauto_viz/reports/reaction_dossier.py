from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from hfauto_viz.core.html import dataframe_to_html, esc, write_html
from hfauto_viz.data.loaders import RunData
from hfauto_viz.renderers.coordinate_plot import write_reaction_coordinate_plot
from hfauto_viz.renderers.mode_view import write_imaginary_mode_view
from hfauto_viz.renderers.mol3d_3dmol import html_for_structure, write_animation_html
from hfauto_viz.renderers.plotly_energy import (
    plotly_energy_div,
    write_energy_profile,
    write_energy_profile_static_svg,
)
from hfauto_viz.renderers.publication_bundle import write_publication_bundle
from hfauto_viz.structure.bond_changes import (
    annotation_shapes_for_frame,
    bond_change_table,
    write_bond_change_outputs,
)
from hfauto_viz.structure.path_sources import reaction_path_frames_from_run
from hfauto_viz.structure.xyz import Frame, frames_to_xyz, read_xyz_frames


def _first_frame(path):
    frames = read_xyz_frames(path)
    return frames[0] if frames else None


def _frame_annotations(frame: Frame | None, reaction_data: dict[str, Any]) -> list[dict[str, Any]]:
    return annotation_shapes_for_frame(frame, reaction_data, include_labels=True, include_bond_lines=True)


def _path_frame_annotations(frames: list[tuple[str, Frame]], reaction_data: dict[str, Any]) -> list[dict[str, Any]]:
    annotations: list[dict[str, Any]] = []
    # Annotate first, TS-like middle, and last frames so the viewer stays readable.
    idxs = sorted({0, max(0, len(frames)//2), max(0, len(frames)-1)})
    for idx in idxs:
        if idx < len(frames):
            annotations.extend({**ann, "frame": idx} for ann in _frame_annotations(frames[idx][1], reaction_data))
    return annotations


def render_reaction_dossier(run: RunData, reaction_id: str, out_dir: str | Path, library_mode: str = 'cdn') -> dict[str, Path]:
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    df = run.read_table('reaction_results')
    if df.empty:
        df = run.read_table('reaction_thermo')
    row: dict[str, Any] = {}
    if not df.empty and 'reaction_id' in df.columns:
        sub = df[df['reaction_id'].astype(str) == str(reaction_id)]
        if not sub.empty:
            row = sub.iloc[0].to_dict()
    rxn = run.reaction_index().get(reaction_id)
    data = rxn.data if rxn else {}
    merged_reaction_data = {**data, **row}
    reactant_id = data.get('reactant_species_id') or row.get('reactant_species_id')
    product_id = data.get('product_species_id') or row.get('product_species_id')
    ts_id = data.get('ts_species_id') or row.get('ts_species_id')
    if not ts_id:
        for a in run.latest_artifacts('reaction_validated'):
            if a.data.get('reaction_id') == reaction_id:
                ts_id = (a.data.get('ts') or {}).get('species_id') or a.data.get('ts_species_id')
                break
    sp_idx = run.species_index()
    reactant = _first_frame(run.artifact_xyz_path(sp_idx.get(str(reactant_id))))
    product = _first_frame(run.artifact_xyz_path(sp_idx.get(str(product_id))))
    ts = _first_frame(run.artifact_xyz_path(sp_idx.get(str(ts_id))))

    asset_prefix = '../../assets'
    paths: dict[str, Path] = {}
    paths['energy_profile'] = write_energy_profile(row, out / 'energy_profile.html') if row else write_html(out / 'energy_profile.html', 'Energy profile', '<p>No thermochemistry row.</p>')
    if row:
        paths.update(write_energy_profile_static_svg(row, out / 'energy_profile.svg'))
    paths['energy'] = paths['energy_profile']

    # 3D structure viewers with chemically meaningful annotations.
    paths['reactant_view'] = write_html(out / 'reactant_view.html', 'Reactant 3D', html_for_structure(reactant, 'reactant_view', 'Reactant', library_mode, asset_prefix=asset_prefix, annotations=_frame_annotations(reactant, merged_reaction_data)))
    paths['ts_view'] = write_html(out / 'ts_view.html', 'TS 3D', html_for_structure(ts, 'ts_view', 'TS', library_mode, asset_prefix=asset_prefix, annotations=_frame_annotations(ts, merged_reaction_data)))
    paths['product_view'] = write_html(out / 'product_view.html', 'Product 3D', html_for_structure(product, 'product_view', 'Product', library_mode, asset_prefix=asset_prefix, annotations=_frame_annotations(product, merged_reaction_data)))

    frames, path_source = reaction_path_frames_from_run(run, reaction_id, reactant, ts, product, interpolation_points=12)
    path_xyz = out / 'path.xyz'
    path_xyz.write_text(frames_to_xyz(frames), encoding='utf-8')
    paths['path_xyz'] = path_xyz
    paths['path_source_json'] = out / 'path_source.json'
    paths['path_source_json'].write_text(__import__('json').dumps(path_source, ensure_ascii=False, indent=2), encoding='utf-8')
    paths['path_animation'] = write_animation_html(frames, out / 'path_animation.html', title=f'Path animation: {reaction_id}', library_mode=library_mode, asset_prefix=asset_prefix, annotations=_path_frame_annotations(frames, merged_reaction_data))

    coord_paths = write_reaction_coordinate_plot(frames, merged_reaction_data, out / 'reaction_coordinate.html', library_mode=library_mode, asset_prefix=asset_prefix)
    paths['reaction_coordinate'] = coord_paths['html']
    paths['reaction_coordinate_csv'] = coord_paths['csv']

    bc_rows = bond_change_table(reactant, ts, product, merged_reaction_data)
    paths.update(write_bond_change_outputs(bc_rows, out / 'bond_change_annotations'))

    paths.update(write_imaginary_mode_view(ts, merged_reaction_data, out / 'imaginary_mode.html', reactant=reactant, product=product, library_mode=library_mode, asset_prefix=asset_prefix, title=f'Imaginary mode: {reaction_id}'))
    paths.update(write_publication_bundle(out))

    quality_bits = {
        'quality_tier': row.get('quality_tier', ''),
        'confidence_score': row.get('confidence_score', ''),
        'delta_G_assoc_pressure_corrected_kcal_mol': row.get('delta_G_assoc_pressure_corrected_kcal_mol', ''),
        'delta_G_act_kcal_mol': row.get('delta_G_act_kcal_mol', ''),
        'k_corrected_s-1': row.get('k_corrected_s-1', ''),
        'path_source': path_source.get('source'),
        'real_path_frames': path_source.get('real_path_frames'),
        'n_path_frames': path_source.get('n_frames'),
    }
    links = [
        ('Energy profile', 'energy_profile.html'),
        ('Static energy SVG', 'energy_profile.svg'),
        ('3D path animation', 'path_animation.html'),
        ('Reaction coordinate plot', 'reaction_coordinate.html'),
        ('Imaginary mode / TS vector', 'imaginary_mode.html'),
        ('Bond-change annotations', 'bond_change_annotations.csv'),
        ('Reactant 3D', 'reactant_view.html'),
        ('TS 3D', 'ts_view.html'),
        ('Product 3D', 'product_view.html'),
        ('Publication figure bundle', 'publication_figures.html'),
    ]
    link_html = ''.join(f"<li><a href='{esc(href)}'>{esc(label)}</a></li>" for label, href in links)
    warn = "" if path_source.get('real_path_frames') else "<div class='warn'><b>Path-frame note:</b> This dossier uses an interpolated reactant→TS→product path unless a real IRC/NEB trajectory artifact is available.</div>"
    body = (
        f"{warn}<div class='card'><h2>Reaction {esc(reaction_id)}</h2>"
        "<p class='muted'>The dossier links energy, annotated 3D structures, path animation, reaction-coordinate diagnostics, imaginary-mode arrows and exportable figure files.</p>"
        f"<ul>{link_html}</ul></div>"
        f"<div class='grid'><div class='card'><h2>Review summary</h2>{dataframe_to_html(pd.DataFrame([quality_bits]))}</div><div class='card'><h2>Bond-change distances</h2>{dataframe_to_html(pd.DataFrame(bc_rows))}</div></div>"
        f"<div class='card'>{plotly_energy_div(row) if row else ''}</div>"
        f"<div class='card'><h2>Reaction row</h2>{dataframe_to_html(pd.DataFrame([row]))}</div>"
    )
    paths['dossier'] = write_html(out / 'dossier.html', f'Reaction dossier: {reaction_id}', body)
    return paths

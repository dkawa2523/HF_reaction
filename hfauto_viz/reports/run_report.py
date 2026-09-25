from __future__ import annotations

from pathlib import Path

from hfauto_viz.core.assets import write_placeholder_assets
from hfauto_viz.core.html import dataframe_to_html, esc, write_html
from hfauto_viz.core.viz_manifest import VizManifest, VizRecord
from hfauto_viz.data.loaders import RunData
from hfauto_viz.renderers.chemiscope_export import write_chemiscope_bundle
from hfauto_viz.renderers.cytoscape_network import (
    write_cytoscape_network,
    write_cytoscape_network_linked,
)
from hfauto_viz.renderers.graphviz_network import build_reaction_network
from hfauto_viz.renderers.plotly_energy import write_energy_comparison
from hfauto_viz.renderers.plotly_screening import (
    write_failure_heatmap,
    write_quality_funnel,
    write_scatter,
)
from hfauto_viz.renderers.ranking_table import write_candidate_thumbnail_table
from hfauto_viz.reports.molecule_dossier import render_molecule_dossier
from hfauto_viz.reports.reaction_dossier import render_reaction_dossier
from hfauto_viz.structure.xyz import read_xyz_frames


def _failure_dataframe(run):
    import pandas as pd
    rows = []
    for a in run.manifest.artifacts:
        if a.status.status == 'failed':
            rows.append({
                'artifact_id': a.artifact_id,
                'artifact_type': a.artifact_type,
                'stage': (a.provenance or {}).get('stage') or a.artifact_id.split('_')[0],
                'category': a.status.category or 'unknown',
                'reason': a.status.reason or '',
            })
    return pd.DataFrame(rows)


def _quality_funnel_counts(run):
    return {
        'molecules': len(run.latest_artifacts('molecule')),
        'sites': len(run.latest_artifacts('site')),
        'conformers': len(run.latest_artifacts('conformer')),
        'HF species': len(run.latest_artifacts('species')),
        'preopt geometries': len(run.latest_artifacts('species_preopt')),
        'DFT optimized': len(run.latest_artifacts('species_optimized')),
        'TS paths': len(run.latest_artifacts('ts_path')),
        'IRC records': len(run.latest_artifacts('irc')),
        'thermo records': len(run.latest_artifacts('thermo')),
        'kinetics records': len(run.latest_artifacts('kinetics')),
    }


def _structures_for_chemiscope(run, df):
    structures = {}; sp_idx = run.species_index()
    for _, row in df.iterrows():
        sid = str(row.get('reactant_species_id') or row.get('species_id') or '')
        sp = sp_idx.get(sid)
        if sp:
            frames = read_xyz_frames(run.artifact_xyz_path(sp))
            if frames:
                structures[sid] = frames[0]
    return structures


def _rel(path: str | Path, out: Path) -> str:
    p = Path(path)
    try:
        if p.is_absolute():
            return p.relative_to(out).as_posix()
    except Exception:
        pass
    return str(path)


def generate_visualizations(run_dir: str | Path, out_dir: str | Path | None = None, config: dict | None = None):
    import pandas as pd

    config = config or {}
    run = RunData(run_dir)
    out = Path(out_dir) if out_dir else Path(run_dir) / '14_viz'
    out.mkdir(parents=True, exist_ok=True)
    viz = VizManifest(run.run_id, out)

    library_mode = str(config.get('library_mode') or config.get('asset_mode') or 'cdn')
    asset_mode = str(config.get('asset_mode') or 'local')
    assets_dir = write_placeholder_assets(out)
    viz.add(VizRecord('viz_assets_bundle', 'asset_bundle', str(assets_dir), renderer='hfauto_viz_assets', data={'asset_mode': asset_mode, 'library_mode': library_mode}))

    if run.has_dummy_or_fallback():
        viz.warn('Fallback/dummy calculation data detected. Use plots for workflow review unless real calculation flags are present.')

    reaction_df = run.read_table('reaction_results')
    candidate_df = run.read_table('candidate_summary')
    act_df = run.read_table('rank_activation')
    cluster_df = run.read_table('cluster_risk')

    screening = out / 'screening'; screening.mkdir(exist_ok=True)
    p = write_scatter(candidate_df, screening / 'scavenger_vs_activation.html', 'Scavenger vs activation', 'scavenger_score', 'activation_score', color='quality_tier', size='confidence_score')
    viz.add(VizRecord('viz_scavenger_vs_activation', 'scatter', str(p), renderer='plotly'))
    if not reaction_df.empty:
        p = write_energy_comparison(reaction_df, screening / 'energy_profile_top_candidates.html', top_n=int(config.get('top_energy_profiles', 12)))
        viz.add(VizRecord('viz_energy_profile_top_candidates', 'energy_profile_comparison', str(p), renderer='plotly_energy'))
    for name, df, x, y, color, size in [
        ('assoc_vs_barrier', reaction_df, 'delta_G_assoc_pressure_corrected_kcal_mol', 'delta_G_act_kcal_mol', 'quality_tier', 'confidence_score'),
        ('activation_map', act_df if not act_df.empty else reaction_df, 'delta_r_HF_A', 'delta_nu_HF_cm1', 'site_type', 'confidence_score'),
        ('cluster_risk_map', cluster_df if not cluster_df.empty else reaction_df, 'delta_G_assoc_pressure_corrected_kcal_mol', 'cluster_growth_score', 'ehs_review_flag', 'confidence_score'),
    ]:
        p = write_scatter(df, screening / f'{name}.html', name, x, y, color=color, size=size)
        viz.add(VizRecord(f'viz_{name}', 'scatter', str(p), renderer='plotly'))
    p = write_quality_funnel(_quality_funnel_counts(run), screening / 'quality_funnel.html')
    viz.add(VizRecord('viz_quality_funnel', 'quality_funnel', str(p), renderer='plotly'))
    p = write_failure_heatmap(_failure_dataframe(run), screening / 'failure_heatmap.html')
    viz.add(VizRecord('viz_failure_heatmap', 'failure_heatmap', str(p), renderer='plotly'))

    mol_rows = candidate_df if not candidate_df.empty else pd.DataFrame([a.data for a in run.latest_artifacts('molecule')])
    max_mols = int(config.get('max_molecule_dossiers', 8))
    for mid in mol_rows.get('mol_id', pd.Series(dtype=str)).dropna().astype(str).drop_duplicates().head(max_mols):
        paths = render_molecule_dossier(run, mid, out / 'molecules' / mid, max_reactions=int(config.get('max_reactions_per_molecule', 2)), library_mode=library_mode)
        viz.add(VizRecord(f'viz_molecule_{mid}', 'molecule_dossier', str(paths['dossier']), parents=[mid], renderer='3dmol+rdkit'))

    p = write_candidate_thumbnail_table(run, candidate_df, screening / 'candidate_ranking_thumbnails.html', top_n=int(config.get('ranking_thumbnail_rows', 30)))
    viz.add(VizRecord('viz_candidate_ranking_thumbnails', 'candidate_ranking_with_thumbnails', str(p), renderer='rdkit+html'))

    max_rxns = int(config.get('max_reaction_dossiers', 8))
    if not reaction_df.empty and 'reaction_id' in reaction_df.columns:
        for rid in reaction_df['reaction_id'].dropna().astype(str).drop_duplicates().head(max_rxns):
            paths = render_reaction_dossier(run, rid, out / 'reactions' / rid, library_mode=library_mode)
            viz.add(VizRecord(f'viz_reaction_{rid}', 'reaction_dossier', str(paths['dossier']), parents=[rid], renderer='3dmol+plotly'))
            if 'reaction_coordinate' in paths:
                viz.add(VizRecord(f'viz_reaction_coordinate_{rid}', 'reaction_coordinate', str(paths['reaction_coordinate']), parents=[rid], renderer='plotly_coordinate'))
            for key, viz_type, renderer in [
                ('imaginary_mode_html', 'imaginary_mode_viewer', '3dmol_mode'),
                ('bond_change_json', 'bond_change_annotations', 'hfauto_bond_change'),
                ('publication_figures_html', 'publication_figure_bundle', 'html_static_export'),
                ('path_source_json', 'reaction_path_source', 'hfauto_path_source'),
            ]:
                if key in paths:
                    viz.add(VizRecord(f'viz_{key}_{rid}', viz_type, str(paths[key]), parents=[rid], renderer=renderer))

    network_paths = build_reaction_network(reaction_df, {mid: a.data for mid, a in run.molecule_index().items()}, out / 'networks', top_n=int(config.get('max_network_reactions', 20)))
    for kind, pth in network_paths.items():
        viz.add(VizRecord(f'viz_reaction_network_{kind}', f'reaction_network_{kind}', str(pth), renderer='graphviz'))
    cyt_paths = write_cytoscape_network(reaction_df, out / 'networks', top_n=int(config.get('max_network_reactions', 50)), library_mode=library_mode)
    for kind, pth in cyt_paths.items():
        viz.add(VizRecord(f'viz_reaction_network_{kind}', f'reaction_network_{kind}', str(pth), renderer='cytoscape'))
    linked_paths = write_cytoscape_network_linked(run, reaction_df, out / 'networks', top_n=int(config.get('max_network_reactions', 50)), library_mode=library_mode)
    for kind, pth in linked_paths.items():
        viz.add(VizRecord(f'viz_reaction_network_{kind}', f'reaction_network_{kind}', str(pth), renderer='cytoscape+3dmol'))

    chem_dir = out / 'chemiscope'
    chem_rows = reaction_df.head(int(config.get('max_chemiscope_structures', 200))).copy()
    if not chem_rows.empty:
        chem_paths = write_chemiscope_bundle(chem_rows, _structures_for_chemiscope(run, chem_rows), chem_dir / 'hf_screening.json.gz')
        for kind, pth in chem_paths.items():
            viz.add(VizRecord(f'viz_{kind}', kind, str(pth), renderer='chemiscope_fallback'))

    warning = '<div class="warn"><b>Fallback/dummy warning:</b> ' + '<br>'.join(esc(w) for w in viz.warnings) + '</div>' if viz.warnings else ''
    links = ''.join(
        f"<li><a href='{esc(_rel(r.path, out))}'>{esc(r.viz_type)}</a> <span class='muted'>({esc(r.renderer or '')})</span></li>"
        for r in viz.records[:80]
    )
    rows = candidate_df.head(20) if not candidate_df.empty else reaction_df.head(20)
    body = (
        f"{warning}"
        f"<div class='card'><h2>Production visualization report</h2><p><b>run_id:</b> {esc(run.run_id)}</p><p><b>latest stage:</b> {esc(run.manifest.stage)}</p><p><b>library_mode:</b> {esc(library_mode)}</p><p class='muted'>Includes linked reaction-network/3D viewer, annotated bond-change structure views, imaginary-mode arrows, path-source tracking, and publication figure bundles.</p></div>"
        f"<div class='grid'><div class='card'><h2>Visualizations</h2><ul>{links}</ul></div><div class='card'><h2>Key review pages</h2><ul><li><a href='screening/candidate_ranking_thumbnails.html'>Ranking with thumbnails</a></li><li><a href='networks/reaction_network_interactive.html'>Interactive network</a></li><li><a href='chemiscope/hf_screening.html'>Chemiscope export</a></li><li><a href='assets/asset_status.html'>Asset status</a></li></ul></div></div>"
        f"<div class='card'><h2>Preview</h2>{dataframe_to_html(rows)}</div>"
        f"<div class='card'><iframe src='networks/reaction_network_linked.html' style='width:100%;height:760px;border:1px solid #ddd'></iframe></div>"
    )
    report_path = write_html(out / 'report.html', f'Production visualization report: {run.run_id}', body)
    viz.add(VizRecord('viz_run_report', 'run_report', str(report_path), renderer='html'))
    viz.write()
    return report_path, viz


def write_run_report(run_dir: str | Path, out_dir: str | Path | None = None, top_n: int = 8, asset_mode: str = 'cdn'):
    return generate_visualizations(run_dir, out_dir=out_dir, config={'max_molecule_dossiers': top_n, 'max_reaction_dossiers': top_n, 'library_mode': asset_mode, 'asset_mode': asset_mode})[0]


def render_run_visual_report(run_data, out_dir, top_n: int = 8, **kwargs):
    run_dir = getattr(run_data, 'run_dir', run_data)
    return write_run_report(run_dir, out_dir=out_dir, top_n=top_n, asset_mode=kwargs.get('asset_mode', kwargs.get('library_mode', 'cdn')))

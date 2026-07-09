from __future__ import annotations

import json
from pathlib import Path

from hfauto.core.config import load_yaml
from hfauto.workflow.runner import run_pipeline
from hfauto_viz.data.loaders import RunData
from hfauto_viz.reports.run_report import generate_visualizations
from hfauto_viz.renderers.coordinate_plot import reaction_coordinate_dataframe
from hfauto_viz.renderers.cytoscape_network import cytoscape_elements
from hfauto_viz.structure.xyz import Atom


def test_phase85_coordinate_dataframe_auto_detects_bhf():
    frame = [Atom('N',0,0,0), Atom('H',1.0,0,0), Atom('F',2.0,0,0)]
    df = reaction_coordinate_dataframe([('rc', frame)], {})
    assert {'r_BH_A','r_HF_A','q_BH_minus_HF_A'} <= set(df.columns)
    assert abs(float(df.iloc[0]['r_BH_A']) - 1.0) < 1e-6
    assert abs(float(df.iloc[0]['r_HF_A']) - 1.0) < 1e-6


def test_phase85_cytoscape_elements_contract():
    import pandas as pd
    df = pd.DataFrame([{'reaction_id':'rxn1','mol_id':'mol1','hf_n':2,'delta_G_act_kcal_mol':5.5,'quality_tier':'Q4'}])
    elements = cytoscape_elements(df)
    assert any(e['data'].get('id') == 'mol1' for e in elements)
    assert any(e['data'].get('source') for e in elements)


def test_phase85_hardened_viz_outputs(tmp_path: Path):
    cfg = load_yaml('configs/pipelines/phase7_full_publicdb_rank_offline.yaml')
    cfg['run_root'] = str(tmp_path / 'runs')
    run_dir = Path(run_pipeline(cfg, 'phase85_base'))
    report, viz = generate_visualizations(run_dir, out_dir=run_dir / '15_viz', config={
        'library_mode': 'local',
        'asset_mode': 'local',
        'max_molecule_dossiers': 2,
        'max_reaction_dossiers': 2,
        'top_energy_profiles': 2,
        'max_network_reactions': 4,
    })
    assert report.exists()
    out = run_dir / '15_viz'
    assert (out / 'assets' / 'asset_manifest.json').exists()
    assert (out / 'assets' / '3Dmol-min.js').exists()
    assert (out / 'screening' / 'candidate_ranking_thumbnails.html').exists()
    assert (out / 'networks' / 'reaction_network_interactive.html').exists()
    assert (out / 'networks' / 'reaction_network_interactive.json').exists()
    assert list((out / 'reactions').glob('*/reaction_coordinate.html'))
    assert list((out / 'reactions').glob('*/energy_profile.svg'))
    assert (out / 'chemiscope' / 'hf_screening.html').exists()
    vm = json.loads((out / 'viz_manifest.json').read_text(encoding='utf-8'))
    types = {v['viz_type'] for v in vm['visualizations']}
    assert 'candidate_ranking_with_thumbnails' in types
    assert 'reaction_network_cytoscape_html' in types or 'reaction_network_cytoscape_json' in types

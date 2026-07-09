from pathlib import Path

from hfauto.workflow.runner import run_pipeline
from hfauto.core.config import load_yaml
from hfauto_viz.data.loaders import RunData
from hfauto_viz.reports.reaction_dossier import render_reaction_dossier
from hfauto_viz.structure.bond_changes import bond_change_table, annotation_shapes_for_frame
from hfauto_viz.structure.modes import approximate_proton_transfer_mode
from hfauto_viz.structure.xyz import Atom


def test_phase11_bond_and_mode_helpers():
    reactant = [Atom('N',0,0,0), Atom('H',1.6,0,0), Atom('F',2.5,0,0)]
    ts = [Atom('N',0,0,0), Atom('H',1.1,0,0), Atom('F',2.0,0,0)]
    product = [Atom('N',0,0,0), Atom('H',1.0,0,0), Atom('F',2.6,0,0)]
    rxn = {'reaction_coordinate': {'atoms': {'base_atom':0,'transfer_h':1,'leaving_f':2}}}
    rows = bond_change_table(reactant, ts, product, rxn)
    assert {r['label'] for r in rows} == {'forming_B_H', 'breaking_H_F'}
    ann = annotation_shapes_for_frame(ts, rxn)
    assert any(a['type'] == 'line' for a in ann)
    vec = approximate_proton_transfer_mode(ts, rxn, reactant=reactant, product=product)
    assert any(v.role == 'transfer_h' and v.norm > 0 for v in vec)


def test_phase11_pipeline_visual_outputs(tmp_path):
    cfg = load_yaml('configs/pipelines/phase11_visualization_refinement_offline.yaml')
    cfg['run_root'] = str(tmp_path / 'runs')
    run_dir = run_pipeline(cfg, 'phase11_test')
    viz = Path(run_dir) / '15_viz'
    assert (viz / 'report.html').exists()
    assert (viz / 'networks' / 'reaction_network_linked.html').exists()
    assert (viz / 'viz_manifest.json').exists()
    reaction_dirs = list((viz / 'reactions').glob('*'))
    assert reaction_dirs, 'expected reaction dossiers'
    first = reaction_dirs[0]
    assert (first / 'imaginary_mode.html').exists()
    assert (first / 'bond_change_annotations.csv').exists()
    assert (first / 'publication_figures.html').exists()

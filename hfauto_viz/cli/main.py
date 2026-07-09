from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer
from rich.console import Console

from hfauto_viz.core.assets import write_placeholder_assets
from hfauto_viz.data.loaders import RunData
from hfauto_viz.renderers.chemiscope_export import write_chemiscope_bundle, write_chemiscope_like
from hfauto_viz.renderers.cytoscape_network import write_cytoscape_network
from hfauto_viz.renderers.graphviz_network import build_reaction_network
from hfauto_viz.reports.molecule_dossier import render_molecule_dossier
from hfauto_viz.reports.reaction_dossier import render_reaction_dossier
from hfauto_viz.reports.run_report import generate_visualizations

app = typer.Typer(help='Visualization companion CLI for completed hfauto runs')
console = Console()


@app.command()
def report(
    run_dir: Path,
    out_dir: Optional[Path] = typer.Option(None, '--out'),
    top_n: int = 8,
    library_mode: str = typer.Option('cdn', help='cdn, local, auto, or none'),
):
    path, viz = generate_visualizations(run_dir, out_dir=out_dir, config={'max_molecule_dossiers': top_n, 'max_reaction_dossiers': top_n, 'library_mode': library_mode, 'asset_mode': library_mode})
    console.print(f'[green]wrote[/green] {path}')
    console.print(f'visualizations: {len(viz.records)}')


@app.command('bundle-assets')
def bundle_assets(out_dir: Path):
    assets = write_placeholder_assets(out_dir)
    console.print(f'[green]wrote[/green] {assets}')


@app.command()
def molecule(run_dir: Path, mol_id: str, out_dir: Optional[Path] = typer.Option(None, '--out'), library_mode: str = 'cdn'):
    run = RunData(run_dir)
    paths = render_molecule_dossier(run, mol_id, out_dir or Path(run_dir) / '14_viz' / 'molecules' / mol_id, library_mode=library_mode)
    console.print_json(data={k: str(v) for k, v in paths.items()})


@app.command()
def reaction(run_dir: Path, reaction_id: str, out_dir: Optional[Path] = typer.Option(None, '--out'), library_mode: str = 'cdn'):
    run = RunData(run_dir)
    paths = render_reaction_dossier(run, reaction_id, out_dir or Path(run_dir) / '14_viz' / 'reactions' / reaction_id, library_mode=library_mode)
    console.print_json(data={k: str(v) for k, v in paths.items()})


@app.command()
def network(run_dir: Path, out_dir: Optional[Path] = typer.Option(None, '--out'), library_mode: str = 'cdn'):
    run = RunData(run_dir)
    out = out_dir or Path(run_dir) / '14_viz' / 'networks'
    paths = build_reaction_network(run.read_table('reaction_results'), {mid: a.data for mid, a in run.molecule_index().items()}, out)
    paths.update(write_cytoscape_network(run.read_table('reaction_results'), out, library_mode=library_mode))
    console.print_json(data={k: str(v) for k, v in paths.items()})


@app.command()
def chemiscope(run_dir: Path, out: Path = typer.Option(...)):
    run = RunData(run_dir)
    rr = run.read_table('reaction_results')
    paths = write_chemiscope_bundle(rr, {}, out)
    console.print_json(data={k: str(v) for k, v in paths.items()})


@app.command()
def status(run_dir: Path):
    run = RunData(run_dir)
    console.print(f'run_id: {run.run_id}')
    console.print(f'latest_stage: {run.manifest.stage}')
    console.print(f'latest_manifest: {run.manifest_path}')
    console.print(f'fallback_or_dummy: {run.has_dummy_or_fallback()}')


if __name__ == '__main__':
    app()

# Phase 8: Visualization package

Phase 8 adds `hfauto_viz`, a read-only visualization companion for completed `hfauto` runs.

## Key design points

- `hfauto_viz` is independent from simulation stages.
- It reads `manifest.json`, CSV/JSONL tables, and XYZ/SDF files.
- It writes HTML/SVG/DOT/JSON visualization artifacts under `<run>/14_viz` or, when used as a pipeline stage, `<run>/15_viz`.
- It does not change calculation results or promote quality tiers.
- It shows fallback/dummy warnings when a run contains non-production calculation artifacts.

## Main commands

```bash
hfauto-viz report runs/<run_id> --out runs/<run_id>/14_viz --top-n 8
hfauto-viz molecule runs/<run_id> mol00012
hfauto-viz reaction runs/<run_id> rxn_...
hfauto-viz network runs/<run_id>
hfauto-viz chemiscope runs/<run_id> --out runs/<run_id>/14_viz/chemiscope/hf_screening.json.gz
```

## Pipeline stage

`configs/pipelines/phase8_full_viz_offline.yaml` appends a final `viz` stage to the Phase 7 offline pipeline. The stage writes a `visualization_bundle` artifact so the latest manifest records visualization outputs.

## Outputs

```text
14_viz or 15_viz/
  report.html
  viz_manifest.json
  screening/
    energy_profile_top_candidates.html
    assoc_vs_barrier.html
    activation_map.html
    cluster_risk_map.html
    quality_funnel.html
    failure_heatmap.html
    scavenger_vs_activation.html
  molecules/<mol_id>/
    dossier.html
    structure_2d.svg
    best_hf_complex.html
  reactions/<reaction_id>/
    dossier.html
    energy_profile.html
    reactant_view.html
    ts_view.html
    product_view.html
    path.xyz
    path_animation.html
  networks/
    reaction_network.dot
    reaction_network.svg
    reaction_network.html
    node_images/*.svg
  chemiscope/
    hf_screening.json.gz
```

## Renderer choices

- 3D structures and path animation: 3Dmol.js HTML snippets with XYZ fallback.
- Energy diagrams and screening maps: Plotly if installed; static/table fallback otherwise.
- 2D molecule depictions: RDKit MolDraw2D if available; SVG placeholder fallback otherwise.
- Reaction network: Graphviz DOT/SVG if `dot` is available; DOT-in-SVG fallback otherwise.
- Structure-property exploration: Chemiscope-like JSON gzip export.

## Scientific review notes

The report emphasizes `quality_tier`, `confidence_score`, reaction energetics, public DB support, and fallback/dummy warnings. Offline dummy data should be used to validate workflow and visualization logic, not as final DFT evidence.

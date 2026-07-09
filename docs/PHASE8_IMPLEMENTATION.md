# Phase 8 implementation: independent visualization layer

Phase 8 adds an independent `hfauto_viz` package for visualization and review.
It consumes an existing hfauto run directory and writes HTML/SVG/JSON outputs
under `14_viz/` without modifying calculation artifacts.

## Implemented visualization outputs

- `14_viz/report.html`: run-level visualization index
- `14_viz/screening/energy_profile_top_candidates.html`: Plotly energy profiles
- `14_viz/screening/assoc_vs_barrier.html`: ΔGassoc(p) vs ΔG‡ map
- `14_viz/screening/activation_map.html`: HF activation descriptor map
- `14_viz/screening/scavenger_map.html`: HF scavenger/capture map
- `14_viz/screening/cluster_risk_map.html`: capture vs cluster-risk proxy
- `14_viz/screening/quality_funnel.html`: workflow quality funnel
- `14_viz/screening/failure_heatmap.html`: failure category/stage heatmap
- `14_viz/molecules/<mol_id>/dossier.html`: molecule dossier with RDKit 2D and 3Dmol.js structure viewer
- `14_viz/reactions/<reaction_id>/dossier.html`: reaction dossier with energy profile, 3D endpoints and path animation
- `14_viz/networks/reaction_network.dot|svg|html`: Graphviz reaction network
- `14_viz/chemiscope/hf_screening.json.gz`: Chemiscope-like structure/property bundle
- `14_viz/viz_manifest.json`: manifest of generated visualization artifacts

## Design rules

1. `hfauto` calculation stages do not import `hfauto_viz`.
2. `hfauto_viz` only reads manifests, CSV/JSONL tables, XYZ/SDF structures and ranking outputs.
3. Visualization outputs are recorded separately in `viz_manifest.json`.
4. Dummy/fallback calculation data are explicitly warned in HTML reports.
5. The first implementation uses Plotly, 3Dmol.js HTML snippets, RDKit SVG and Graphviz.
6. The visualization package keeps compatibility aliases (`write_run_report`, `render_run_visual_report`, `energy_points(mode="process")`) for downstream scripts.

## CLI examples

```bash
hfauto-viz report runs/r001 --out runs/r001/14_viz --top-n 10
hfauto-viz molecule runs/r001 mol00012
hfauto-viz reaction runs/r001 rxn_mol00012_siteN3_hf2_PT
hfauto-viz network runs/r001
hfauto-viz chemiscope runs/r001
hfauto-viz status runs/r001
```

## Remaining visualization work

- Native Chemiscope export when optional `chemiscope` and `ase` are installed
- Cytoscape.js interactive network with right-panel 3D structure viewer
- Local/offline JS bundle manager for 3Dmol.js and Plotly instead of CDN mode
- Static high-quality PNG/SVG snapshots for paper-style figures
- TS imaginary-mode displacement arrows
- Full NEB/IRC energy-path animation using real path frames
- Dossier-level retry recommendations linked to failure artifacts
- Method/backend comparison dashboard

# Phase 8.5 implementation: visualization production hardening

Phase 8 introduced the independent `hfauto_viz` package. Phase 8.5 hardens it for review and production-style use without changing the scientific calculation artifacts.

## Added capabilities

- Offline/local asset bundle scaffold under `15_viz/assets/`.
- Asset modes: `cdn`, `local`, `auto`, `none`.
- Interactive Cytoscape-style reaction network with Graphviz SVG fallback.
- Graphviz nodes can carry molecule thumbnails and dossier links.
- Candidate ranking table with RDKit 2D thumbnails and molecule dossier links.
- Reaction dossiers now include:
  - interactive energy diagram,
  - static SVG/CSV energy diagram,
  - 3D reactant/TS/product viewers,
  - smoothed reaction path animation,
  - reaction-coordinate plot and CSV.
- Reaction-coordinate diagnostics compute `r(B-H)`, `r(H-F)`, and `q = r(B-H)-r(H-F)` using ReactionRecord atom indices when available, or a conservative B-H-F auto-detection fallback for visualization.
- Chemiscope-like export now writes a companion HTML index.
- `viz_manifest.json` records asset bundle, ranking thumbnails, Cytoscape network, coordinate plots, and static energy outputs.

## Scope and non-goals

`hfauto_viz` remains a read-only consumer of an hfauto run directory. It never changes energies, QC flags, ranking scores, or calculation artifacts. Placeholder local JavaScript files are intentionally small. For a fully offline interactive report, replace the placeholder files in `15_viz/assets/` with approved vendor copies of Plotly, 3Dmol.js, and Cytoscape.js.

## Example

```bash
PYTHONPATH=. python -m hfauto.cli.main pipeline \
  --config configs/pipelines/phase85_viz_hardened_offline.yaml \
  --run-id phase85_demo

# or generate visualization only for an existing run
PYTHONPATH=. python -m hfauto_viz.cli.main report runs/phase85_demo \
  --out runs/phase85_demo/15_viz \
  --library-mode local
```

## Main outputs

```text
15_viz/report.html
15_viz/viz_manifest.json
15_viz/assets/asset_manifest.json
15_viz/screening/candidate_ranking_thumbnails.html
15_viz/networks/reaction_network_interactive.html
15_viz/networks/reaction_network_interactive.json
15_viz/reactions/<reaction_id>/reaction_coordinate.html
15_viz/reactions/<reaction_id>/reaction_coordinate.csv
15_viz/reactions/<reaction_id>/energy_profile.svg
15_viz/reactions/<reaction_id>/energy_profile.csv
15_viz/chemiscope/hf_screening.html
```

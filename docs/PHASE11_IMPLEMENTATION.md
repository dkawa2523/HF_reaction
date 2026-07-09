# Phase 11: Production visualization refinement

Phase 11 strengthens `hfauto_viz` without changing the scientific calculation artifacts.
The visualization package remains a read-only consumer of completed `hfauto` runs.

## Added capabilities

- Annotated B-H-F structure viewers for reactant, TS, and product.
- Bond-change tables for forming B-H and breaking H-F distances.
- Approximate imaginary-mode viewer for offline/fallback runs, with data files that can be replaced by true normal-mode vectors from production ORCA parsing.
- Path-source detection that prefers real IRC/NEB/scan XYZ artifacts and falls back to interpolated reactant → TS → product frames.
- Linked Cytoscape-style reaction network and 3Dmol viewer.
- Publication/review figure bundle index for each reaction dossier.
- Phase 11 pipeline configuration.

## Important scientific note

If true Hessian-derived normal-mode vectors are not available, the imaginary-mode viewer uses a chemically interpretable approximation based on the canonical proton-transfer coordinate q = r(B-H) - r(H-F). The generated HTML and CSV explicitly mark this as an approximation and must not be treated as a substitute for verified frequency/IRC analysis.

## Key outputs

```text
15_viz/reactions/<reaction_id>/bond_change_annotations.csv
15_viz/reactions/<reaction_id>/imaginary_mode.html
15_viz/reactions/<reaction_id>/imaginary_mode_vectors.csv
15_viz/reactions/<reaction_id>/path_source.json
15_viz/reactions/<reaction_id>/publication_figures.html
15_viz/networks/reaction_network_linked.html
```

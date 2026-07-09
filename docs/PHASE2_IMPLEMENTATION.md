# Phase 2 implementation notes

This phase turns the initial skeleton into a more useful HF gas-phase reactivity screening prototype while keeping the same design rule: every stage reads a `manifest.json` and writes a new `manifest.json`.

## Implemented in this phase

### Scientific output improvements

- Added explicit candidate and standalone `(HF)n` species in `build-hf` so association energies can be computed from separated components.
- Added reaction endpoint generation for `B...(HF)n -> BH+...F(HF)n-1-` with atom-order consistency for NEB/IRC workflows.
- Added HF activation descriptors: `r_HF_A`, `delta_r_HF_A`, `nu_HF_cm1`, `delta_nu_HF_cm1`, `B_H_distance_A`, and `B_H_F_angle_deg` where available.
- Added pressure-aware thermochemistry fields in `thermo_records.jsonl` and `reaction_thermo.csv`.
- Added quality tiers and confidence scores in final thermochemistry and ranking outputs.
- Added separate output tables for:
  - `rank_scavenger.csv`
  - `rank_activation.csv`
  - `cluster_risk.csv`
  - `candidate_summary.csv`
  - `reaction_results.csv` or `reaction_results.parquet`
  - `failure_report.csv`

### Architecture and usability improvements

- Added partial pipeline options: `--from`, `--to`, and `--start-manifest`.
- Kept all stages independently executable through `hfauto run-stage`.
- Added `hfauto failures` and `hfauto inspect` for result review.
- Expanded HTML reporting with artifact counts, quality tier counts, failure summary, ranking previews, and candidate summary.
- Added a cache-first PubChem provider skeleton with network disabled by default.

### Package design improvements

- Added typed record contracts under `hfauto/core/schemas/records.py` while retaining generic `Artifact` transport.
- Kept chemistry-specific logic under `hfauto/chemistry/`.
- Kept execution-specific logic under `hfauto/backends/`.
- Kept stage orchestration under `hfauto/stages/`.
- Added deterministic dummy QM behavior so full DAG tests produce meaningful relative signs for association, ion-pair formation, and TS barriers.

## Validation commands

```bash
PYTHONPATH=. pytest -q
PYTHONPATH=. python -m hfauto.cli.main pipeline --config configs/pipelines/dummy.yaml --run-id demo_phase2
PYTHONPATH=. python -m hfauto.cli.main status runs/demo_phase2
PYTHONPATH=. python -m hfauto.cli.main report runs/demo_phase2 --out runs/demo_phase2/report.html
```

Expected current status for the bundled example is zero failures and populated ranking tables.

## Remaining production work

- Replace dummy QM with real xTB/CREST execution and parser.
- Replace dummy DFT with ORCA/Psi4/PySCF/NWChem execution and parser.
- Replace dummy TS midpoint with ORCA NEB-TS, constrained scan -> OptTS, and pysisyphus fallback.
- Replace dummy IRC with real IRC endpoint graph matching.
- Replace dummy thermal model with GoodVibes/ORCA Quasi-RRHO/Arkane.
- Add molecular dossier figures and energy diagrams.

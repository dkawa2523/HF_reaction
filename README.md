# hfauto finalized baseline

`hfauto` is a modular workflow for gas-phase HF reactivity screening of SDF candidate molecules.  It focuses on HF association, proton transfer, HF-cluster assisted ion-pair formation, thermochemistry, kinetics, ranking, visualization, and production operations.

This finalized baseline keeps the Phase 12.5 safety improvements and adds documentation polish.  It does **not** claim that offline or dummy/fallback results are chemically valid.  Use offline runs to verify wiring, reports, and operations outputs only.

## Package boundaries

```text
hfauto      calculation workflow, public-data enrichment, thermochemistry, kinetics, ranking
hfauto_viz  read-only visualization and molecule/reaction dossiers
hfauto_ops  read-only HPC planning, retry, reuse, QCArchive payloads, backend comparison
```

`hfauto_viz` and `hfauto_ops` read artifacts from a run directory.  They do not rewrite energies, structures, rankings, or QC results.

## Quick start: offline smoke run

```bash
PYTHONPATH=. python -m hfauto.cli.main pipeline \
  --config configs/pipelines/final_offline_smoke.yaml \
  --run-id demo_final

PYTHONPATH=. python -m hfauto.cli.main science-status runs/demo_final
PYTHONPATH=. python -m hfauto.cli.main next-actions runs/demo_final
```

Main review entry points:

```text
runs/<run_id>/15_viz/report.html
runs/<run_id>/14_rank/candidate_summary.csv
runs/<run_id>/14_rank/rank_screening.csv
runs/<run_id>/14_rank/rank_production.csv
runs/<run_id>/14_rank/ranking_summary.json
runs/<run_id>/13_connector-audit/production_connector_audit.html
runs/<run_id>/16_ops/operations_report.html
```

## Which ranking should I use?

```text
rank_screening.csv    development/exploration; may include fallback values
rank_scientific.csv   rows with real DFT minima or better
rank_production.csv   validated TS/IRC plus production thermochemistry
rank_scavenger.csv    HF scavenger-oriented score, with eligibility flags
rank_activation.csv   HF activation-oriented score, with eligibility flags
candidate_summary.csv one row per candidate with science/process/ops next actions
```

Production decisions should only use rows where:

```text
production_rank_eligible = true
main_values_are_dummy = false
main_values_are_fallback = false
quality_tier >= Q4
production_thermo_ready = true
real_irc_executed = true
```

## Independent stage execution

Every stage reads a `manifest.json` and writes a new `manifest.json`.

```bash
PYTHONPATH=. python -m hfauto.cli.main run-stage thermo \
  --in runs/<run_id>/09_sp/manifest.json \
  --out runs/<run_id>/10_thermo_rerun \
  --config configs/stages/thermo_goodvibes_phase10.yaml \
  --global-config configs/pipelines/final_offline_smoke.yaml
```

The `--global-config` option passes temperature, pressure, and mode settings into independent stage runs.

## Quality tiers

```text
Q0  structure/data only
Q1  screening or fallback values; no real DFT evidence
Q2  real DFT minima/frequency available
Q3  TS frequency available, IRC not validated
Q4  TS + IRC validated
Q5  Q4 + high-level single point / calibration support
```

Confidence scores are capped by evidence tier so Q0/Q1 dummy or fallback values cannot look production-ready.

## Documentation map

Start with:

```text
docs/current/index.md                 documentation entry point
docs/current/hf_gas_reactivity_technical_report.md  technical report for chemists and reaction/process engineers
docs/current/user_guide.md            user workflows and outputs
docs/current/science_gate.md          scientific readiness and ranking gates
docs/current/stage_io_reference.md    compact stage input/output map
docs/current/output_reference.md      result tables and how to interpret them
docs/current/extension_points.md      backend and stage extension points
docs/current/future_methods.md        useful future chemistry/simulation methods
docs/current/stage_future_extensions.md stage-by-stage future extensions
docs/current/production_checklist.md  before using results for candidate decisions
docs/current/developer_architecture.md package boundaries and code organization
```

Historical phase notes are kept under `docs/PHASE*.md` for traceability, but day-to-day users should read `docs/current/` first.

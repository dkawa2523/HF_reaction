# User guide

## 1. Offline smoke run

Use this to confirm installation, workflow wiring, reports, ranking tables, and operations outputs.  The numbers are not DFT evidence.

```bash
PYTHONPATH=. python -m hfauto.cli.main pipeline \
  --config configs/pipelines/final_offline_smoke.yaml \
  --run-id demo_final
```

Then inspect:

```bash
PYTHONPATH=. python -m hfauto.cli.main science-status runs/demo_final
PYTHONPATH=. python -m hfauto.cli.main next-actions runs/demo_final
```

## 2. Main files to open

```text
15_viz/report.html                              visual review entry point
14_rank/candidate_summary.csv                   one-row-per-candidate summary
14_rank/ranking_summary.json                    screening/scientific/production row counts
14_rank/rank_screening.csv                      exploratory ranking
14_rank/rank_production.csv                     production candidate ranking, often empty in offline runs
13_connector-audit/production_connector_audit.html connector readiness
16_ops/operations_report.html                   retry, reuse, and HPC planning
```

## 3. How to read candidate summary

The key columns are:

```text
quality_tier                  evidence level, Q0 to Q5
confidence_score              tier-capped confidence, not a probability
scientific_rank_eligible       true only when real DFT minima evidence is present
production_rank_eligible       true only when TS/IRC and production thermo gates pass
main_values_are_dummy          true when the main numeric values are development placeholders
main_values_are_fallback       true when fallback values are used
science_next_action            next chemistry/simulation step
process_next_action            next process or EHS review step
ops_next_action                next execution/operations step
```

## 4. Independent stage execution

Stages can be rerun independently with a prior manifest.

```bash
PYTHONPATH=. python -m hfauto.cli.main run-stage dft-minima \
  --in runs/<run_id>/05_preopt/manifest.json \
  --out runs/<run_id>/06_dft_minima_rerun \
  --config configs/stages/dft_minima_orca_offline.yaml \
  --global-config configs/pipelines/final_offline_smoke.yaml
```

`--global-config` keeps temperature, pressure, and run mode consistent with the original campaign.

## 5. Recommended production workflow

1. Run an offline smoke campaign.
2. Replace fallback preopt with real xTB/CREST.
3. Run real DFT minima/frequency for selected species.
4. Run real TS/IRC for selected reactions.
5. Run GoodVibes or ORCA quasi-RRHO thermal corrections.
6. Run high-level single-point corrections for final candidates.
7. Use `rank_production.csv` and reaction dossiers for decision review.

# Phase 9 Implementation: HPC / Large-scale operations

Phase 9 adds a read-only operations layer for scaling HF gas-phase reactivity screening from a local research workflow to campaign-scale execution.

## Scope

Phase 9 does **not** modify scientific artifacts. It generates operational artifacts from existing manifests and pipeline YAML:

- artifact and stage indexes
- resource plans
- SLURM submit templates
- retry plans
- duplicate/reuse candidate indexes
- calculation cache indexes
- run comparison reports
- backend comparison tables
- Snakemake profile / Snakefile templates
- operations dashboard HTML

## Main packages

```text
hfauto_ops/
  core/run_index.py       Artifact/stage index, duplicate signatures, SQLite export
  core/resource.py        Stage resource estimator and resource_plan.csv
  core/retry.py           Failure/fallback retry plan
  core/compare.py         Run comparison tables and HTML
  schedulers/slurm.py     SLURM submit script bundle
  reports/operations_report.py
  bundle.py               One-call operations bundle

hfauto/hpc/
  resources.py            Lightweight resource contract and estimates
  job_plan.py             JobRecord, job_plan.csv/json, retry job records
  schedulers.py           SLURM/PBS/local submit script templates, Snakemake profile
  cache.py                Calculation cache-key and duplicate detection
  compare.py              Generic manifest-level comparison utilities
  dashboard.py            HPC dashboard HTML
```

`hfauto_ops` is the recommended user-facing Phase 9 operations package. `hfauto/hpc` contains reusable lower-level utilities used by the optional `hpc-plan` stage and CLI helpers.

## New stages

### `ops`

Read-only operations bundle stage.

Inputs:

- latest workflow manifest
- optional `pipeline_config`
- optional scheduler/resource settings

Outputs:

```text
artifact_index.csv/jsonl
stage_index.csv
duplicate_candidates.csv
ops_index.sqlite
resource_plan.csv/json
retry_plan.csv/json
retry_commands.sh
slurm/*.slurm
slurm/submit_all.sh
operations_report.html
operations_manifest.json
```

Artifact:

```text
artifact_type = operations_bundle
```

### `hpc-plan`

Lower-level execution-plan stage that produces:

```text
job_plan.csv/json
retry_plan.csv/json
calculation_cache_index.csv/json
snakemake/Snakefile
snakemake_profile/config.yaml
submit/*.slurm.sh
hpc_dashboard.html
```

Artifact types:

```text
hpc_plan
cache_index
hpc_bundle
```

## CLI examples

```bash
# Full offline Phase 9 demo pipeline
PYTHONPATH=. python -m hfauto.cli.main pipeline \
  --config configs/pipelines/phase9_hpc_ops_offline.yaml \
  --run-id phase9_demo

# Build operation bundle on an existing run
PYTHONPATH=. python -m hfauto_ops.cli.main bundle runs/phase9_demo \
  --out runs/phase9_demo/16_ops \
  --pipeline-config configs/pipelines/phase9_hpc_ops_offline.yaml

# Build generic HPC plan on an existing run
PYTHONPATH=. python -m hfauto.cli.main hpc-plan runs/phase9_demo \
  --out runs/phase9_demo/ops_hpc \
  --config configs/pipelines/phase9_hpc_ops_offline.yaml \
  --scheduler slurm

# Compare two runs
PYTHONPATH=. python -m hfauto_ops.cli.main compare-runs runs/run_a runs/run_b \
  --out comparisons/run_a_vs_run_b
```

## Review notes

- Generated submit scripts are templates. They do not load site-specific modules, activate environments, or set accounts.
- Retry plans include failed artifacts and dummy/fallback artifacts because fallback results are useful for software validation but weak as production scientific evidence.
- Duplicate signatures are conservative reuse candidates, not automatic substitution rules.
- Operations artifacts are deliberately independent of scientific scoring and do not change rankings.

## Remaining production work

- Site-specific SLURM/PBS profiles and module/environment templates
- Snakemake executor integration tests on an actual HPC scheduler
- QCArchive/QCFractal connector
- robust array-job splitting by molecule/site/reaction
- centralized calculation cache shared across runs
- resource-estimate calibration from historical walltimes

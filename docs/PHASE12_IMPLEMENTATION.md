# Phase 12: HPC production integration

Phase 12 adds the production operations layer for running HF reactivity campaigns at scale.  The scientific workflow remains unchanged; this phase adds a read-only operations package that indexes completed runs and produces reviewable execution, reuse, archive, and comparison plans.

## Scope

Added capabilities:

- stage-level resource plans
- artifact-level job-array plans for heavy calculations
- SLURM/PBS/LSF/local scheduler templates
- scheduler submit/status/cancel command wrappers, dry-run by default
- retry plans for failed or fallback artifacts
- conservative calculation reuse planning
- QCArchive/QCFractal export payloads and import plan templates
- backend/method comparison tables
- resource auto-tuning suggestions
- operations HTML report and manifest

## Scientific intent

Phase 12 does not change energies, geometries, thermochemistry, kinetics, rankings, or visualizations.  It helps campaign operators answer:

- which HF complexes, ion pairs, TS searches, IRCs, and SP calculations should be run as separate HPC jobs;
- which dummy/fallback results need production replacement;
- which prior calculations are possible reuse candidates;
- which backends or methods differ for the same species/reaction;
- which stages are likely to need larger walltime or memory;
- what should be exported to QCArchive after schema review.

## Main package boundaries

```text
hfauto      : chemistry, QM, thermochemistry, kinetics, rankings
hfauto_viz  : reports, 3D/2D visualization, reaction dossiers
hfauto_ops  : HPC, retry, reuse, archive, run/backend comparison
```

`hfauto_ops` reads run artifacts and writes operations artifacts. It is a read-only consumer with respect to scientific results.

## Key outputs

```text
16_ops/
  artifact_index.csv
  stage_index.csv
  resource_plan.csv
  array_job_plan.csv
  array_items.jsonl
  artifact_job_plan.csv
  artifact_job_arrays.csv
  scheduler/
    slurm/ submit_all.sh, submit_array_slurm.sh, status.sh, cancel.sh
    pbs/   submit_all.sh, submit_array_pbs.sh, status.sh, cancel.sh
    lsf/   submit_all.sh, submit_array_lsf.sh, status.sh, cancel.sh
  retry_plan.csv
  reuse_plan.csv
  calculation_reuse_review.csv
  qcarchive_submit_payload.jsonl
  qcarchive_summary.json
  backend_comparison.csv
  backend_comparison_dashboard.html
  resource_autotune_recommendations.csv
  operations_report.html
  operations_manifest.json
```

## Example commands

Full offline pipeline:

```bash
PYTHONPATH=. python -m hfauto.cli.main pipeline \
  --config configs/pipelines/phase12_hpc_production_offline.yaml \
  --run-id phase12_demo
```

Build only operations bundle for an existing run:

```bash
PYTHONPATH=. python -m hfauto_ops.cli.main bundle runs/phase12_demo \
  --pipeline-config configs/pipelines/phase12_hpc_production_offline.yaml \
  --out runs/phase12_demo/16_ops
```

Generate artifact-level array plan:

```bash
PYTHONPATH=. python -m hfauto_ops.cli.main artifact-plan runs/phase12_demo \
  --pipeline-config configs/pipelines/phase12_hpc_production_offline.yaml \
  --out runs/phase12_demo/16_ops
```

Generate portable scheduler scripts:

```bash
PYTHONPATH=. python -m hfauto_ops.cli.main scheduler-bundle runs/phase12_demo \
  --pipeline-config configs/pipelines/phase12_hpc_production_offline.yaml \
  --scheduler slurm \
  --out runs/phase12_demo/16_ops
```

Create reuse and QCArchive plans:

```bash
PYTHONPATH=. python -m hfauto_ops.cli.main reuse-plan runs/phase12_demo \
  --out runs/phase12_demo/16_ops

PYTHONPATH=. python -m hfauto_ops.cli.main qcarchive-plan \
  runs/phase12_demo/16_ops/artifact_job_plan.csv \
  --out runs/phase12_demo/16_ops/qcarchive
```

## Production notes

- Submit/status/cancel commands are dry-run/review-first by default.
- Reuse plans are conservative. They should not substitute calculations automatically without geometry, charge/multiplicity, method, basis, and QC review.
- QCArchive export is a bridge format. Native QCSchema/QCFractal submission requires additional site-specific schema and credential validation.
- Artifact-level array commands rely on sliced manifests and are intended for production hardening on the target HPC environment before large campaigns.

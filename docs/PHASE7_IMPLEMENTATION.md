# Phase 7 implementation: public DB enrichment and method-validation scaffolding

Phase 7 adds a public-data and calibration layer on top of the Phase 6 HF gas-phase
reactivity workflow.  The design goal is to make computed HF association,
ion-pair, TS, thermo, and kinetics results easier to review scientifically while
preserving independent stage execution and offline reproducibility.

## Added capabilities

### Public DB providers

The following providers are wired into `hfauto enrich` through the common provider
registry:

- `pubchem`: CID, identity cross-check, basic molecular properties, synonyms, and
  optional 3D SDF seed through PubChem PUG-REST. Network access is opt-in and
  cache-first.
- `nist_webbook`: PA/GB and spectrum-oriented reference data. Online WebBook
  parsing is conservative because WebBook content is HTML; fixture and bundled
  reference-set matches are supported for reproducible QA.
- `atct`: thermochemical reference fields from fixture/local tables.
- `cccbdb`: vibrational/dipole/geometry reference fields from fixture/local
  tables.
- `comptox`: process/EHS property placeholder with fixture and opt-in CTX API
  skeleton.
- `niosh`: EHS review flags from fixture or bundled flags.
- `cas_common_chemistry`: optional CAS RN candidate extraction/cross-check.

All providers return machine-readable status under `db_provider_status` and never
stop the workflow because of cache misses or unavailable public records.

### Public-data audit outputs

`hfauto enrich` now writes:

- `enriched_molecule_records.jsonl`
- `public_data_coverage.csv`

Each `molecule_enriched` artifact contains:

- `identity` with PubChem/CAS/identity confidence flags;
- `public_data.basic_props` for molecular properties;
- `public_data.reference_values` for PA/GB, thermochemistry, IR/frequency, etc.;
- `public_data.gas_process` for gas-process feasibility and EHS review flags;
- `public_data.db_support` for public DB support scores.

### Method-validation / calibration stage

A new independent stage is available:

```bash
hfauto run-stage calibrate \
  --in runs/<run_id>/12_kinetics/manifest.json \
  --out runs/<run_id>/13_calibrate
```

It produces:

- `public_data_coverage.csv`
- `public_db_audit.csv`
- `calibration_reference_records.csv/jsonl`
- `method_validation_metrics.csv`
- `method_validation_report.html`
- `calibration_summary.jsonl`

The metrics include a basicity-rank diagnostic (`basicity_rank_spearman_proxy`) that
compares public PA values with computed HF association strength proxy when both are
available.  This is a diagnostic only; PA/GB values are not substituted for HF
association or TS free energies.

### Ranking integration

`rank` tables carry forward public-data and calibration columns where available:

- `reference_support_level`
- `calibration_action`
- `gas_process_feasibility`
- `ehs_review_flag`
- `public_db_support_score`

Process/EHS flags are treated as review penalties and gates, not as scientific
energy corrections.

## Pipeline configs

- `configs/pipelines/phase7_publicdb_calibration_offline.yaml`: fully offline
  review pipeline using fixtures/bundled references.
- `configs/pipelines/phase7_publicdb_calibration_network.yaml`: network-enabled
  production-style skeleton.  Review API/license terms before use.
- `configs/pipelines/phase7_full_publicdb_rank_offline.yaml`: offline pipeline
  with rank outputs and public-data audit.

## Validation

Validated with:

```bash
PYTHONPATH=. pytest -q
```

Result: `21 passed`.

Offline pipeline smoke test:

```bash
PYTHONPATH=. python -m hfauto.cli.main pipeline \
  --config configs/pipelines/phase7_publicdb_calibration_offline.yaml \
  --run-id phase7_verify
```

The latest manifest stage is `calibrate`; ranking artifacts are carried forward
from the previous stage.

## Remaining limitations

- Production public DB use should pin API versions, license/terms, and curated
  snapshots.
- NIST WebBook parsing remains conservative because the source is HTML.
- ATcT/CCCBDB are currently fixture/local-table oriented; production API/client
  integration is a future enhancement.
- CompTox CTX production use needs an EPA API key and endpoint/version config.
- Calibration metrics are coverage/diagnostic metrics until real high-level QM
  calculations replace dummy/fallback values.

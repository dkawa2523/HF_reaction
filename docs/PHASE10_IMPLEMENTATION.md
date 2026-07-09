# Phase 10: Production connector hardening

Phase 10 strengthens the external connector layer without changing the core
Artifact/Manifest contracts.

## Scope

Implemented:

- GoodVibes external execution adapter and tolerant CSV parser.
- GoodVibes per-species production-readiness flags.
- Arkane skeleton export plus optional opt-in external execution wrapper.
- Cantera draft mechanism validation, optional real Cantera smoke-test reactor,
  and connector QC JSON outputs.
- `connector-audit` stage that summarizes GoodVibes, Arkane, Cantera, and public
  DB provider readiness.
- Calibration dashboard additions: reference residual/proxy table and dashboard
  HTML.
- Phase 10 offline and production-style connector templates.

Not implemented as guaranteed production science:

- GoodVibes/Arkane/Cantera scientific correctness for a specific project; users
  must supply real ORCA logs, reviewed Arkane species/TS blocks, NASA thermo, and
  project-specific mechanism definitions.
- Live public DB stability guarantees. Providers remain cache-first and
  network-opt-in.

## Key output files

```text
10_thermo/species_thermo.csv
10_thermo/reaction_thermo.csv
12_kinetics/cantera_validation.json
12_kinetics/cantera_reactor_results.csv
12_kinetics/arkane/arkane_export_summary.json
13_calibrate/reference_residuals.csv
13_calibrate/method_calibration_dashboard.html
14_connector-audit/production_connector_audit.csv
14_connector-audit/production_connector_summary.csv
14_connector-audit/production_connector_audit.html
```

## Production flags

Important record fields:

```text
production_thermo_ready
external_goodvibes_executed
external_goodvibes_status
production_cantera_ready
cantera_validation_status
production_arkane_ready
arkane_status
connector_quality
```

These flags are deliberately conservative. A connector can run successfully but
still not be scientifically production-ready if the mechanism/species definitions
are placeholders.

## Example

Offline regression:

```bash
PYTHONPATH=. python -m hfauto.cli.main pipeline \
  --config configs/pipelines/phase10_production_connectors_offline.yaml \
  --run-id phase10_offline
```

Production-style template:

```bash
HFAUTO_ALLOW_GOODVIBES=1 HFAUTO_ALLOW_CANTERA=1 PYTHONPATH=. python -m hfauto.cli.main pipeline \
  --config configs/pipelines/phase10_real_connectors_template.yaml \
  --run-id phase10_real
```

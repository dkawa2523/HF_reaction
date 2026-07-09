# Code maintenance guide

## Keep the package split clean

```text
hfauto      may create and modify scientific artifacts
hfauto_viz  may only read scientific artifacts and write visualization outputs
hfauto_ops  may only read scientific artifacts and write operations outputs
```

## Prefer small adapters

When adding a backend, keep the stage interface stable and put tool-specific logic inside the backend module.  Avoid adding tool-specific columns unless they are useful for QC or interpretation.

## Do not over-expand contracts

The current code uses lightweight conventions:

```text
manifest.json in, manifest.json out
artifact_type string
common data fields where possible
science gates in ranking outputs
```

Add runtime enforcement only when it prevents a real class of errors.  Prefer documentation and small helper constants over heavy schema layers unless production incidents justify stricter contracts.

## When adding new outputs

Add:

```text
1. a human-readable CSV/HTML when useful
2. a machine-readable JSON/JSONL when downstream code needs it
3. a short note in docs/current/output_reference.md
4. one focused test
```

Avoid adding multiple overlapping diagnostics for the same question.

## Testing priorities

1. Offline smoke pipeline still completes.
2. Dummy/fallback values do not pass production gates.
3. Parsers have golden-output tests for every production connector.
4. Ranking tables include science/process/ops next actions.
5. Visualization and operations remain read-only with respect to scientific artifacts.

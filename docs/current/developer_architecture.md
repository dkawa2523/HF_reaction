# Developer architecture

The code is split into three packages:

```text
hfauto      compute workflow and artifacts
hfauto_viz  visualization read-only consumer
hfauto_ops  operations/HPC read-only consumer
```

## Core rule

A stage receives:

```text
input manifest + stage config + global config -> output manifest + files
```

Previous stage outputs are not overwritten.

## Where code should go

```text
hfauto/core        shared units, hashing, manifest/artifact schemas, QC helpers
hfauto/chemistry   chemistry logic: site detection, HF templates, descriptors
hfauto/backends    external tools and replaceable implementations
hfauto/stages      stage orchestration and stage-level output files
hfauto_viz         report and structure visualization
hfauto_ops         HPC/retry/reuse/QCArchive planning
```

## Keep extensions simple

Prefer:

- one clear CSV/JSON output over several diagnostic fragments;
- failure artifacts over uncaught exceptions;
- explicit dummy/fallback flags over implicit assumptions;
- backend-local parsing logic over stage-specific parsing logic.

## Do not do this

- Do not let visualization or operations packages modify calculation artifacts.
- Do not let dummy/fallback data enter production rankings.
- Do not add new artifact types without documenting which stage creates and consumes them.
- Do not make public database live access mandatory for pipeline execution.

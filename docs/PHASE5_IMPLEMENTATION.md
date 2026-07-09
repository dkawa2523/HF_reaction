# Phase 5 Implementation: TS Search and IRC Validation

Phase 5 adds the transition-state and IRC layer on top of the Phase 4 DFT opt/freq/SP infrastructure.

## Scientific scope

The implemented reaction-path contract focuses on proton-transfer reactions of the form:

```text
B···(HF)n -> BH+···F(HF)n-1-
```

For each reaction seed, Phase 5 can now produce:

- `ts_path` artifacts containing ORCA NEB-TS input and endpoint files.
- `species` artifacts with `state = transition_state`.
- `calculation` artifacts for TS frequency / OptTS results.
- `reaction_validated` artifacts with TS frequency and mode QC.
- `irc_attempt` and `irc` artifacts with forward/backward endpoint checks.
- `reaction_path_validated` artifacts summarizing TS + IRC path status.

The quality tier logic remains conservative: dummy/fallback TS or IRC artifacts are useful for software review and pipeline validation, but they do not upgrade a result to real Q3/Q4/Q5 scientific quality.

## Implemented modules

```text
hfauto/backends/ts/base.py
  Backend-neutral result objects and endpoint/reaction-coordinate helpers.

hfauto/backends/ts/dummy.py
  Deterministic TS/IRC fallback for offline review and CI.

hfauto/backends/ts/orca_nebts.py
  ORCA NEB-TS, OptTS, and IRC input renderers; opt-in subprocess execution;
  ORCA output parsing; dummy fallback when explicitly configured.

hfauto/backends/ts/scan_optts.py
  Interface-compatible constrained-scan/OptTS placeholder.

hfauto/backends/ts/pysisyphus.py
  Interface-compatible pysisyphus placeholder.

hfauto/chemistry/reaction_path_qc.py
hfauto/chemistry/reaction_path.py
  Geometry and endpoint QC helpers.

hfauto/stages/ts_search.py
  Engine-order TS search stage.

hfauto/stages/irc.py
  IRC execution / endpoint validation stage.
```

## Offline review pipeline

```bash
PYTHONPATH=. python -m hfauto.cli.main pipeline \
  --config configs/pipelines/phase5_ts_irc_offline.yaml \
  --run-id phase5_offline
```

This renders ORCA NEB-TS / IRC inputs but does not execute ORCA. It uses explicit dummy fallback and labels outputs with:

```text
fallback_dummy: true
scientific_use: software_test_only_not_ts / software_test_only_not_irc
```

## Production-style pipeline

```bash
HFAUTO_ALLOW_ORCA=1 PYTHONPATH=. python -m hfauto.cli.main pipeline \
  --config configs/pipelines/phase5_real_ts_irc.yaml \
  --run-id phase5_real
```

Production-style configs disable dummy fallback for ORCA DFT/TS/IRC so failed chemistry remains visible as failure artifacts.

## Key QC fields

TS artifacts include:

```text
n_imag
imag_freq_cm1
mode_overlap_score
mode_overlap_method
imag_mode_matches_reaction_coordinate
ts_validated_by_frequency
reaction_coordinate_progress_score
reaction_coordinate_between_endpoints
```

IRC artifacts include:

```text
forward_ok
backward_ok
reactant_endpoint_match
product_endpoint_match
endpoint_graph_match / endpoint_pair QC
endpoint_rmsd_A
irc_validated
real_irc_executed
fallback_dummy
```

## Remaining limitations

- Normal-mode displacement projection is currently heuristic unless backend output provides an explicit marker; future work should parse ORCA normal mode vectors for true projection onto `r(B-H)-r(H-F)`.
- `scan_optts` and `pysisyphus` satisfy the backend interface but still use deterministic fallback behavior.
- ORCA IRC endpoint discovery is tolerant and version-agnostic; production use should be validated with representative ORCA 6.x outputs.
- GoodVibes / Arkane / Cantera integration is still Phase 6.
- Public DB production providers and method-calibration reports are still Phase 7.

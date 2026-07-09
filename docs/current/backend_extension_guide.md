# Backend extension guide

Backends are replaceable implementations behind stable stage contracts.  New backends should not change artifact schemas unless the stage output contract genuinely needs new data.

## Recommended pattern

1. Add the backend under `hfauto/backends/<category>/`.
2. Keep external executable calls inside the backend.
3. Store raw stdout/stderr and command metadata when subprocesses run.
4. Return existing artifact types whenever possible.
5. Add a parser test with a small golden output.
6. Add a fallback or clear failure artifact, not an uncaught exception.

## QM backend checklist

A QM backend should report:

```text
engine, engine_version, method_id, task
scf_converged, geometry_converged
real_qm_executed, fallback_dummy
electronic_energy_hartree
zpe_hartree, gibbs_298K_hartree when available
frequencies_cm1, n_imag when available
final_xyz
```

## TS backend checklist

A TS backend should report:

```text
reaction_id
TS species id
n_imag, imag_freq_cm1
mode_overlap_score when available
real_ts_search_executed
fallback_dummy
```

## Thermochemistry backend checklist

A thermochemistry backend should report:

```text
production_thermo_ready
low_frequency_count
quasi_rrho_applied or equivalent
external tool status
G_standard_hartree
G_process_hartree
```

## Public database backend checklist

A public-data backend should be cache-first and network opt-in.  Store:

```text
provider
status
matched
source identifier
retrieval mode: fixture/cache/live
retrieved_at when live
notes or parser warnings
```

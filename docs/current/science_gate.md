# Production science gate

The workflow produces development, screening, scientific, and production outputs.  To avoid misuse, candidate decisions should be made from rows where `production_rank_eligible` is true.

## Minimal gates

Species-level:

```text
real QM backend used
fallback_dummy = false
SCF converged
geometry converged
minima have n_imag = 0
```

TS-level:

```text
n_imag = 1
imag_freq_cm1 < -100
mode_overlap_score >= 0.7
```

IRC-level:

```text
real_irc_executed = true
irc_validated = true
reactant_endpoint_match = true
product_endpoint_match = true
```

Thermochemistry:

```text
production_thermo_ready = true
low-frequency treatment recorded
temperature and pressure correction recorded
```

Ranking:

```text
rank_screening.csv   exploratory only
rank_scientific.csv  real DFT minima or better
rank_production.csv  validated TS/IRC and production thermo
```

Q0/Q1 confidence is capped.  Dummy/fallback values can appear in screening tables but are excluded from scientific and production rankings.

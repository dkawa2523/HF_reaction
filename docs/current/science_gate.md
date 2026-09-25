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
n_imag = exactly 1 below imaginary_frequency_cutoff_cm1
imag_freq_cm1 < imaginary_frequency_cutoff_cm1
mode_overlap_score >= configured path/coordinate threshold
```

`imaginary_frequency_cutoff_cm1` is the single numerical-noise boundary used
by frequency parsing, minimum QC, and TS QC. Its default is `-1.0 cm-1` with a
strict `<` comparison, matching the established NWChem minimum parser: for
example, `-5 cm-1` is counted and `-0.8 cm-1` is ignored. This cutoff only
separates reported numerical noise from modes to inspect; it is not by itself
evidence for a chemical saddle. A literal Cartesian mode projection onto the
target reaction coordinate is still required, and a reaction-path claim still
requires real bidirectional IRC connectivity to the validated endpoints.
The current reference pipelines use 0.50 as the authorization threshold for
starting IRC; that value is not a universal proof of connectivity and IRC
remains decisive.

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

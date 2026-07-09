# Production science validation plan

This plan defines how to validate the workflow before using it for candidate selection.

## Reference molecules

Start with:

```text
HF
NH3
methylamine
dimethylamine
trimethylamine
pyridine
aniline
water
methanol
isopropanol
```

## Quantities to validate

```text
proton affinity / gas basicity trend
HF association free energy trend
r_HF and HF stretch shift in B···HF
minima frequency sanity
TS frequency and IRC endpoint match for small acid-base systems
thermal correction sensitivity
```

## Method ladder

```text
GFN2-xTB / CREST
r2SCAN-3c opt/freq
ωB97X-D4/def2-TZVPPD SP or opt/freq subset
high-level SP for selected references
Gaussian same-method comparison where available
```

## Outputs

```text
reference_campaign_summary.csv
method_bias_summary.csv
frequency_error_summary.csv
thermo_sensitivity.csv
production_science_validation_report.html
```

## Decision criteria

Do not claim production accuracy until:

- the reference campaign completes with real QM outputs;
- method-tier errors are summarized;
- TS/IRC validation works on known small systems;
- production ranking uses only production-eligible rows.

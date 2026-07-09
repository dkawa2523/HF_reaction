# Stage-by-stage future extensions

This file lists useful future extensions by processing stage.  It is intentionally practical: each item should either improve chemical reliability, user review, or production execution.

## ingest / enrich

Future additions:

```text
salt/solvate stripping policy for gas-process candidates
formal charge and multiplicity review report
PubChem 3D seed prioritization by InChIKey identity confidence
process-property enrichment for vapor pressure / boiling point / hazard flags
```

Do not block proprietary candidates when public DB lookup fails.  Record `db_hit=false` and continue.

## detect-sites

Future additions:

```text
basicity-site ML classifier trained on curated PA/GB data
site steric accessibility score
multi-site bridge candidates for HF chains
explicit exclusion flags for amide, pyrrole, nitro, sulfonamide, quaternary ammonium
```

Useful output:

```text
site_confidence
site_priority
site_exclusion_reason
```

## conformers / build-hf

Future additions:

```text
multi-direction HF placement
HF cluster chain/cyclic templates
B···(HF)n CREST ensemble refinement
Boltzmann-weighted complex ensemble
failed-template diversity report
```

Useful output:

```text
n_complex_conformers
lowest_complex_energy
complex_ensemble_free_energy
hf_template_id
```

## preopt

Future additions:

```text
xTB batch runner with array jobs
GFN1/GFN2 sensitivity comparison
automatic detection of unintended proton transfer
fragment distance / graph sanity report
```

Useful output:

```text
hf_dissociated
proton_transferred_unintentionally
fragment_graph_ok
selected_for_dft
```

## dft-minima / sp

Future additions:

```text
method ladder automation
r2SCAN-3c -> wB97X-D4 -> double-hybrid -> DLPNO-CCSD(T)
basis sensitivity on final candidates
program-backend comparison ORCA/Psi4/PySCF/NWChem
QCSchema-compatible result export
```

Useful output:

```text
method_tier
method_bias_estimate
backend_comparison_delta
high_level_correction_applied
```

## ts-search / irc

Future additions:

```text
ORCA NEB-TS production parser hardening
scan -> OptTS fallback
pysisyphus GSM fallback
true imaginary mode overlap from normal-mode vectors
real IRC/NEB frame ingestion for visualization
```

Useful output:

```text
mode_overlap_score
endpoint_graph_match
irc_endpoint_rmsd_A
path_source
```

## thermo / kinetics

Future additions:

```text
external GoodVibes production parser hardening
ORCA quasi-RRHO direct extraction
hindered rotor correction for selected cases
Arkane TST/RRKM/Master Equation
Eckart tunneling for final candidates
Cantera production mechanism with NASA thermo
```

Useful output:

```text
production_thermo_ready
thermal_model
frequency_scale_factor
k_TST
k_corrected
kinetics_uncertainty_band
```

## rank / calibrate

Future additions:

```text
reference campaign bias correction
method-tier uncertainty annotation
separate activation/scavenger/process feasibility decision views
rank stability under method and conformer perturbation
```

Useful output:

```text
rank_stability
method_uncertainty_kcal_mol
production_rank_eligible
science_next_action
```

## viz / ops

Future additions:

```text
true IRC/NEB animation
normal-mode arrows from Hessian parser
interactive Cytoscape + 3Dmol + energy-profile synchronization
HPC job history ingestion
resource actuals and auto-tuning
QCArchive live synchronization
```

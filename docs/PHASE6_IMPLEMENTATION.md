# Phase 6 Implementation: thermochemistry, pressure corrections, TST kinetics, and mechanism export

Phase 6 adds a production-shaped thermochemistry and kinetics layer while keeping the Phase 5 Artifact/Manifest contract intact.  It is designed so the workflow remains executable without external GoodVibes/Arkane/Cantera installations, while making the output fields and quality flags ready for those production backends.

## Scientific additions

### Species-level thermochemistry

The `thermo` stage now creates `species_thermo` artifacts and `species_thermo.csv` with:

- `G_standard_hartree`
- `G_process_hartree`
- `thermal_correction_source_hartree`
- `thermal_correction_used_hartree`
- `pressure_correction_hartree`
- `quasi_rrho_applied`
- `low_frequency_count`
- `quasi_rrho_correction_kcal_mol`
- `energy_source_calc_id`
- `thermal_source_calc_id`
- `is_scientific_energy`
- `is_scientific_thermal`

The internal low-frequency correction is intentionally transparent and conservative.  It is useful for workflow review and sensitivity screening, but final reporting should use GoodVibes or ORCA quasi-RRHO thermochemistry from validated frequency calculations.

### Reaction-level thermochemistry

The `thermo` stage now creates `reaction_thermo.csv` and `thermo` artifacts with:

- `delta_G_assoc_standard_kcal_mol`
- `delta_G_assoc_pressure_corrected_kcal_mol`
- `delta_G_assoc_process_kcal_mol`
- `delta_G_ionpair_kcal_mol`
- `delta_G_act_kcal_mol`
- `delta_G_act_process_kcal_mol`
- `K_assoc_standard`
- `K_assoc_process_adjusted`
- `K_ionpair`
- `imag_freq_cm1`
- `quality_tier`
- `confidence_score`

Association uses the separated candidate and `(HF)n` reference species:

```text
ΔG_assoc° = G°[B···(HF)n] - G°[B] - G°[(HF)n]
```

The pressure-adjusted association value applies ideal-gas `RT ln(p/p°)` corrections to configured reactant partial pressures.

### Kinetics and mechanism export

The `kinetics` stage now creates:

- `kinetics_records.csv`
- `kinetics_records.jsonl`
- `hfauto_kinetics.yaml`
- `hfauto_kinetics.json`
- `cantera_mechanism.yaml`
- `arkane/arkane_tst_input.py`

The default `tst` backend computes Eyring TST rates and optional Wigner tunneling correction from the TS imaginary frequency.  The Cantera output is a review/skeleton mechanism with pseudo species; production Cantera simulation still requires Arkane or curated NASA thermochemistry for all species.

## Architecture additions

- `hfauto/core/thermo_models.py`: small auditable thermochemistry and TST helper functions.
- `hfauto/core/thermochemistry.py`: compatibility helpers used by tests and legacy Phase 6 code paths.
- `hfauto/backends/thermo/goodvibes.py`: GoodVibes-compatible thermochemistry backend with internal fallback.
- `hfauto/backends/thermo/arkane.py`: Arkane input skeleton exporter.
- `hfauto/backends/kinetics/tst.py`: Eyring/Wigner TST backend.
- `hfauto/backends/kinetics/cantera.py`: hfauto kinetics YAML/JSON and Cantera draft export.
- `hfauto/stages/thermo.py`: species-level and reaction-level thermochemistry.
- `hfauto/stages/kinetics.py`: kinetics records and mechanism exports.

## New pipeline configs

- `configs/pipelines/phase6_thermo_kinetics_offline.yaml`
- `configs/pipelines/phase6_real_thermo_kinetics.yaml`
- `configs/stages/thermo_phase6.yaml`
- `configs/stages/kinetics_phase6.yaml`

## Run examples

```bash
PYTHONPATH=. python -m hfauto.cli.main pipeline \
  --config configs/pipelines/phase6_thermo_kinetics_offline.yaml \
  --run-id phase6_offline
```

```bash
PYTHONPATH=. python -m hfauto.cli.main run-stage thermo \
  --in runs/<run_id>/09_sp/manifest.json \
  --out runs/<run_id>/10_thermo \
  --config configs/stages/thermo_phase6.yaml
```

```bash
PYTHONPATH=. python -m hfauto.cli.main run-stage kinetics \
  --in runs/<run_id>/10_thermo/manifest.json \
  --out runs/<run_id>/12_kinetics \
  --config configs/stages/kinetics_phase6.yaml
```

## Validation

```bash
PYTHONPATH=. pytest -q
PYTHONPATH=. python -m hfauto.cli.main pipeline \
  --config configs/pipelines/phase6_thermo_kinetics_offline.yaml \
  --run-id validation_phase6
```

Expected output includes populated thermochemistry, kinetics, ranking, and report tables with zero failure artifacts for the bundled example.

## Remaining production work

- External GoodVibes execution and parser tests against real GoodVibes output.
- ORCA quasi-RRHO thermochemistry option.
- Arkane production execution and parsing of TST/RRKM outputs.
- True Cantera gas mechanism generation after Arkane/NASA thermo is available.
- Public database production providers and method calibration report.
- Structure-rich molecule dossier and energy-diagram visualization.

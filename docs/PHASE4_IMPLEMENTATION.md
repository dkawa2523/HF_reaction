# Phase 4 implementation: ORCA DFT opt/freq and single-point backend

## Purpose

Phase 4 moves the workflow from software-contract/dummy energies toward Gaussian-class DFT calculations by adding a production-oriented ORCA backend for:

- minima optimization and frequency calculations, especially `r2SCAN-3c Opt Freq`;
- high-level single-point calculations, especially `wB97X-D4/def2-TZVPPD`;
- ORCA input rendering, execution, parsing, QC, and optimized-geometry propagation.

The implementation remains safe for developer laptops and CI: ORCA subprocess execution is disabled by default and must be explicitly enabled.

## Scientific output added in Phase 4

Each ORCA calculation can produce a `calculation` Artifact containing:

```text
electronic_energy_hartree
zpe_hartree
enthalpy_298K_hartree
gibbs_298K_hartree
thermal_correction_gibbs_hartree
thermal_correction_enthalpy_hartree
frequencies_cm1
n_imag
lowest_freq_cm1
imag_freq_cm1
hf_stretch_cm1
dipole_D
scf_converged
geometry_converged
normal_termination
geometry_sane
hf_dissociated
proton_transferred_unintentionally
```

The `dft-minima` stage now creates `species_optimized` Artifacts. Downstream `sp` and `descriptors` prefer these optimized geometries, while preserving canonical species IDs so thermochemistry joins remain stable.

## Architecture decisions

### Execution is opt-in

ORCA can be invoked only when one of the following is true:

```text
settings.allow_subprocess: true
HFAUTO_ALLOW_ORCA=1
HFAUTO_ALLOW_SUBPROCESS=1
```

If ORCA is not run, the backend either returns a recoverable `orca_not_run` failure Artifact or, if `fallback_to_dummy: true`, returns a deterministic dummy calculation marked with:

```text
qc.fallback_dummy = true
qc.real_orca_executed = false
qc.scientific_use = software_test_only_not_dft
```

### Stage code is backend-agnostic

`dft-minima` and `sp` call the QMEngine interface only:

```python
engine.optimize_frequency(species, method, workdir)
engine.single_point(species, method, workdir)
```

ORCA-specific input syntax, subprocess behavior, parser details, and fallback behavior are isolated in `hfauto/backends/qm/orca.py`.

### Failures are first-class artifacts

If ORCA input rendering, execution, or parsing fails, the result is a structured failed Artifact with category, reason, recommended fallback, input/output paths, and command-result metadata.

## New/updated files

```text
hfauto/backends/qm/orca.py
  ORCA input renderer, subprocess runner, parser, dummy fallback, optimized XYZ extraction

hfauto/stages/dft_minima.py
  creates species_optimized artifacts and optimized_species_records.jsonl

hfauto/stages/sp.py
  prefers species_optimized, then species_preopt, then original species

hfauto/stages/descriptors.py
  prefers optimized/preoptimized reactant-complex geometries and ORCA HF-stretch values

configs/pipelines/phase4_orca_offline.yaml
  renders ORCA input and uses dummy fallback for review/CI

configs/pipelines/phase4_real_orca.yaml
  production-style ORCA opt/freq/SP configuration

configs/stages/dft_minima_orca_offline.yaml
configs/stages/sp_orca_offline.yaml
configs/stages/dft_minima_orca.yaml
configs/stages/sp_orca.yaml
  independent-stage configs for ORCA review or execution

tests/test_phase4_orca.py
  parser, input rendering, mocked ORCA subprocess, and offline pipeline tests
```

## How to run

Offline/review mode:

```bash
PYTHONPATH=. python -m hfauto.cli.main pipeline \
  --config configs/pipelines/phase4_orca_offline.yaml \
  --run-id phase4_offline
```

Real ORCA mode:

```bash
PYTHONPATH=. python -m hfauto.cli.main pipeline \
  --config configs/pipelines/phase4_real_orca.yaml \
  --run-id phase4_real
```

For a single stage, use:

```bash
PYTHONPATH=. python -m hfauto.cli.main run-stage dft-minima \
  --in runs/<run_id>/05_preopt/manifest.json \
  --out runs/<run_id>/06_dft-minima \
  --config configs/stages/dft_minima_orca_offline.yaml
```

## Validation

The Phase 4 test suite validates:

- ORCA parser fields from fixture output;
- ORCA Cartesian-coordinate extraction;
- ORCA input keyword rendering;
- mocked ORCA execution without a real binary;
- full offline Phase 4 pipeline with ORCA input rendering and dummy fallback;
- Phase 3 xTB/CREST fallback preprocessing;
- ranking/report contracts.

Validation command:

```bash
PYTHONPATH=. pytest -q
```

## Known limitations

- ORCA normal-mode projection for a true HF-stretch assignment is not implemented; `hf_stretch_cm1` is currently a high-frequency proxy.
- Thermochemistry uses ORCA parsed values or DFT thermal correction plus high-level SP. GoodVibes/quasi-RRHO is planned for Phase 6.
- TS search is still dummy/midpoint-based; real NEB-TS/OptTS/IRC is Phase 5.
- Psi4/PySCF/NWChem are not yet implemented as real DFT backends.

## Recommended next phase

Phase 5 should implement real TS search and IRC validation:

```text
ORCA NEB-TS input renderer and runner
FAST-NEB-TS / ZOOM-NEB-TS options
constrained proton-transfer scan -> OptTS fallback
IRC forward/backward execution
endpoint graph/RMSD matching
imaginary-mode overlap score
validated ΔG‡PT from real TS calculations
```

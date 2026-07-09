# Phase 3 Implementation: xTB / CREST backend integration

This phase converts the Phase 2 data-contract workflow into a practical pre-DFT structure generation workflow.

## Implemented

### Chemistry / simulation outputs

- CREST conformer backend with opt-in subprocess execution.
- RDKit fallback for CREST-disabled or CREST-missing environments.
- xTB pre-optimization backend with opt-in subprocess execution.
- Dummy fallback for xTB-disabled or xTB-missing environments.
- `species_preopt` artifacts that preserve canonical species IDs while pointing to optimized XYZ geometries.
- Geometry sanity QC for HF-containing species:
  - `r_HF_A`
  - `delta_r_HF_A`
  - `B_H_distance_A`
  - `B_H_F_angle_deg`
  - `hf_dissociated`
  - `proton_transferred_unintentionally`
- DFT and SP stages now prefer `species_preopt` geometries while keeping original reaction/species IDs for thermo aggregation.
- Descriptor extraction now prefers preoptimized reactant-complex geometries.

### Architecture

- `ConformersStage` now calls a conformer backend registry.
- `get_conformer_engine()` added to `hfauto.backends.registry`.
- CREST parsing utilities:
  - `parse_multi_xyz()`
  - `parse_crest_energies()`
- xTB parsing utility:
  - `parse_xtb_total_energy()`
- Subprocess execution is explicit and safe:
  - `allow_subprocess: true`
  - or `HFAUTO_ALLOW_XTB=1`
  - or `HFAUTO_ALLOW_CREST=1`
  - or `HFAUTO_ALLOW_SUBPROCESS=1`

## New pipeline configs

- `configs/pipelines/phase3_xtb_crest_offline.yaml`
  - CREST backend with RDKit fallback.
  - xTB backend with dummy fallback.
  - Runs without external binaries.

- `configs/pipelines/phase3_real_xtb_crest.yaml`
  - Real CREST/xTB subprocess mode.
  - Requires `crest` and `xtb` on PATH or explicit executable paths.

## Important Artifact contracts

### Conformer Artifact

```json
{
  "artifact_type": "conformer",
  "data": {
    "conformer_id": "mol00001_conf0000",
    "mol_id": "mol00001",
    "source": "crest | rdkit | rdkit_fallback_for_crest",
    "relative_energy_kcal_mol": 0.0,
    "boltzmann_weight_298K": 1.0,
    "xyz_path": ".../conf0000.xyz"
  }
}
```

### Preoptimized Species Artifact

```json
{
  "artifact_type": "species_preopt",
  "data": {
    "species_id": "spc_mol00001_...",
    "source_species_id": "spc_mol00001_...",
    "source_species_artifact_id": "spc_mol00001_...",
    "preopt_calc_id": "calc_...",
    "xyz_path": ".../xtbopt.xyz"
  },
  "qc": {
    "fallback_dummy": false,
    "geometry_qc": {
      "hf_dissociated": false,
      "proton_transferred_unintentionally": false
    }
  }
}
```

## Scientific review notes

Phase 3 does **not** yet provide final DFT-quality energies. Its purpose is to improve initial geometry quality before DFT and TS search. The scientifically useful outputs at this phase are:

- robust candidate/HF endpoint generation,
- conformer coverage metadata,
- preoptimized HF-complex geometry checks,
- early flags for HF dissociation or unintended proton transfer,
- production-compatible manifests for DFT and TS stages.

## Next phase recommendation

Proceed to **Phase 4: ORCA DFT opt/freq real backend**.

Recommended implementation order:

1. ORCA input renderer for `Opt Freq` and `SP`.
2. ORCA subprocess runner with opt-in execution.
3. ORCA parser for final energy, Gibbs correction, frequencies, optimized XYZ, dipole.
4. Scientific QC:
   - SCF convergence,
   - geometry convergence,
   - minima `n_imag == 0`,
   - TS `n_imag == 1`,
   - HF stretch extraction.
5. Add small scientific regression tests with mocked ORCA outputs.

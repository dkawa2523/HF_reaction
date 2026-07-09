# Stage input/output reference

This is a compact map, not a heavy runtime contract.  Runtime remains manifest based: each stage reads a `manifest.json` and writes a new `manifest.json`.

| Stage | Main inputs | Main outputs | Scientific purpose |
|---|---|---|---|
| ingest | SDF | molecule | Read and normalize candidates |
| enrich | molecule | molecule_enriched | Public-data support and identity metadata |
| detect-sites | molecule / molecule_enriched | site | Identify N/O/S/π basic sites |
| conformers | molecule, site | conformer | Candidate conformer ensemble |
| build-hf | conformer, site | species, reaction | Build B···(HF)n and ion-pair endpoints |
| preopt | species | species_preopt, calculation | xTB/fallback preoptimization and geometry QC |
| dft-minima | species_preopt/species | species_optimized, calculation | Real DFT minima/frequency evidence |
| ts-search | reaction, species | reaction_validated, TS species, calculation | TS candidate and frequency QC |
| irc | reaction_validated | reaction_path_validated | Endpoint connectivity check |
| sp | optimized/preopt species | calculation | Higher-level electronic correction |
| thermo | calculation, reaction path | thermo | G(T,p), ΔGassoc, ΔG‡, K |
| kinetics | thermo | kinetics | TST/Wigner rates and mechanism skeleton |
| calibrate | public data, thermo | method_validation/table | Public-reference support and method diagnostics |
| connector-audit | connector artifacts | audit tables | Production connector readiness |
| rank | thermo, kinetics, descriptors | rankings/tables | Screening, scientific, production ranking |
| viz | latest manifest | visualization bundle | Read-only review reports |
| ops | latest manifest | operations bundle | Read-only HPC/retry/reuse planning |

For developers, the lightweight list is also available in `hfauto/core/artifact_types.py` as `STAGE_CONTRACTS`.

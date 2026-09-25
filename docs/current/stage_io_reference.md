# Stage input/output reference

| Stage | Inputs | Outputs | Responsibility |
|---|---|---|---|
| ingest | SDF | molecule | Normalize supplied molecules |
| enumerate-states | molecule | molecule_state | Bounded protomer/tautomer hypotheses |
| conformers | molecule_state | conformer | CREST/RDKit intramolecular ensembles |
| build-complexes | conformer/species selectors | species | Config-bounded multicomponent seeds and CREST NCI complexes; exact-fingerprint definition checkpoints |
| preopt | species | species_preopt, preopt_geometry, calculation | Unbiased xTB relaxation; never overwrite the requested state after a production connectivity change |
| relaxation-discovery | preopt_geometry or CREST topology-stop trajectory | reaction_trial, reaction_discovery_attempt, product species_preopt, reaction_candidate | Normalize an unbiased connectivity-changing relaxation as a low-level, effectively-barrierless hypothesis; never claim a barrier |
| generate-reactions | species/preopt | reaction_trial | Product-free proposals; global budget is round-robin across selected source conformers |
| explore-reactions | reaction_trial | reaction_discovery_attempt, reaction_candidate, product species_preopt | ReaDuct NT2/AFIR discovery |
| discovery-audit | trials, attempts, candidates, basins, TS/IRC; optional comparison manifests | discovery_coverage, optional discovery_saturation | Measure finite-budget attrition and observed cross-budget plateaus; never claim exhaustive chemistry |
| dft-minima | species/preopt, reaction_candidate | species_optimized, calculation | Candidate endpoints plus bounded state anchors, exact trial sources, real same-method opt/freq evidence, resumable checkpoints; multi-formula budgets are round-robin by elemental composition/charge/multiplicity and use energy ordering only within one surface |
| minimum-registry | all accepted minima | minimum_basin, minimum_registry | Global same-PES basin identity |
| connect-minima | registry, accepted minima, reaction_trial | reaction_candidate | Check distinct DFT basins against declared continuous coordinates without inventing product geometry |
| reaction-plan | candidates, registry | reaction, basin_pair_assessment, reaction_case | Promote only different validated basins |
| recover-path | interrupted raw path, reaction, validated endpoints | reaction_path, path_attempt, saddle_attempt, reaction_case | Recover geometry-only TS seeds; path energy stays non-publishable |
| ts-search | reaction, case, endpoints | reaction_path, path_attempt, saddle_attempt, reaction_validated | Chemistry-aware XYZ_PATH where supported, budgeted/checkpointed NEB/string, and first-order saddle validation |
| irc | reaction_validated | reaction_path_validated | Bidirectional basin connectivity |
| reaction-classify | all path evidence | reaction_classification | Exclusive evidence-based path class |
| sp | selected species, optionally validated TS calculations | calculation | One fixed-geometry electronic energy; validated-TS mode excludes unvalidated saddle seeds |
| method-panel | validated reaction and real SP calculations | table, method_uncertainty | Fixed-geometry electronic sensitivity/spread; never claims stationary points on the alternative PES |
| thermo | validated calculations/path | thermo | Reaction and activation thermochemistry |
| thermo-sensitivity | validated frequency calculations and TS | thermo_sensitivity | Re-run external GoodVibes over declared low-frequency cutoffs; no internal fallback is accepted as a complete row |
| basin-populations | minimum_basin and production species_thermo | basin_population, table | Conditional Boltzmann populations over observed same-PES basins; never claims exhaustive coverage |
| reaction-rank | production thermo, optional method_uncertainty and basin_population | ranking | Like-for-like point, conservative-bound, or explicitly conditional population-adjusted ranking |

The lightweight machine-readable summary is `hfauto/core/artifact_types.py`.

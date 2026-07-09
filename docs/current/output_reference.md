# Output reference

## Ranking outputs

```text
rank_screening.csv    all rows, including exploratory and fallback rows
rank_scientific.csv   rows passing the real-DFT scientific gate
rank_production.csv   rows passing validated TS/IRC and production thermo gate
rank_scavenger.csv    HF capture-oriented score with eligibility flags
rank_activation.csv   HF activation-oriented score with eligibility flags
cluster_risk.csv      HF cluster growth / particle-risk proxy
candidate_summary.csv best row per candidate with next actions
ranking_summary.json  run-level ranking readiness summary
```

## Thermochemistry outputs

```text
species_thermo.csv     species-level G(T,p), pressure correction, quasi-RRHO metadata
reaction_thermo.csv    ΔGassoc, ΔGionpair, ΔG‡, K, quality tier, confidence
thermo_records.jsonl   machine-readable thermo artifacts
```

Important distinction:

```text
ΔG_assoc_standard_kcal_mol
  standard-state association free energy

ΔG_assoc_pressure_corrected_kcal_mol
  process partial-pressure corrected association free energy
```

For gas processes, the pressure-corrected value is often the more relevant decision variable.

## Kinetics outputs

```text
kinetics_records.csv       TST/Wigner rate constants and quality metadata
reactor_screening.csv      simple residence-time conversion proxy
cantera_mechanism.yaml     draft mechanism or production mechanism, depending on connector status
cantera_validation.json    Cantera connector readiness and validation status
```

## Visualization outputs

```text
15_viz/report.html                         campaign report
15_viz/molecules/<mol_id>/dossier.html     molecule dossier
15_viz/reactions/<reaction_id>/dossier.html reaction dossier
15_viz/networks/reaction_network*.html     reaction network views
```

## Operations outputs

```text
16_ops/array_job_plan.csv                 artifact-level job planning
16_ops/retry_plan.csv                     failed/fallback retry candidates
16_ops/reuse_plan.csv                     conservative reuse candidates
16_ops/qcarchive_payloads.jsonl           QCArchive review payloads
16_ops/operations_report.html             operations dashboard
```

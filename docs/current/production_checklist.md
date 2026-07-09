# Production checklist

Use this checklist before candidate decisions or external reporting.

## Scientific readiness

- [ ] `rank_production.csv` contains the rows being discussed.
- [ ] `production_rank_eligible = true` for selected rows.
- [ ] `main_values_are_dummy = false`.
- [ ] `main_values_are_fallback = false`.
- [ ] DFT minima have converged SCF and geometry optimization.
- [ ] Minima have `n_imag = 0`.
- [ ] TS has `n_imag = 1` and an appropriate imaginary mode.
- [ ] IRC connects expected reactant and product endpoints.
- [ ] Low-frequency thermal correction method is recorded.
- [ ] Temperature and partial pressure assumptions are recorded.
- [ ] Top candidates have high-level single-point correction or a documented reason not to.

## Chemical review

- [ ] HF cluster size `n` is appropriate for the process question.
- [ ] Competing sites in the molecule were considered.
- [ ] Conformer coverage was sufficient for flexible candidates.
- [ ] Strong HF capture and cluster growth risk were reviewed separately.
- [ ] Activation and scavenger rankings were not conflated.

## Process review

- [ ] Vapor pressure / boiling point / delivery feasibility was reviewed.
- [ ] EHS flags were reviewed by an appropriate owner.
- [ ] Candidate stability under process temperature was considered.

## Operations review

- [ ] Connector audit is reviewed.
- [ ] Retry plan is empty or justified.
- [ ] Reuse candidates are reviewed before applying.
- [ ] HPC logs are archived for production runs.

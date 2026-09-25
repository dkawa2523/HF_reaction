# Prioritized extension backlog

This file lists scientific work that remains. The generic artifact contracts,
bounded discovery audit, fixed-geometry method uncertainty, path recovery, and
resumable TS execution are implemented and tested; they are not backlog items.

## P0: finish the TMA(HF)2 stationary-point proof

- Complete the active same-PES saddle refinement from the preserved NEB internal
  maximum. The stopped NEB profile is a geometry seed only and its energy rise
  is not an activation barrier.
- If the bounded string reaches walltime first, retain only a method-identity-
  checked, atom-order-valid single-maximum geometry seed. The implemented
  interrupted-path contract explicitly leaves normal termination and path
  energy publication false.
- Accept a TS only with one significant imaginary frequency aligned with the
  path or declared coordinate and above both registered minima.
- Run bidirectional pysisyphus/NWChem IRC and independently reoptimize both
  termini into the registered basins.
- If any gate fails, preserve the failure and route the case to another bounded
  path family or an intermediate search. Do not publish an activation barrier.

## P1: expand discovery coverage

- The bounded pilot contract covers NH3, TMA, aniline, and pyridine with
  HF1/HF2/HF3 plus a TMA/HF1-3/explicit-water microsolvation series in
  `configs/pipelines/amine_hf_pilot_hf1_hf3.yaml` (15 compositions, at most 30
  trials and 60 ReaDuct attempts). Input ingestion, neutral-state enumeration,
  gas-phase conformers, 15 NCI compositions, and exact atom counts are
  validated. The v2 ID-collision evidence remains preserved, while the clean
  fingerprint-checkpointed v3 campaign completed 17/17 definitions and its
  downstream audit completed with 30 automatic trials, 58 ReaDuct attempts,
  and 5 ReaDuct products. One additional TMA(HF)3 candidate came from an
  unbiased connectivity-changing relaxation.
- Complete the active checkpointed DFT endpoint follow-up in
  `configs/pipelines/amine_hf_pilot_dft_followup.yaml`. Its four-candidate
  budget is distributed across molecular composition/electronic-state surfaces
  rather than by scientifically invalid cross-formula absolute energies. Only
  distinct frequency-validated basin pairs may proceed to TS/IRC.
- Repeat selected cases at larger NCI-seed, trial, NT2/AFIR, and DFT-anchor
  budgets. The cross-manifest comparison contract now records count increments
  and observed plateaus, but the repeated real campaigns remain to be run.
  Search priority and a finite plateau are never chemical evidence of
  exhaustiveness.
- Add at least two chemically different positive TS/IRC regressions beyond
  HCN/HNC and one expected same-basin or effectively barrierless case.
  The first generic conformational regression is configured as trans/cis
  HONO with an explicit dihedral trial in
  `configs/pipelines/hono_isomerization_nwchem.yaml`. Both distinct minima are
  frequency validated and its connectivity-preserving v3 NEB is active; TS/IRC and a second
  non-HCN completion remain required. A chemically different intramolecular
  hydrogen-migration contract (formaldehyde -> trans-hydroxymethylene) is now
  validated through trial generation in
  `configs/pipelines/formaldehyde_hydroxymethylene_nwchem.yaml`; its real
  minima/TS/IRC campaign has not started.
  An independent two-start water regression is configured in
  `configs/pipelines/same_basin_water_nwchem.yaml`; it passed with two
  frequency-validated optimizations collapsed into one registry basin and
  therefore never launched TS/IRC.

## P2: extend electronic-structure validation

- The HCN/HNC fixed-geometry method panel is complete and records a 0.9824
  kcal/mol electronic-barrier range. This measures electronic sensitivity but
  does not establish stationary points on another PES.
- Reoptimize and frequency-check a representative alternative-method endpoint
  and TS subset. The fixed-geometry integration-grid sensitivity is complete:
  fine versus xfine shifts the HCN/HNC electronic barrier by only 0.0000586
  kcal/mol. The
  B3LYP-D3/def2-TZVPD HCN/HNC endpoints are now reoptimized and frequency
  validated in `configs/pipelines/hcn_alternative_pes_validation.yaml`; its
  same-PES v3 NEB/TS/IRC remains.
  The external-GoodVibes 50/100/150 cm-1 panel is complete with a 0.00311
  kcal/mol Gibbs activation range for HCN/HNC.
- Extend the implemented conditional basin-population model with discovery
  saturation evidence before treating any ensemble as complete. Ranking can
  report or explicitly apply the finite-ensemble population penalty; it never
  silently assumes exhaustive conformer coverage.

## P3: scale production execution

- ReaDuct now has per-attempt input/raw-hash checkpoints, total-attempt and
  stage-boundary walltime budgets, and integrity-checked resume tests. Connect
  these contracts to a real Slurm/PBS campaign and validate restart after an
  actual scheduler preemption; neither scheduler exists in the current WSL
  environment.
- Add solvent or surface models only as separate, explicit model layers. The
  molecular gas-phase core must not silently imply wet-etch equipment behavior.

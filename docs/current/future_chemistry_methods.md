# Future chemistry and simulation methods

This document lists useful future chemistry extensions.  They are deliberately kept outside the current minimal baseline so the code remains maintainable.

## 1. More robust HF complex generation

Current templates are sufficient as seeds.  Future improvements:

- multi-direction HF placement around each basic site
- steric-grid placement before xTB/CREST filtering
- dual-site bridged HF motifs
- explicit `(HF)n` chain and cyclic motifs
- automatic pruning by H-F distance, B-H distance, and clash score

Expected benefit: fewer missed low-energy `B···(HF)n` complexes and better DFT/TS success rates.

## 2. HF cluster-size convergence

Add systematic `n = 1, 2, 3, 4` convergence checks:

```text
ΔG_assoc(n)
ΔG_ionpair(n)
ΔG_addHF(n -> n+1)
ΔG‡PT(n)
```

Expected benefit: distinguish true HF activation from irreversible HF clustering or particle-risk behavior.

## 3. Si-O probe reactions for semiconductor relevance

For HF activation additives, add gas-phase probe reactions:

```text
B···HF + Si(OH)4 -> fluorinated silanol model products
B···HF + (HO)3Si-O-Si(OH)3 -> siloxane cleavage / fluorination models
```

Expected benefit: connect HF activation descriptors to a process-relevant Si-O reaction proxy without immediately moving to periodic surface DFT.

## 4. Method ladder and uncertainty

Add a formal method ladder:

```text
xTB/CREST screening
r2SCAN-3c opt/freq
ωB97X-D4 or ωB97M-V single point
Double-hybrid or DLPNO-CCSD(T) final SP for top candidates
```

Track method sensitivity per candidate.  Do not hide method disagreements behind one score.

## 5. Automated TS fallback hierarchy

Keep a small, explicit fallback sequence:

```text
NEB-TS
constrained scan -> OptTS
pysisyphus GSM/NEB
manual review queue
```

Expected benefit: improved TS success without turning the workflow into a black-box reaction discovery system.

## 6. True normal-mode and path parsing

Current visualization supports approximate proton-transfer arrows.  Future work:

- parse ORCA Hessian/frequency normal modes
- export true imaginary-mode vectors
- parse full IRC/NEB multi-frame paths
- compute mode overlap from real normal-mode displacement

Expected benefit: stronger TS validation and better reaction dossier figures.

## 7. Microkinetic and reactor integration

After production thermochemistry and validated TS data are available:

- Arkane TST/RRKM/Master Equation for selected reactions
- Cantera mechanisms with reviewed thermo data
- batch/PFR-like process screening over temperature, pressure, and residence time

Expected benefit: convert ΔG and ΔG‡ into process-relevant rates and sensitivity trends.

## 8. ML and surrogate models

Use ML only as acceleration, not as final evidence:

- ML potentials for conformer/path seed generation
- active learning for failed or high-uncertainty candidates
- surrogate ranking after DFT-calibrated reference campaigns

Expected benefit: scale to larger SDF libraries while retaining DFT validation for top candidates.

## 9. Periodic or cluster-surface models

For final process relevance, consider optional higher-level models:

- silanol cluster models
- amorphous silica cluster fragments
- periodic slab DFT for final candidates

Expected benefit: bridge gas-phase HF activation to surface process chemistry.

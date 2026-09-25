# Production science validation plan

This plan separates workflow correctness from claims of chemical coverage or
predictive accuracy.

## Molecular panel

Use ammonia, methylamine, dimethylamine, trimethylamine, triethylamine,
aniline, pyridine, imidazole, DABCO, and dimethyl ether. For the bounded pilot,
start with NH3, TMA, aniline, and pyridine. Evaluate HF1, HF2, HF3, and an
explicit water-assisted branch as separate compositions; never merge their
standard states or electronic surfaces implicitly.

## Evidence to collect

- CREST NCI seed and conformer counts by composition.
- Generated trial count and executed NT2/AFIR count.
- Low-level product yield, accepted DFT minima, and global minimum basins.
- First-order TS and bidirectional IRC completion.
- Electronic, ZPE, enthalpy, quasi-harmonic Gibbs, and method-sensitivity ranges.
- Proton-affinity and gas-basicity anchors as context only, not surrogate
  reaction barriers.

Run `discovery-audit` after each bounded branch. Its zero-candidate conclusion
is `no_distinct_product_observed_within_budget`, never "no reaction exists."

The completed v3 pilot contains 15 NCI compositions, 30 automatic trials, 58
ReaDuct attempts, and five ReaDuct products. One further TMA(HF)3 product was
observed during unbiased xTB relaxation. The active DFT follow-up samples four
different molecular surfaces; these counts describe finite observed coverage,
not an exhaustive reaction catalogue.

## Method ladder

1. GFN2-xTB/CREST and ReaDuct for bounded discovery.
2. PBE0-D3/def2-SVPD opt/freq, path, TS, and IRC as the current common PES.
3. Diffuse triple-zeta fixed-geometry single-point sensitivity.
4. Alternative-method opt/freq and TS revalidation on a representative subset.
5. Higher-level single points for selected small references when resources
   permit.

## Decision criteria

Candidate ranking is production-eligible only when both endpoint minima, one
first-order saddle, bidirectional IRC connectivity, and external
thermochemistry pass. Method spread and conformer/basin populations must
accompany cross-candidate ranking once more than one candidate is compared.

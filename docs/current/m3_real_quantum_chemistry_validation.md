# M3 real quantum-chemistry validation

## Result

HCN -> HNC is the positive end-to-end benchmark. Real calculations cover unknown
product discovery, two DFT minima, a DFT path, a first-order saddle, bidirectional
IRC, endpoint reoptimization, and activation thermochemistry. No dummy or fallback
result is accepted as evidence.

## Common PES and scientific gates

- ReaDuct 6.1.0 + GFN2-xTB for discovery;
- NWChem 7.2.3 PBE0-D3(0)/def2-SVPD for minima, path, TS, frequencies, and IRC
  gradients;
- pysisyphus 1.0.0 + QCEngine 0.50.0 for bidirectional IRC;
- charge 0, multiplicity 1 throughout.

| gate | evidence | result |
|---|---:|---|
| HCN/HNC minima | both `n_imag = 0` | pass |
| basin identity | graph + permutation RMSD distinct | pass |
| TS energy | -93.1669942844 Eh | recorded |
| TS imaginary mode | -1131.57 cm-1, exactly one | pass |
| mode/path alignment | mass-weighted overlap 0.547 | pass (limit 0.50) |
| bidirectional IRC | real NWChem gradients, both branches | pass |
| optimized IRC endpoints | registered HCN/HNC basins | pass |
| maximum endpoint RMSD | 2.27e-5 A | pass |

## Independent thermochemistry cross-check

The GoodVibes 4.3 NWChem adapter creates a hash-bound parser-compatibility copy
inside the run. The raw NWChem output is unchanged. The copy only normalizes
NWChem metadata spellings that GoodVibes cannot parse for a linear molecule;
electronic energies, frequencies, rotational values, masses, and geometries are
not changed. Source/prepared SHA-256 hashes and applied rules are recorded.

The external GoodVibes result in `runs/hcn_goodvibes_v4` is production-ready:

| forward HCN -> HNC barrier at 298.15 K | kcal/mol |
|---|---:|
| electronic | 46.6672 |
| E + scaled ZPE | 43.3008 |
| quasi-harmonic enthalpy | 43.4756 |
| quasi-harmonic Gibbs | 42.2809 |

The external Gibbs barrier differs from the earlier internal NWChem-derived
value (42.2296 kcal/mol) by about 0.0513 kcal/mol. This closes the independent
thermochemistry implementation check; it does not replace method/basis
sensitivity analysis.

## Fixed-geometry electronic-structure sensitivity

`runs/hcn_method_sensitivity_v1` completed three real NWChem single-point
panels on the same validated HCN, HNC, and TS geometries:

| method | electronic barrier (kcal/mol) | reaction energy (kcal/mol) |
|---|---:|---:|
| PBE0-D3/def2-SVPD | 46.6672 | 13.2353 |
| PBE0-D3/def2-TZVPD | 46.4714 | 13.8780 |
| B3LYP-D3/def2-TZVPD | 47.4539 | 14.1325 |

The electronic barrier range is 0.9824 kcal/mol (half-range 0.4912; maximum
absolute shift from the reference 0.7866 kcal/mol). The reaction-energy range
is 0.8972 kcal/mol. These are fixed-geometry sensitivities only: the
alternative methods have not reoptimized or frequency-validated all three
stationary points.

With `uncertainty_policy: conservative_bound`, the production ranking keeps
the GoodVibes Gibbs barrier as the 42.2809 kcal/mol point estimate and adds the
0.7866 kcal/mol maximum electronic shift, yielding a transparent conservative
ranking value of 43.0675 kcal/mol. The point estimate and uncertainty remain
separate columns in the artifact.

The finite-ensemble population stage applied to the two observed basins gives
conditional 298.15 K populations of 0.99999999954 for HCN and 4.57e-10 for HNC
from the 12.7418 kcal/mol standard Gibbs difference. This is conditional on the
observed basin set and is not an exhaustive conformer claim.

## Evidence locations

- discovery and minima: `runs/readuct_hcn_validation`,
  `runs/hcn_readuct_nwchem_validation`;
- DFT path and TS: `runs/hcn_ts_irc_validation`, `runs/hcn_saddle_v3`;
- IRC: `runs/hcn_irc_v5`;
- external thermochemistry: `runs/hcn_goodvibes_v4`;
- method sensitivity: `runs/hcn_method_sensitivity_v1`;
- conditional populations: `runs/hcn_goodvibes_v4/12_basin-populations`.

## Remaining validation

- Reoptimize and frequency-check a representative alternative-method subset;
  evaluate grid and low-frequency sensitivity.
- Add at least two chemically different positive TS/IRC regressions and one
  same-basin or effectively barrierless regression.
- Treat HCN/HNC as workflow validation, not as an accuracy guarantee for amine/HF
  chemistry.

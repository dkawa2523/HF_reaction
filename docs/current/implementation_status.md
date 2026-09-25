# Implementation and evidence status

This status deliberately separates software completion from a completed
chemical claim.

| Area | Software status | Real evidence status |
|---|---|---|
| config-bounded multicomponent complexes | implemented | TMA + 2 HF: 4 CREST NCI complexes; amine/HF pilot: 15/15 NCI compositions including HF1-3 and explicit water, with exact atom counts and no ID collisions |
| product-free reaction trials | implemented | TMA(HF)2: 12 trials; HCN/HNC: 1 trial; amine/HF pilot: 30 generated trials plus 1 unbiased-relaxation trial |
| ReaDuct NT2/AFIR discovery | implemented, no dummy fallback | TMA(HF)2: 24 attempts and 0 low-level products; HCN/HNC: 1 accepted product; amine/HF pilot: 58 ReaDuct attempts and 5 accepted low-level products |
| DFT minimum registry | implemented globally | TMA(HF)2: 2 distinct minima; HCN/HNC: 2 distinct minima |
| NEB/string and TS frequency gate | implemented | HCN/HNC passed; TMA(HF)2 recovered seed had 4 imaginary modes and was rejected; adaptive string is active |
| bidirectional IRC and endpoint reoptimization | implemented | HCN/HNC passed; TMA(HF)2 awaits a valid TS |
| activation thermochemistry | implemented with external GoodVibes | HCN/HNC passed; TMA(HF)2 cannot be reported yet |
| finite-budget discovery audit | implemented | TMA and HCN branches evaluated |
| fixed-geometry method uncertainty | implemented | HCN/HNC: 3 methods complete; electronic barrier range 0.9824 kcal/mol |
| finite-ensemble basin populations | implemented | HCN/HNC: conditional 298 K populations computed over 2 observed basins |
| TS attempt budgets/checkpoints/recovery | implemented and unit-tested | one interrupted TMA NEB profile recovered as a geometry-only seed |
| interrupted NEB/string seed recovery | method identity, normal termination, and path-energy publication are separate and regression-tested | future TMA continuation can preserve a timed-out same-PES path as geometry only; it still requires independent saddle/frequency/IRC validation |
| ReaDuct attempt budgets/checkpoints/recovery | implemented with input/raw hash verification | production pilot uses per-attempt checkpointing and a 4 h stage-boundary budget |
| explicit conformational reaction coordinates | implemented for distance, angle, and dihedral trials | trans/cis HONO: 2 distinct PBE0-D3/def2-SVPD minima; NEB active |
| chemistry-aware initial paths | implemented for acyclic dihedral rotation and bonded-center single-atom transfer | HONO v3 and alternative-PES HCN v3 use validated NWChem XYZ_PATH inputs |
| expected same-basin regression | implemented | two water starts: 2 real minima calculations assigned to 1 basin; TS/IRC correctly suppressed |
| alternative-PES stationary-point validation | configured and running | B3LYP-D3/def2-TZVPD HCN/HNC minima both passed; NEB active |
| integration-grid sensitivity | implemented with explicit numerical provenance | HCN/HNC fixed-geometry fine/xfine panel complete; activation shift 0.0000586 kcal/mol and reaction-energy shift 0.000000203 kcal/mol |
| low-frequency thermochemistry sensitivity | implemented with external GoodVibes only | HCN/HNC 50/100/150 cm-1 panel complete; activation range 0.00311 kcal/mol |
| cross-budget discovery comparison | implemented | observed yield plateaus and count increments are reportable; exhaustive claims remain forbidden |
| second non-HCN positive regression | contract and chemistry-aware H-migration path implemented | formaldehyde/trans-hydroxymethylene trial generation passed; real minima/TS/IRC not started |
| definition-level complex checkpoint/resume | implemented with request fingerprints and non-overwriting raw attempts | amine/HF v3 completed 17/17 definitions with 0 failures |
| preopt job checkpoint/resume | implemented with input/final hashes and non-overwriting attempts | amine/HF v3 completed 15/15 NCI jobs; 14 retained topology and 1 was routed to relaxation discovery |
| ensemble-aware preopt selection | implemented with generic required source-data keys | amine/HF and TMA(HF)2 optimize only CREST NCI members, not their placement seeds |
| connectivity-changing relaxation discovery | implemented for xTB preopt and CREST topology safety stops | TMA(HF)3 produced one composition-preserving H-F break/N-H form candidate; it is an effectively-barrierless low-level hypothesis, not a barrier |
| composition-balanced DFT promotion | implemented and regression-tested | four amine/HF molecular surfaces selected round-robin; PBE0-D3/def2-SVPD endpoint follow-up is active |

## Current TMA(HF)2 interpretation

The two accepted PBE0-D3/def2-SVPD minima are separated by 0.9367 kcal/mol in
electronic energy. This is a minimum-to-minimum energy difference, not an
activation barrier. The preserved stopped NEB profile has one internal maximum
at image 5 and an unrelaxed rise of 10.5014 kcal/mol, but
`path_energy_publishable=false`. It may only seed an independently optimized
first-order saddle.

The bounded low-level search conclusion is
`no_distinct_product_observed_within_budget`. The DFT anchor ensemble separately
observed one distinct-basin candidate. These statements are not contradictory
and neither proves exhaustive absence or presence of other products.

The recovered image-5 seed has four significant imaginary frequencies
(-1112.83, -1016.93, -102.85, and -64.64 cm-1). The -102.85 cm-1 mode has a
0.7728 overlap with the local path tangent, but a four-negative-mode point is
not a first-order saddle. It was rejected and routed to an adaptive same-PES
string calculation; no activation barrier was created from this seed.

The adaptive string has completed multiple geometry-update cycles. Its path remains
geometry evidence only until NWChem reports convergence and the independently
refined maximum passes the one-imaginary-mode and IRC gates.

## Additional real regressions

- trans/cis HONO optimized to two distinct PBE0-D3/def2-SVPD minima with
  `n_imag=0`; their electronic energy separation is about 0.1552 kcal/mol.
  The declared H-O-N-O dihedral changes by pi radians and the pair was promoted
  to a reaction. Naive Cartesian and wrong-fragment torsional paths were
  preserved and rejected. The v3 NEB uses a connectivity-preserving torsional
  XYZ path; TS/IRC remain in progress and its path energy is not yet a barrier.
- Two independently distorted water inputs both optimized with `n_imag=0` and
  electronic energies differing by about 1.7e-7 kcal/mol. The minimum registry
  assigned both to one basin (`minimum_count=2`, `basin_count=1`).
- B3LYP-D3/def2-TZVPD reoptimization gave frequency-validated HCN and HNC
  minima at -93.464698107813 and -93.442111713185 hartree. The corrected C,N,H
  atom mapping passed the connection gate. Cartesian and distance-compressing
  paths were preserved and rejected; v3 uses a bonded-center polar transfer
  arc and remains under same-PES NEB validation.

The external GoodVibes cutoff panel reused the three accepted HCN/HNC/TS raw
frequency outputs. At 50, 100, and 150 cm-1 it gave forward Gibbs barriers of
42.2803, 42.2809, and 42.2834 kcal/mol, respectively. This 0.00311 kcal/mol
range is evidence for this benchmark only, not a universal cutoff error.

The self-contained fine/xfine integration-grid panel records `grid` and the SCF
energy threshold in every NWChem calculation artifact. At fixed validated
HCN/HNC/TS geometries, changing xfine to fine shifted the electronic activation
energy by -0.0000586 kcal/mol and the reaction energy by +0.000000203 kcal/mol.
These are numerical-integration sensitivities only; they do not revalidate
stationary points on another PES.

## Pilot discovery correction

The first pilot exposed repeated-fragment seed collisions. The placement engine
now uses the covalent neighbour graph and deterministic rigid-body collision
alternatives. A subsequent v2 run exposed a separate artifact-identity defect:
an NCI conformer reused the starting molecule's `mol_id`, overwrote that
molecule's conformer revision, and contaminated later HF2/HF3 compositions.
For example, NH3/HF3 contained 16 atoms instead of 10 and TMA/HF2 contained 19
instead of 17. Therefore v2 HF2/HF3 structures are explicitly invalid and must
not enter preopt or reaction discovery; their raw files remain negative
provenance. Species CREST IDs are now species-scoped, multicomponent states do
not inherit one fragment's `mol_id`, and tests enforce the 6/8/10-atom NH3/HF
series. A clean, definition-checkpointed v3 campaign completed all 17
definitions and all 15 requested NCI compositions with zero failures. The
exact expected atom counts were recovered for NH3/HF1-3 (6/8/10), TMA/HF1-3
(15/17/19), aniline/HF1-3 (16/18/20), pyridine/HF1-3 (13/15/17), and
TMA/HF1-3/H2O (18/20/22). All 15 NCI artifact IDs are unique.

The downstream finite campaign preoptimized exactly these 15 NCI members,
generated 30 automatic trials, and recorded 59 discovery attempts including
the unbiased relaxation. Six low-level candidates were accepted: five from
NT2/AFIR and one from the TMA(HF)3 unbiased proton-transfer relaxation. The
coverage audit explicitly reports zero DFT minima, zero validated TSs, and zero
IRC-connected reactions at this point; candidate count is not a reaction-rate
or activation-barrier claim.

The DFT follow-up distributes its four-candidate budget round-robin across
elemental composition, charge, and multiplicity surfaces. Discovery energies
are used for ordering only within one surface; absolute energies of different
formulas are never compared. Eight exact endpoints are now under checkpointed
PBE0-D3/def2-SVPD optimization. The production configuration includes the
explicit-water series without introducing a wet-etch equipment model.

## Current positive reference

HCN -> HNC has two frequency-validated minima, one first-order saddle, real
bidirectional IRC with reoptimized endpoint matches, and external GoodVibes
thermochemistry. At 298.15 K its PBE0-D3/def2-SVPD forward values are 46.6672
kcal/mol electronic, 43.3008 kcal/mol including scaled ZPE, 43.4756 kcal/mol
quasi-harmonic enthalpy, and 42.2809 kcal/mol quasi-harmonic Gibbs activation
energy.

## Evidence locations

- `runs/tma_hf2_no_product_v2/08_discovery-audit`
- `runs/tma_hf2_discovery_audit_v2/00_discovery-audit`
- `runs/tma_hf2_recovered_path_v1`
- `runs/hcn_goodvibes_v4`
- `runs/hcn_goodvibes_v4/12_basin-populations`
- `runs/hcn_method_sensitivity_v1`
- `runs/hono_isomerization_v1`
- `runs/same_basin_water_v1`
- `runs/hcn_alternative_pes_v1`
- `runs/hcn_alternative_pes_v3`
- `runs/hcn_thermo_sensitivity_v1`
- `runs/hcn_grid_sensitivity_v2`
- `runs/hono_isomerization_v3`
- `runs/formaldehyde_hydroxymethylene_contract_v1`
- `runs/amine_hf_pilot_hf1_hf3_v1`
- `runs/amine_hf_pilot_hf1_hf3_v2`
- `runs/amine_hf_pilot_hf1_hf3_v3`
- `runs/amine_hf_pilot_hf1_hf3_v3_discovery`
- `runs/amine_hf_pilot_dft_followup_v1`

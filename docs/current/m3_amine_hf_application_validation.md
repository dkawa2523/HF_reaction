# TMA(HF)2 application validation

## Scope

This case studies a finite, neutral trimethylamine + 2 HF complex on a gas-phase
molecular PES. It does not model liquid transport, a surface, reactor operation,
or an etch rate. TMA, the HF count, and proton transfer are configuration inputs;
they are not hard-coded into the discovery, basin, or path backends.

## Bounded low-level discovery

`configs/pipelines/m3_trimethylamine_hf2_nwchem.yaml` generated four CREST NCI
complexes and distributed twelve product-free trials across them. The final
evidence in `runs/tma_hf2_explore_v11` contains:

- 4 CREST NCI complexes within 1.519 kcal/mol;
- 6 unbiased xTB-preoptimized structures;
- 12 reaction trials and 24 ReaDuct attempts;
- 0 chemically distinct low-level product candidates;
- 11 equivalent-atom permutations, 8 unrealized requested changes, and 5 NT2
  attempts without a TS guess.

Every AFIR endpoint was optimized again after removing the artificial force.
AFIR-biased energies are excluded from barrier fields. This remains a valid
negative result for the bounded low-level search, but it is no longer treated as
proof that the DFT surface contains only one basin.

## DFT anchor ensemble and basin result

The `dft-minima` stage now selects a bounded `anchor_minima_per_state` ensemble
independently of low-level product discovery. This removes the previous blind
spot in which `0 xTB products` prevented DFT examination of chemically different
NCI/conformer seeds.

Two pre-existing real NWChem 7.2.3 PBE0-D3(0)/def2-SVPD opt/freq results were
re-registered in `runs/tma_hf2_minima_connection_v2`:

| minimum | electronic energy (Eh) | n_imag | result |
|---|---:|---:|---|
| neutral-like complex | -374.753442395264 | 0 | accepted DFT minimum |
| shared-proton/NCI complex | -374.754935121294 | 0 | accepted DFT minimum |

They have the same composition, charge 0, multiplicity 1, and method lineage,
but different registry basin identities. The more recent NCI_0 result at about
-374.754932 Eh is a duplicate of the shared-proton basin (permutation RMSD about
0.003 A), not a third basin. The neutral-like minimum lies 0.9367 kcal/mol above
the shared-proton minimum at this electronic-structure level; this endpoint
difference is not an activation barrier.

## Connecting minima without inventing a product

The new `connect-minima` stage compares every eligible DFT basin pair against
pre-existing, product-free `reaction_trial` coordinates. It checks both
directions and emits a `reaction_candidate` only when:

- both endpoints are real frequency-validated minima on exactly one PES;
- the basins are distinct;
- composition and atom ordering are preserved;
- the declared continuous coordinate moves in a consistent direction.

For this TMA(HF)2 pair, the N-H/H-F coordinated distance progress is 0.0844 A.
No individual distance changes by the 0.15 A required for a full bond-rearrangement
claim, so the candidate is conservatively labeled `coordinate_reorganization`
with no asserted bond changes. It is not called a completed proton transfer.
`reaction-plan` nevertheless promotes it as a scientifically valid distinct-basin
path hypothesis.

## Current path status and acceptance rule

`configs/pipelines/tma_hf2_validated_minima_path.yaml` drives the same-PES path,
TS, frequency, IRC, classification, thermochemistry, and ranking stages. Raw path
work is stored under `runs/tma_hf2_path_v1`; it must not be interpreted as a
barrier until the following sequence completes:

1. a converged string/NEB path contains a resolved internal maximum;
2. saddle refinement yields exactly one significant imaginary frequency;
3. the imaginary mode agrees with the path or declared coordinate;
4. bidirectional IRC termini independently reoptimize into the two registered
   TMA(HF)2 basins.

Failure at any gate produces a typed unresolved or multistep result, not a false
activation barrier. Two distinct minima alone do not prove that one elementary
first-order saddle connects them.

The recovered maximum from the stopped 11-image NEB was tested by a real
PBE0-D3/def2-SVPD analytic Hessian. It has four significant imaginary
frequencies (-1112.83, -1016.93, -102.85, and -64.64 cm-1). Although the
-102.85 cm-1 mode overlaps the local path tangent by 0.7728, the geometry is a
higher-order saddle seed, not a first-order TS. The workflow therefore rejected
it and routed the case to a 9-image adaptive NWChem string calculation. No
barrier or IRC claim is made from the rejected seed.

## Remaining application work

- Finish the TMA(HF)2 path/TS/IRC proof above.
- Repeat bounded discovery at several NCI-anchor and reaction-trial budgets to
  quantify search incompleteness.
- Extend composition panels to HF1/HF2/HF3, water-assisted clusters, and multiple
  amine classes.
- Evaluate method, basis, dispersion, cluster-size, and low-frequency sensitivity
  before using results to rank candidate gases.

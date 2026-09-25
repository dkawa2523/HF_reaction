# Generic molecular reaction discovery platform

## Purpose and scope

The platform starts from supplied molecules, enumerates molecular states and
conformers, builds non-covalent complexes, discovers previously unknown product
minima, and validates an elementary reaction with NEB/string, a first-order
saddle, and bidirectional IRC.  It supports amine/HF candidate studies but does
not encode TMA, HF2, or proton transfer as a core assumption.

The modeled object is a gas-phase molecule or finite molecular complex.  Solid
surfaces, liquid transport, reactor flow, and wafer-scale wet-etch physics are
out of scope.  Therefore results rank intrinsic molecular reaction hypotheses;
they do not directly predict an etch rate or selectivity in equipment.

## One-way artifact flow

```text
ingest
 -> enumerate-states
 -> conformers
 -> build-complexes
 -> preopt
 -> relaxation-discovery
 -> generate-reactions
 -> explore-reactions
 -> dft-minima
 -> minimum-registry
 -> connect-minima
 -> reaction-plan
 -> ts-search
 -> irc
 -> reaction-classify
 -> thermo
 -> basin-populations
 -> reaction-rank
```

`generate-reactions` creates product-free `reaction_trial` requests.  Automatic
features propose proton transfer, relays of at most two donor steps,
substitution, association, and dissociation coordinates.  Reviewed explicit
coordinates are normalized into exactly the same contract.

`explore-reactions` executes bounded ReaDuct NT2 and AFIR attempts.  NT2 tries a
small ordered set of TS optimizers.  If the low-level saddle fails but the final
driven frame relaxes without a force into a distinct xTB endpoint, that endpoint
may become a candidate, explicitly marked as having no low-level TS/IRC claim.
AFIR endpoints are always reoptimized after removing the artificial force.
Biased AFIR energy is excluded from every barrier field.

`preopt` never silently relabels a connectivity-changing production relaxation
as the input state. `relaxation-discovery` converts such a result (including a
CREST initial optimization that safety-stopped after a topology change) into a
standard low-level candidate. It may support an effectively-barrierless
hypothesis on that low-level surface, but it does not supply a DFT barrier, TS,
or IRC claim.

## Public evidence contracts

- `reaction_trial`: source state, charge/multiplicity, atom-pair associations
  and dissociations, driver order, rationale, priority, and hard attempt budget.
- `reaction_discovery_attempt`: backend versions, convergence, trajectory,
  hashes, raw paths, unbiased endpoint status, and explicit failure reason.
- `reaction_candidate`: composition-preserving low-level structural change only;
  this is not a reaction or barrier.
- `minimum_basin`: one or more real frequency-validated DFT minima with a common
  composition, electronic state, method lineage, geometry, and energy identity.
- `minimum_registry`: global species-to-basin mapping across every trial.
- `reaction`: emitted only for two different registry basins on exactly one PES.

`dft-minima` has two independent admission routes. Low-level candidates bring
their exact trial source and product endpoint; a bounded
`anchor_minima_per_state` ensemble also samples conformer/NCI inputs even when
low-level discovery returns no product. `minimum-registry` remains the only
authority for final DFT basin identity.

`connect-minima` does not guess product coordinates. It checks distinct,
same-PES DFT basin pairs against existing product-free trial coordinates in
both directions. A pair may be promoted as a conservative coordinate
reorganization without claiming a complete bond rearrangement.

## Scientific gates

Reactant and product acceptance requires:

- identical elemental composition, charge, multiplicity, atom order, and method
  lineage;
- real DFT optimization and complete frequency analysis;
- `n_imag = 0` for both endpoints;
- assignment to two distinct, non-ambiguous registry basins.

A transition state requires a converged stationary point, exactly one
significant imaginary frequency, alignment of that mode with the declared
reaction coordinate or path tangent, and an energy above both validated minima.
An activation barrier is final only after forward and backward IRC endpoints are
unconstrainedly reoptimized and matched to the declared registry basins.

The classifier distinguishes a validated elementary step, same-basin
relaxation, effectively barrierless behavior within the declared energy
resolution, and a multistep path containing validated intermediate basins.
Unresolved evidence remains unresolved; it is not forced into a class.

## Bounded amine/HF exploration

`build-complexes` accepts:

```yaml
compositions:
  - composition_id: tma_hf2
    max_seeds: 4
    max_nci_complexes: 4
    components:
      - selector: {artifact_id: spc_tma}
        count: 1
      - selector: {artifact_id: spc_hf}
        count: 2
```

Selectors and counts are expanded up to the declared `max_components` by the same
code.  Multiple bounded distance/direction/twist seeds are passed to CREST
`--nci`.  CREST NCI never falls back to RDKit.  The production TMA(HF)2 config
limits discovery to twelve trials distributed across four NCI conformers and
DFT to at most three low-level-diverse products plus each product's exact trial
source.  These limits control cost, not chemical likelihood.  Obvious low-level
duplicates may be removed before DFT, but only the frequency-validated DFT
registry decides final basin identity.

If all TMA(HF)2 candidates collapse to the same DFT basin or no product minimum
is found, the correct result is `same_basin` or `invalid_endpoints`; NEB, TS,
IRC, and a barrier are then not fabricated.

## Implemented validation

The HCN benchmark has been executed with ReaDuct 6.1.0/GFN2.  NT2 produced a
drive trajectory; its direct GFN2 saddle optimizations failed and were retained
as negative evidence.  Removing the drive and optimizing the final frame gave
an HNC endpoint (N-H about 1.00 A and C-H about 2.16 A).  NWChem 7.2.3
PBE0-D3/def2-SVPD opt/freq then confirmed HCN and HNC as different minima with
zero imaginary frequencies and identical numerical/method lineage.  Raw files
are under `runs/readuct_hcn_validation` and
`runs/hcn_readuct_nwchem_validation`.

The subsequent DFT string/saddle calculation and bidirectional IRC also passed
the scientific gates.  The saddle has one significant imaginary frequency
(-1131.57 cm-1); its mass-weighted mode overlap with the path tangent is 0.547.
pysisyphus supplied the two IRC branches and independent NWChem optimizations
of their terminal structures returned the declared HCN and HNC basins (maximum
permutation RMSD 2.27e-5 A).  The forward HCN -> HNC barriers at 298.15 K are
46.6672 kcal/mol (electronic), 43.2517 kcal/mol (E+ZPE), 43.4305 kcal/mol
(enthalpy), and 42.2296 kcal/mol (Gibbs).  Evidence is retained in
`runs/hcn_saddle_v3`, `runs/hcn_irc_v5`, and `runs/hcn_thermo_v6_complete`.
GoodVibes 4.3 independently reproduced the 298.15 K Gibbs barrier as
42.2809 kcal/mol, within about 0.0513 kcal/mol of the internal value. The
hash-bound compatibility copy and raw unchanged NWChem source are retained in
`runs/hcn_goodvibes_v4`; production thermochemistry is therefore enabled for
this benchmark.

The bounded TMA(HF)2 application was then rerun across four CREST NCI
conformers with three trials per conformer.  Twenty-four ReaDuct attempts,
including twelve association-targeted AFIR calculations for mixed
association/dissociation trials, produced no chemically distinct low-level
product.  Eleven endpoints were only equivalent-H/F permutations, eight did
not realize the requested bond changes, and five produced no NT2 TS guess.
The low-level result still records `no_distinct_low_level_product`. A separate
DFT anchor audit subsequently registered two distinct same-PES TMA(HF)2 minima.
`connect-minima` linked them conservatively as coordinate reorganization and
`reaction-plan` authorized a path calculation. Evidence is under
`runs/tma_hf2_minima_connection_v2`; the path branch is
`runs/tma_hf2_path_v1`. No TMA(HF)2 barrier is final until TS/frequency/IRC gates
all pass.

## Remaining scientific work

- Add regression fixtures from completed ReaDuct/NWChem/pysisyphus raw outputs
  without replacing live production capability checks.
- Add method/basis sensitivity checks before treating the HCN number as a
  reference value.
- Calibrate method and finite-cluster limitations before comparing amine gases
  for a process decision.

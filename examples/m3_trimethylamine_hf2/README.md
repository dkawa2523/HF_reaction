# TMA(HF)2 M3 endpoint pair

This pair is a molecular, gas-phase cluster benchmark for an actual amine/HF
candidate. It contains trimethylamine and two HF molecules; no surface or wet
etch transport model is implied.

Both source structures were independently accepted as zero-imaginary-mode
PBE0-D3(0)/def2-SVPD minima in
`runs/hf2_endpoint_reassessment_002/06_proton-transfer-endpoints/cases/02_trimethylamine_hf2`.
The neutral and shared-proton minima have identity-invariant proton coordinates
of -0.21975 and -0.13541 A, respectively. The latter lies close to the
classification boundary, so the path is called proton reorganization rather
than a completed neutral-to-ion-pair proton transfer.

The shared-proton XYZ was reordered only within symmetry-equivalent methyl
groups and rigidly Kabsch-aligned to the neutral endpoint. The reviewed product
index used for each reactant index was:

```text
[2, 1, 3, 0, 9, 8, 7, 12, 10, 11, 6, 5, 4, 13, 14, 15, 16]
```

The reactive atoms remain fixed: N=1, transfer H=13, accepting F=14,
spectator H=15, terminal F=16. This avoids an artificial path caused by
permuting equivalent methyl groups while preserving the declared chemistry.

The completed 5-bead validation attempt is retained under
`runs/m3_real_qm_tma_hf2_002`. It relaxed to a monotonic electronic-energy
profile and stalled just above NWChem's maximum-gradient convergence
criterion. The runnable configuration therefore uses 9 beads for the next
near-reactant barrier-resolution attempt; the 5-bead result is not reported as
a transition state.

# Final baseline review

## What is ready

```text
workflow wiring
manifest/artifact data flow
HF association / ion-pair endpoint model
ranking separation into screening/scientific/production
visual review package
operations and HPC planning package
connector audit scaffolding
```

## What is not yet a production scientific claim

```text
offline dummy/fallback energies
fallback TS/IRC paths
internal-only quasi-RRHO approximations without production audit
pseudo Cantera mechanisms
unvalidated method-ladder comparisons
```

## Main risk now controlled

The main Phase 12 risk was that dummy/fallback screening values could be mistaken for production ranking values.  The finalized baseline controls this by:

```text
confidence caps
scientific_rank_eligible
production_rank_eligible
rank_screening / rank_scientific / rank_production separation
science-status and next-actions CLI
ranking_summary.json
```

## Recommended next scientific milestone

Run a small reference campaign with real xTB/CREST, ORCA DFT, TS/IRC, GoodVibes, and high-level SP for a curated HF/amine set.  Use the result to establish method bias and uncertainty before screening large proprietary SDF libraries.

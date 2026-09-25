# hfauto technical overview

`hfauto` is a manifest-based platform for molecular complex and reaction-path
discovery. Its scientific scope is a molecular potential-energy surface with
explicit composition, charge, multiplicity, and electronic-structure method.
It does not model transport, liquid hydrodynamics, a solid surface, or wet-etch
equipment.

## Canonical workflow

```text
ingest -> enumerate-states -> conformers -> build-complexes -> preopt
       -> generate-reactions -> explore-reactions -> dft-minima
       -> minimum-registry -> connect-minima -> reaction-plan
       -> NEB/string -> TS/frequency -> IRC -> reaction-classify
       -> thermo -> basin-populations -> reaction-rank -> discovery-audit
```

Discovery trials do not assume a product. ReaDuct NT2/AFIR may propose a
candidate, but AFIR-biased energies are never barriers. A candidate becomes a
`reaction` only after real opt/freq calculations place both endpoints in
different basins on the same configured PES.

An activation barrier is production-eligible only when all of the following
hold:

- both endpoints are frequency-validated minima (`n_imag = 0`);
- the TS has exactly one significant imaginary mode aligned with the path or
  declared coordinate;
- bidirectional IRC endpoints independently reoptimize into the declared
  registry basins;
- thermochemistry records its temperature, standard state, frequency scaling,
  and low-frequency treatment.

An interrupted path may provide a saddle geometry seed after endpoint and
trajectory validation. Its unconverged path energy remains nonpublishable.
Finite search budgets are reported by `discovery-audit`; zero observed products
never means that no reaction exists.

## Code boundaries

- `hfauto/chemistry`: pure chemical and numerical decisions.
- `hfauto/backends`: external-program adapters and raw evidence parsing.
- `hfauto/stages`: one manifest transformation per workflow stage.
- `hfauto/workflow`: shared orchestration and reaction-state transitions.
- `hfauto_ops`: scheduler-neutral planning, retry, reuse, and run operations.
- `hfauto_viz`: read-only visualization and reporting.

The artifact contract is centralized in `hfauto/core/artifact_types.py`; the
lazy stage registry is `hfauto/stages/registry.py`. Production failures create
typed failure artifacts. Dummy or fallback calculations cannot satisfy science
gates.

For maintained details, see `docs/current/index.md`, especially
`reaction_discovery_platform.md`, `stage_io_reference.md`,
`developer_architecture.md`, and `science_gate.md`.

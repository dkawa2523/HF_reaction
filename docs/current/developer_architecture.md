# Developer architecture

Each stage has one responsibility and one manifest boundary:

```text
chemistry/reaction_trials.py       product-free chemical proposals
stages/build_complexes.py          multicomponent assembly and definition checkpoints
backends/reaction_discovery/       ReaDuct execution and raw evidence
stages/preopt.py                   unbiased relaxation and connectivity-retention gate
stages/relaxation_discovery.py     topology-changing relaxation -> standard candidate
stages/generate_reactions.py       trial artifact selection and budgets
stages/explore_reactions.py        attempts and low-level candidates
workflow/discovery_coverage.py     finite-budget evidence accounting
stages/discovery_audit.py          coverage and cross-budget comparison publication
stages/dft_minima.py               real opt/freq minimum promotion
stages/minimum_registry.py         the only basin deduplication authority
chemistry/minimum_connections.py   continuous trial-coordinate comparison
stages/connect_minima.py           distinct-basin candidate publication
stages/reaction_plan.py            candidate classification/reaction promotion
stages/recover_path.py             interrupted path -> geometry-only saddle seed
backends/ts/                       NEB, string, saddle execution
chemistry/path_initialization.py   connectivity-preserving initial path geometry
backends/ts/nwchem_path_support.py endpoint/path validation shared by NEB/string
stages/ts_search.py                budgets, checkpoints, and TS orchestration
stages/irc.py                      independent connectivity validation
stages/sp.py                       fixed-geometry energy execution
stages/method_panel.py             method sensitivity only; no PES promotion
stages/thermo_sensitivity.py       external GoodVibes cutoff sensitivity only
stages/basin_populations.py        conditional populations over observed basins
stages/reaction_rank.py            evidence and uncertainty-aware comparison
chemistry/xyz.py                   generic XYZ model and I/O
```

Backend exceptions become failure artifacts.  Stages do not parse external
program output, backends do not select chemical candidates, and reporting does
not alter scientific artifacts.  Add a reaction-discovery backend by returning
`DiscoveryResult`; do not copy stage orchestration or create another basin
registry.

Long external commands stream stdout/stderr through `core/executables.py`.
Backends own their raw work directories and stages allocate a new retry
directory instead of overwriting interrupted path evidence. `build-complexes`
checkpoints after each definition and reuses it only when the definition,
source-geometry hashes, and CREST settings have the same request fingerprint.
Species-scoped CREST IDs cannot reuse a fragment `mol_id`. A second manifest
guard rejects any cross-lineage artifact-ID collision.

Production preoptimization compares atom-mapped covalent graphs before and
after relaxation. A changed graph is not written back as a revision of the
requested state. `relaxation-discovery` instead normalizes the endpoint into
the same trial/attempt/candidate contract as ReaDuct and labels it only as a
low-level effectively-barrierless hypothesis. CREST-specific safety-stop text
is parsed by the CREST backend, not by the stage.

The standard HCN, TMA(HF)2, and amine/HF pipelines all use this same route:
real xTB is required, each job is fingerprint-checkpointed, retained states may
be revised, and topology-changing relaxations are normalized before driven
trial generation.

`dft-minima` owns the expensive-candidate budget. It derives a molecular-surface
key from elemental composition, charge, and multiplicity, orders low-level
energies only inside that key, and distributes both candidate and anchor
budgets round-robin across keys. Selection logic therefore remains valid when
one manifest contains differently sized amines, HF counts, or explicit water.
The low-level pilot and its DFT follow-up use separate pipeline files so the
discovery campaign can be repeated without accidentally launching DFT/TS jobs.

Stages that consume an ensemble may declare `required_source_data_keys`.
Amine/HF and TMA(HF)2 preoptimization require
`crest_nci_conformer_id`, so deterministic placement seeds remain provenance
but are not redundantly optimized or mistaken for sampled NCI complexes.

`preopt` uses the same restart rule at job granularity: engine/method plus input
XYZ hash must match, and the saved final-XYZ hash must still verify. A retry
always receives a new `attempt_NN` directory. Production promotion additionally
requires an explicitly real, non-dummy optimization.

An unconverged NEB/string profile may supply a saddle *geometry* only when the
endpoint geometry and one internal maximum are validated. Interrupted paths
also require method-identity evidence from the rendered input, program version,
and requested dispersion, while normal termination remains explicitly false.
Their path energy is marked non-publishable. The independently optimized
first-order saddle and IRC, not path convergence, are the final chemical
evidence.

NEB and string histories are normalized by `chemistry/path_diagnostics.py`.
String stagnation is based on the reported convergence tolerance and the
best-XRMS improvement across two complete windows; it is not inferred from one
oscillating step. `recover-path` dispatches preserved NEB and string histories
through the same backend-neutral boundary.

Double-ended paths use a supplied `XYZ_PATH` when available.  The chemistry
layer can generate one for an acyclic dihedral rotation or a single-atom
transfer across bonded donor/acceptor centers.  It rejects nonbonded collisions,
unintended covalent-bond compression/stretching, cyclic dihedral cuts, and
ambiguous transfers. Unsupported coordinates remain explicitly unsupported;
they are not mislabeled as chemistry-aware interpolation.

Do not create a `reaction` from a guessed product, a biased AFIR endpoint, or a
geometry-only distinction.  Do not add reaction-family-specific stage routes.
New chemistry features should emit the existing `ReactionTrialRecord`.

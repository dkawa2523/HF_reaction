> 移行中(ブランチ refactor/fundamental-2026-09): 本 README は Wave 9 で書き直すまで最新ではない。設計は docs/current/refactor_design.md。
# hfauto molecular reaction discovery

`hfauto` is a manifest-based workflow for discovering molecular products,
validating distinct minima, and calculating reaction paths.  Its scientific
scope is gas-phase molecules and non-covalent molecular complexes.  A wet-etch
reactor or solid-surface physical model is intentionally outside this package.

The production workflow is:

```text
ingest -> enumerate-states -> conformers -> build-complexes
       -> preopt -> relaxation-discovery
       -> generate-reactions -> explore-reactions
       -> dft-minima -> minimum-registry -> connect-minima -> reaction-plan
       -> ts-search -> irc -> reaction-classify -> thermo
       -> basin-populations -> reaction-rank
       -> discovery-audit
```

`method-panel` and `thermo-sensitivity` are optional validation branches. They
quantify fixed-geometry electronic/grid sensitivity and external GoodVibes
low-frequency-cutoff sensitivity without changing the production reaction.

The central rule is that discovery does not create a reaction.  ReaDuct AFIR
or NT2 first publishes a `reaction_candidate`.  A `reaction` is created only
after both endpoints are real opt/freq minima on one configured PES and the
global `minimum_registry` assigns them to different basins.  Same-basin pairs,
missing minima, ambiguous basin identity, and method-lineage mismatches never
enter NEB/TS.

## Production environment

The validated environment is WSL/Linux, Python 3.12, CREST, xTB, NWChem 7.2.3,
pysisyphus 1.0, ReaDuct 6.1.0, and the SCINE xTB wrapper:

```bash
uv run --extra production --frozen hfauto doctor \
  --config configs/pipelines/molecular_reaction_discovery_nwchem.yaml --strict
```

ReaDuct wheels are Linux/Python-3.12-specific in this project.  A missing or
failed external program creates a failure artifact; no dummy calculation is
accepted as scientific evidence.

## Reference pipelines

HCN -> HNC is the compact positive benchmark:

```bash
hfauto pipeline \
  --config configs/pipelines/molecular_reaction_discovery_nwchem.yaml \
  --run-id hcn_validation
```

TMA(HF)2 uses four bounded CREST NCI complexes, twelve reaction trials shared
across them, low-level-diverse product endpoints, and a bounded DFT anchor
ensemble that remains active when low-level discovery finds no product:

```bash
hfauto pipeline \
  --config configs/pipelines/m3_trimethylamine_hf2_nwchem.yaml \
  --run-id tma_hf2_validation
```

The broader amine/HF pilot deliberately separates bounded low-level discovery
from expensive evidence promotion. After the first command reaches its final
`discovery-audit`, the second pipeline selects candidates round-robin across
composition/charge/multiplicity surfaces and can be stopped at
`reaction-plan` before any TS work is authorized:

```bash
hfauto pipeline \
  --config configs/pipelines/amine_hf_pilot_hf1_hf3.yaml \
  --run-id amine_hf_pilot

hfauto pipeline \
  --config configs/pipelines/amine_hf_pilot_dft_followup.yaml \
  --run-id amine_hf_pilot_dft \
  --from dft-minima --to reaction-plan \
  --start-manifest runs/amine_hf_pilot/08_discovery-audit/manifest.json
```

AFIR artificial-force energies are never reported as barriers.  A production
activation barrier requires a first-order saddle with exactly one significant
imaginary mode and bidirectional IRC endpoints that reoptimize into the
declared reactant and product registry basins.

The completed HCN/HNC benchmark also passes an external GoodVibes 4.3
thermochemistry cross-check. `configs/pipelines/hcn_method_sensitivity.yaml`
adds a separate fixed-geometry electronic method panel; its values are labeled
as sensitivity and do not revalidate stationary points on another PES.
`configs/pipelines/hcn_numerical_sensitivity.yaml` adds a fine/xfine grid panel
and an external-GoodVibes 50/100/150 cm-1 cutoff panel.

Long TS searches checkpoint after every attempt and support explicit reaction,
attempt, and wall-time budgets. `recover-path` may recover a geometry seed from
a preserved interrupted path, but it never publishes the unconverged path
energy as a barrier.

See [the canonical workflow design](docs/current/reaction_discovery_platform.md)
and [stage I/O](docs/current/stage_io_reference.md).

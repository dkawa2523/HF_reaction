# Production environment

The supported production path is Linux or WSL with Python 3.12:

```text
CREST NCI -> xTB preoptimization -> ReaDuct NT2/AFIR
-> NWChem 7.2.3 PBE0-D3/def2-SVPD minima and TS
-> pysisyphus 1.0 + QCEngine + NWChem bidirectional IRC
```

Install/synchronize from the locked project environment. The `production`
extra pins ReaDuct 6.1.0 and its SCINE dependencies on Linux Python 3.12.
CREST, xTB and NWChem may be on `PATH` or declared by an explicit executable
path in the stage that owns the program. Production preflight validates those
same explicit paths, including the MPI launcher.

Validate the installed capabilities with the real HCN/HNC workflow:

```bash
uv run --extra production --frozen hfauto doctor \
  --config configs/pipelines/molecular_reaction_discovery_nwchem.yaml \
  --strict
```

The JSON report exposes both `required_python` and `required_external`. A
ReaDuct discovery stage requires `scine_readuct` and `scine_xtb_wrapper` in the
former; CREST and xTB remain external executable requirements.

Using a plain pre-existing `.venv` without the `production` extra is not a
capability check; `uv run --extra production --frozen` is the canonical,
reproducible invocation.

The bounded TMA(HF)2 application is defined in
`configs/pipelines/m3_trimethylamine_hf2_nwchem.yaml`.
The multi-formula amine/HF campaign is split deliberately: bounded discovery
uses `configs/pipelines/amine_hf_pilot_hf1_hf3.yaml`, and checkpointed DFT/TS
promotion uses `configs/pipelines/amine_hf_pilot_dft_followup.yaml` with a
completed discovery manifest as input.

Production rules are fail-closed:

- artificial AFIR energies are never reported as barriers;
- absolute low-level energies are never compared across different elemental
  compositions, charges, or multiplicities when allocating the DFT budget;
- dummy fallback is disabled;
- only frequency-validated DFT minima enter the minimum registry;
- only two distinct minima on one numerical PES are promoted to a reaction;
- a TS needs exactly one significant imaginary mode aligned with the path;
- connectivity is accepted only after bidirectional IRC and endpoint
  reoptimization recover the registered reactant/product basins.

Raw external-program inputs, outputs, trajectories and hashes remain under
`runs/` and are not cleaned by the workflow.

External stdout and stderr are streamed directly into each backend work
directory. Long NEB/IRC jobs therefore expose live raw progress, do not retain
the complete output in process memory, and preserve output produced before a
timeout.

Production DFT minima are checkpointed after each completed external job.  On
an interrupted rerun, a result is reused only when the input geometry, source
artifact, engine/method request, real-QM status, successful command record, and
all recorded raw-file SHA-256 hashes match.  An incomplete or failed retry is
written to a new `attempt_NN` directory so that the earlier evidence is not
overwritten.

Path/TS attempts apply the same non-overwrite rule to interrupted raw work
directories. A rerun receives a `_retry_NN` directory; partial NWChem databases,
trajectories, and inputs remain audit evidence rather than an implicit restart.

`ts-search` additionally checkpoints after each completed attempt and skips an
integrity-valid frequency-validated TS when resumed. `max_reactions`,
`max_total_attempts`, and `max_stage_walltime_s` defer work with an
`execution_budget_exhausted` artifact; this category is operational and is not
evidence that a reaction or saddle does not exist.

For a preserved interrupted path, `recover-path` validates the endpoint atom
mapping and selects a unique internal maximum as a geometry-only seed. It always
sets `path_energy_publishable=false`; only the later saddle/frequency/IRC chain
can publish an activation barrier.

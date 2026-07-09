# Extension points

The system is intentionally split into three packages.

```text
hfauto      calculation workflow and ranking
hfauto_viz  read-only visualization
hfauto_ops  read-only operations
```

## Add a new QM backend

Implement a backend under `hfauto/backends/qm/` and register it in `hfauto/backends/registry.py` or the relevant stage config.  Keep the stage output as `calculation` artifacts with the same core fields:

```text
electronic_energy_hartree
zpe_hartree
gibbs_298K_hartree
frequencies_cm1
n_imag
scf_converged
geometry_converged
fallback_dummy
real_<engine>_executed
```

## Add a new reaction family

Add a builder under `hfauto/chemistry/` and emit normal `species` and `reaction` artifacts.  Avoid special one-off downstream logic.  A new reaction family should define:

```text
reaction_type
reactant_species_id
product_species_id
optional ts_species_id
reaction_coordinate definition
expected endpoint graph check
```

## Add a new visualization

Add a renderer under `hfauto_viz/renderers/` and register it from the relevant dossier/report function.  Visualization should read existing artifacts and write a `visualization` entry in `viz_manifest.json`.

## Add a new operations action

Add it under `hfauto_ops/`.  Operations should not rewrite scientific artifacts.  If a workflow needs to reuse or replace a calculation, create a review plan first rather than mutating results silently.

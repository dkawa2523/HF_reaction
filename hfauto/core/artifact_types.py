from __future__ import annotations

# Central, lightweight names used in stage contracts and documentation.
# The workflow still accepts string artifact_type values for backward compatibility.

MOLECULE = "molecule"
MOLECULE_ENRICHED = "molecule_enriched"
SITE = "site"
CONFORMER = "conformer"
SPECIES = "species"
SPECIES_PREOPT = "species_preopt"
SPECIES_OPTIMIZED = "species_optimized"
REACTION = "reaction"
REACTION_VALIDATED = "reaction_validated"
REACTION_PATH_VALIDATED = "reaction_path_validated"
CALCULATION = "calculation"
THERMO = "thermo"
KINETICS = "kinetics"
DESCRIPTOR = "descriptor"
RANKING = "ranking"
TABLE = "table"
VISUALIZATION = "visualization"
OPS = "ops"
FAILURE = "failure"

# Human-oriented contract summary.  This is intentionally not enforced at runtime
# to avoid introducing a heavy contract layer.
STAGE_CONTRACTS = {
    "ingest": {"in": ["SDF"], "out": [MOLECULE]},
    "enrich": {"in": [MOLECULE], "out": [MOLECULE_ENRICHED]},
    "detect-sites": {"in": [MOLECULE_ENRICHED, MOLECULE], "out": [SITE]},
    "conformers": {"in": [MOLECULE, SITE], "out": [CONFORMER]},
    "build-hf": {"in": [CONFORMER, SITE], "out": [SPECIES, REACTION]},
    "preopt": {"in": [SPECIES], "out": [SPECIES_PREOPT, CALCULATION]},
    "dft-minima": {"in": [SPECIES_PREOPT, SPECIES], "out": [SPECIES_OPTIMIZED, CALCULATION]},
    "ts-search": {"in": [REACTION, SPECIES_OPTIMIZED, SPECIES_PREOPT], "out": [REACTION_VALIDATED, SPECIES, CALCULATION]},
    "irc": {"in": [REACTION_VALIDATED], "out": [REACTION_PATH_VALIDATED]},
    "sp": {"in": [SPECIES_OPTIMIZED, SPECIES_PREOPT, SPECIES], "out": [CALCULATION]},
    "thermo": {"in": [CALCULATION, REACTION_PATH_VALIDATED], "out": [THERMO]},
    "kinetics": {"in": [THERMO], "out": [KINETICS]},
    "calibrate": {"in": [MOLECULE_ENRICHED, THERMO], "out": [TABLE, "method_validation"]},
    "rank": {"in": [THERMO, KINETICS, DESCRIPTOR], "out": [RANKING, TABLE]},
    "viz": {"in": [RANKING, TABLE, SPECIES, REACTION], "out": [VISUALIZATION]},
    "ops": {"in": ["latest manifest"], "out": [OPS, TABLE]},
}

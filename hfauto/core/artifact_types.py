from __future__ import annotations

# Central, lightweight names used in stage contracts and documentation.
# Artifact types remain strings because manifests are language-neutral JSON.

MOLECULE = "molecule"
MOLECULE_ENRICHED = "molecule_enriched"
MOLECULE_STATE = "molecule_state"
SITE = "site"
CONFORMER = "conformer"
SPECIES = "species"
SPECIES_PREOPT = "species_preopt"
PREOPT_GEOMETRY = "preopt_geometry"
SPECIES_OPTIMIZED = "species_optimized"
REACTION = "reaction"
REACTION_TRIAL = "reaction_trial"
REACTION_DISCOVERY_ATTEMPT = "reaction_discovery_attempt"
DISCOVERY_COVERAGE = "discovery_coverage"
REACTION_CANDIDATE = "reaction_candidate"
MINIMUM_BASIN = "minimum_basin"
MINIMUM_REGISTRY = "minimum_registry"
REACTION_VALIDATED = "reaction_validated"
REACTION_PATH_VALIDATED = "reaction_path_validated"
REACTION_SCAN = "reaction_scan"
REACTION_PATH = "reaction_path"
REACTION_CASE = "reaction_case"
REACTION_CLASSIFICATION = "reaction_classification"
REACTION_SEGMENTATION_PLAN = "reaction_segmentation_plan"
PATH_ATTEMPT = "path_attempt"
PATH_ENSEMBLE_MEMBER = "path_ensemble_member"
PATH_ENSEMBLE_ASSESSMENT = "path_ensemble_assessment"
MINIMUM_MODE_FOLLOWING_PLAN = "minimum_mode_following_plan"
MINIMUM_MODE_FOLLOWING_ASSESSMENT = "minimum_mode_following_assessment"
ENDPOINT_PAIR_SELECTION = "endpoint_pair_selection"
ENDPOINT_SEED_ENSEMBLE_PLAN = "endpoint_seed_ensemble_plan"
ENDPOINT_SEED_SELECTION = "endpoint_seed_selection"
SADDLE_SEED_SCAN = "saddle_seed_scan"
SADDLE_ATTEMPT = "saddle_attempt"
CALCULATION = "calculation"
THERMO = "thermo"
KINETICS = "kinetics"
DESCRIPTOR = "descriptor"
RANKING = "ranking"
TABLE = "table"
METHOD_UNCERTAINTY = "method_uncertainty"
BASIN_POPULATION = "basin_population"
VISUALIZATION = "visualization"
OPS = "ops"
FAILURE = "failure"

# Human-oriented contract summary.  This is intentionally not enforced at runtime
# to avoid introducing a heavy contract layer.
STAGE_CONTRACTS = {
    "ingest": {"in": ["SDF"], "out": [MOLECULE]},
    "enrich": {"in": [MOLECULE], "out": [MOLECULE_ENRICHED]},
    "enumerate-states": {"in": [MOLECULE_ENRICHED, MOLECULE], "out": [MOLECULE_STATE]},
    "detect-sites": {"in": [MOLECULE_ENRICHED, MOLECULE], "out": [SITE]},
    "conformers": {"in": [MOLECULE_STATE, MOLECULE, SITE], "out": [CONFORMER]},
    "generate-reactions": {
        "in": [SPECIES, SPECIES_PREOPT],
        "out": [REACTION_TRIAL],
    },
    "explore-reactions": {
        "in": [REACTION_TRIAL, SPECIES_PREOPT],
        "out": [REACTION_DISCOVERY_ATTEMPT, REACTION_CANDIDATE, SPECIES_PREOPT],
    },
    "discovery-audit": {
        "in": [
            REACTION_TRIAL,
            REACTION_DISCOVERY_ATTEMPT,
            REACTION_CANDIDATE,
            MINIMUM_REGISTRY,
            REACTION_VALIDATED,
            REACTION_PATH_VALIDATED,
        ],
        "out": [DISCOVERY_COVERAGE],
    },
    "minimum-registry": {
        "in": [REACTION_CANDIDATE, SPECIES_OPTIMIZED, CALCULATION],
        "out": [MINIMUM_BASIN, MINIMUM_REGISTRY],
    },
    "connect-minima": {
        "in": [MINIMUM_REGISTRY, MINIMUM_BASIN, SPECIES_OPTIMIZED, REACTION_TRIAL],
        "out": [REACTION_CANDIDATE],
    },
    "build-complexes": {
        "in": [CONFORMER, SPECIES],
        "out": [SPECIES],
    },
    "preopt": {
        "in": [SPECIES],
        "out": [SPECIES_PREOPT, PREOPT_GEOMETRY, CALCULATION],
    },
    "relaxation-discovery": {
        "in": [PREOPT_GEOMETRY, CONFORMER, SPECIES],
        "out": [
            REACTION_TRIAL,
            REACTION_DISCOVERY_ATTEMPT,
            REACTION_CANDIDATE,
            SPECIES_PREOPT,
        ],
    },
    "dft-minima": {
        "in": [SPECIES_PREOPT, SPECIES, REACTION_CANDIDATE],
        "out": [SPECIES_OPTIMIZED, CALCULATION],
    },
    "method-panel": {
        "in": [REACTION_VALIDATED, REACTION_PATH_VALIDATED, CALCULATION],
        "out": [TABLE, METHOD_UNCERTAINTY],
    },
    "basin-populations": {
        "in": [MINIMUM_BASIN, "species_thermo"],
        "out": [BASIN_POPULATION, TABLE],
    },
    "minimum-mode-follow": {
        "in": [CALCULATION, SPECIES],
        "out": [MINIMUM_MODE_FOLLOWING_PLAN, SPECIES],
    },
    "minimum-mode-assess": {
        "in": [MINIMUM_MODE_FOLLOWING_PLAN, CALCULATION, SPECIES_OPTIMIZED],
        "out": [MINIMUM_MODE_FOLLOWING_ASSESSMENT],
    },
    "endpoint-seeds": {
        "in": [REACTION, ENDPOINT_PAIR_SELECTION, SPECIES_OPTIMIZED, CALCULATION],
        "out": [ENDPOINT_SEED_ENSEMBLE_PLAN, SPECIES],
    },
    "endpoint-seed-screen": {
        "in": [ENDPOINT_SEED_ENSEMBLE_PLAN, SPECIES_PREOPT, CALCULATION],
        "out": [ENDPOINT_SEED_SELECTION],
    },
    "reaction-plan": {
        "in": [
            REACTION_CANDIDATE,
            MINIMUM_REGISTRY,
            MINIMUM_BASIN,
            SPECIES_OPTIMIZED,
            PATH_ATTEMPT,
            SADDLE_ATTEMPT,
        ],
        "out": [REACTION, "basin_pair_assessment", REACTION_CASE],
    },
    "recover-path": {
        "in": [REACTION, SPECIES_OPTIMIZED, REACTION_CASE],
        "out": [
            REACTION_PATH,
            PATH_ATTEMPT,
            SADDLE_ATTEMPT,
            SPECIES,
            REACTION_CASE,
        ],
    },
    "reaction-segments": {
        "in": [REACTION, REACTION_CLASSIFICATION, SPECIES_OPTIMIZED],
        "out": [REACTION_SEGMENTATION_PLAN, REACTION],
    },
    "ts-search": {
        "in": [
            REACTION,
            SPECIES_OPTIMIZED,
            SPECIES_PREOPT,
        ],
        "out": [
            REACTION_PATH,
            PATH_ATTEMPT,
            SADDLE_SEED_SCAN,
            SADDLE_ATTEMPT,
            REACTION_CASE,
            REACTION_VALIDATED,
            SPECIES,
            CALCULATION,
        ],
    },
    "path-intermediates": {
        "in": [
            REACTION,
            REACTION_PATH,
            PATH_ATTEMPT,
            PATH_ENSEMBLE_ASSESSMENT,
        ],
        "out": [SPECIES],
    },
    "path-ensemble": {
        "in": [REACTION, SPECIES_OPTIMIZED, REACTION_PATH, PATH_ATTEMPT],
        "out": [
            REACTION_PATH,
            PATH_ATTEMPT,
            PATH_ENSEMBLE_MEMBER,
        ],
    },
    "path-ensemble-assess": {
        "in": [
            REACTION,
            SPECIES_OPTIMIZED,
            PATH_ENSEMBLE_MEMBER,
            REACTION_PATH,
            PATH_ATTEMPT,
        ],
        "out": [PATH_ENSEMBLE_MEMBER, PATH_ENSEMBLE_ASSESSMENT],
    },
    "irc": {"in": [REACTION_VALIDATED], "out": [REACTION_PATH_VALIDATED]},
    "reaction-classify": {
        "in": [
            REACTION,
            SPECIES_OPTIMIZED,
            REACTION_PATH,
            PATH_ATTEMPT,
            REACTION_VALIDATED,
            REACTION_PATH_VALIDATED,
            PATH_ENSEMBLE_ASSESSMENT,
            MINIMUM_MODE_FOLLOWING_ASSESSMENT,
        ],
        "out": [REACTION_CLASSIFICATION],
    },
    "sp": {"in": [SPECIES_OPTIMIZED, SPECIES_PREOPT, SPECIES], "out": [CALCULATION]},
    "thermo": {"in": [CALCULATION, REACTION_PATH_VALIDATED], "out": [THERMO]},
    "kinetics": {"in": [THERMO], "out": [KINETICS]},
    "calibrate": {"in": [MOLECULE_ENRICHED, THERMO], "out": [TABLE, "method_validation"]},
    "rank": {"in": [THERMO, KINETICS, DESCRIPTOR], "out": [RANKING, TABLE]},
    "reaction-rank": {
        "in": [THERMO, KINETICS, METHOD_UNCERTAINTY],
        "out": [RANKING],
    },
    "viz": {"in": [RANKING, TABLE, SPECIES, REACTION], "out": [VISUALIZATION]},
    "ops": {"in": ["latest manifest"], "out": [OPS, TABLE]},
}

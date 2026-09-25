"""Lazy stage registry.

The generic workflow no longer imports every optional or domain-specific stage
at process start.  A stage is imported only when it is executed.
"""

from __future__ import annotations

STAGES: dict[str, str] = {
    "ingest": "hfauto.stages.ingest:IngestStage",
    "enumerate-states": "hfauto.stages.enumerate_states:EnumerateStatesStage",
    "detect-sites": "hfauto.stages.detect_sites:DetectSitesStage",
    "conformers": "hfauto.stages.conformers:ConformersStage",
    "build-complexes": "hfauto.stages.build_complexes:BuildComplexesStage",
    "generate-reactions": "hfauto.stages.generate_reactions:GenerateReactionsStage",
    "explore-reactions": "hfauto.stages.explore_reactions:ExploreReactionsStage",
    "discovery-audit": "hfauto.stages.discovery_audit:DiscoveryAuditStage",
    "minimum-registry": "hfauto.stages.minimum_registry:MinimumRegistryStage",
    "connect-minima": "hfauto.stages.connect_minima:ConnectMinimaStage",
    "preopt": "hfauto.stages.preopt:PreoptStage",
    "relaxation-discovery": (
        "hfauto.stages.relaxation_discovery:RelaxationDiscoveryStage"
    ),
    "dft-minima": "hfauto.stages.dft_minima:DFTMinimaStage",
    "minimum-mode-follow": (
        "hfauto.stages.minimum_mode_follow:MinimumModeFollowStage"
    ),
    "minimum-mode-assess": (
        "hfauto.stages.minimum_mode_assess:MinimumModeAssessStage"
    ),
    "endpoint-seeds": "hfauto.stages.endpoint_seeds:EndpointSeedsStage",
    "endpoint-seed-screen": (
        "hfauto.stages.endpoint_seed_screen:EndpointSeedScreenStage"
    ),
    "reaction-plan": "hfauto.stages.reaction_plan:ReactionPlanStage",
    "ts-search": "hfauto.stages.ts_search:TSSearchStage",
    "path-ensemble": "hfauto.stages.path_ensemble:PathEnsembleStage",
    "path-ensemble-assess": (
        "hfauto.stages.path_ensemble_assess:PathEnsembleAssessStage"
    ),
    "path-intermediates": (
        "hfauto.stages.path_intermediates:PathIntermediatesStage"
    ),
    "irc": "hfauto.stages.irc:IRCStage",
    "reaction-classify": (
        "hfauto.stages.reaction_classify:ReactionClassifyStage"
    ),
    "reaction-segments": "hfauto.stages.reaction_segments:ReactionSegmentsStage",
    "sp": "hfauto.stages.sp:SinglePointStage",
    "method-panel": "hfauto.stages.method_panel:MethodPanelStage",
    "basin-populations": "hfauto.stages.basin_populations:BasinPopulationsStage",
    "thermo": "hfauto.stages.thermo:ThermoStage",
    "thermo-sensitivity": (
        "hfauto.stages.thermo_sensitivity:ThermoSensitivityStage"
    ),
    "reaction-rank": "hfauto.stages.reaction_rank:ReactionRankStage",
}

def get_stage(name: str):
    from importlib import import_module

    try:
        reference = STAGES[str(name)]
    except KeyError as exc:
        raise KeyError(f"Unknown stage: {name}") from exc
    module_name, class_name = reference.split(":", maxsplit=1)
    stage_type = getattr(import_module(module_name), class_name)
    return stage_type()

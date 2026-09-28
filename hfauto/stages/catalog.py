"""Stage catalog (design §7.4): pipeline ``stage:`` names to lazily imported classes."""

from __future__ import annotations

import importlib

from hfauto.stages.spec import Stage

STAGES: dict[str, str] = {
    "structures": "hfauto.stages.structures:StructuresStage",
    "conformers": "hfauto.stages.conformer_search:ConformersStage",
    "minima": "hfauto.stages.minima:MinimaStage",
    "explore": "hfauto.stages.explore:ExploreStage",
    "reaction-paths": "hfauto.stages.reaction_paths:ReactionPathsStage",
    "sp": "hfauto.stages.single_point:SinglePointStage",
    "thermo": "hfauto.stages.thermochemistry:ThermoStage",
    "report": "hfauto.stages.report:ReportStage",
}

def get(name: str) -> type[Stage]:
    """The stage class for ``name``; an unknown name raises KeyError."""
    if name not in STAGES:
        raise KeyError(f"unknown stage {name!r}; known: {sorted(STAGES)}")
    module, _, attr = STAGES[name].partition(":")
    return getattr(importlib.import_module(module), attr)

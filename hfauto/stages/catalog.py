"""Stage catalog (design §7.4): pipeline ``stage:`` names to lazily imported classes."""

from __future__ import annotations

import importlib
from collections.abc import Iterator
from contextlib import contextmanager

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

_OVERRIDES: dict[str, type[Stage]] = {}


def get(name: str) -> type[Stage]:
    """The stage class for ``name``; an unknown name raises KeyError."""
    if name in _OVERRIDES:
        return _OVERRIDES[name]
    if name not in STAGES:
        raise KeyError(f"unknown stage {name!r}; known: {sorted(STAGES)}")
    module, _, attr = STAGES[name].partition(":")
    return getattr(importlib.import_module(module), attr)


@contextmanager
def override(name: str, cls: type[Stage]) -> Iterator[None]:
    """Test only: resolve ``name`` to ``cls`` inside the block."""
    previous = _OVERRIDES.get(name)
    _OVERRIDES[name] = cls
    try:
        yield
    finally:
        if previous is None:
            del _OVERRIDES[name]
        else:
            _OVERRIDES[name] = previous

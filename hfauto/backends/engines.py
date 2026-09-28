"""Typed lazy engine registry (design §6.2). Only hfauto.pipeline may import it.

The table is final: seven entries, no aliases and no test keys. Tests swap an entry with
``override``; an unknown (capability, name) is a KeyError everywhere.
"""

from __future__ import annotations

import importlib
from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from typing import TYPE_CHECKING, Any

from hfauto.backends.protocols import (
    Capability,
    ConformerEngine,
    DiscoveryEngine,
    Engine,
    EngineFactory,
    PathEngine,
    QMEngine,
    Requirements,
    SaddleRefiner,
)
from hfauto.core.method import EngineSite

if TYPE_CHECKING:
    from hfauto.execution.jobs import JobRunner

_TABLE: dict[Capability, dict[str, str]] = {
    Capability.QM: {"nwchem": "hfauto.backends.nwchem.engine:NWChemEngine",
                    "xtb": "hfauto.backends.xtb:XTBEngine"},
    Capability.PATH: {"nwchem_string": "hfauto.backends.nwchem.engine:NWChemString",
                      "pysis_neb": "hfauto.backends.pysis.engine:PysisNEB"},
    Capability.SADDLE: {"nwchem_saddle": "hfauto.backends.nwchem.engine:NWChemSaddle"},
    Capability.CONFORMERS: {"crest": "hfauto.backends.crest:CRESTEngine"},
    Capability.DISCOVERY: {"readuct": "hfauto.backends.readuct.engine:ReaDuctEngine"},
}

_PROTOCOLS: dict[Capability, type] = {  # runtime_checkable Protocols
    Capability.QM: QMEngine,
    Capability.PATH: PathEngine,
    Capability.SADDLE: SaddleRefiner,
    Capability.CONFORMERS: ConformerEngine,
    Capability.DISCOVERY: DiscoveryEngine,
}

_OVERRIDES: dict[tuple[Capability, str], Callable[..., Engine]] = {}


def _target(capability: Capability, name: str) -> str:
    try:
        return _TABLE[capability][name]
    except KeyError:
        raise KeyError(f"no {capability!s} engine named {name!r}") from None


def _load(target: str) -> Any:
    module, _, attr = target.partition(":")
    return getattr(importlib.import_module(module), attr)


def create(capability: Capability, name: str, *, jobs: JobRunner, site: EngineSite) -> Engine:
    """Import and build the engine, checking it implements the capability's Protocol."""
    target = _target(capability, name)
    factory: EngineFactory = _OVERRIDES.get((Capability(capability), name)) or _load(target)
    engine = factory(jobs=jobs, site=site)
    protocol = _PROTOCOLS[Capability(capability)]
    if not isinstance(engine, protocol):
        raise TypeError(f"{name!r} does not implement {protocol.__name__}")
    return engine


def requirements_for(selection: Mapping[Capability, Sequence[str]]) -> dict[str, Requirements]:
    """Requirements per engine name; overridden (test) engines require nothing."""
    out: dict[str, Requirements] = {}
    for capability, names in selection.items():
        for name in names:
            target = _target(capability, name)
            if (Capability(capability), name) in _OVERRIDES:
                out[name] = Requirements()
            else:
                out[name] = _load(target).requirements()
    return out


@contextmanager
def override(
    capability: Capability, name: str, factory: Callable[..., Engine]
) -> Iterator[None]:
    """Test only: ``create`` calls ``factory(jobs=..., site=...)`` for this existing key."""
    _target(capability, name)
    key = (Capability(capability), name)
    previous = _OVERRIDES.get(key)
    _OVERRIDES[key] = factory
    try:
        yield
    finally:
        if previous is None:
            _OVERRIDES.pop(key, None)
        else:
            _OVERRIDES[key] = previous

"""Capability Protocols and their typed settings / results (design §6.1).

Every engine returns a typed result or a ``Failure``; it never builds artifacts and never
creates other engines. Engines are built as ``Engine(jobs=JobRunner, site=EngineSite)`` by
``hfauto.backends.engines.create``, the only registry.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
from typing import TYPE_CHECKING, ClassVar, Literal, Protocol, runtime_checkable

from pydantic import BaseModel, ConfigDict

from hfauto.chemistry.xyz import Molecule
from hfauto.core.evidence import Evidence, Failure, FileRef, Geometry, PathProfile
from hfauto.core.method import Deadline, EngineSite, MethodSpec
from hfauto.core.records import ReactionTrial

if TYPE_CHECKING:
    from hfauto.execution.jobs import JobRunner

_FROZEN = ConfigDict(frozen=True, extra="forbid")


class Capability(StrEnum):
    QM = "qm"
    PATH = "path"
    SADDLE = "saddle"
    CONFORMERS = "conformers"
    DISCOVERY = "discovery"


@dataclass(frozen=True)
class Requirements:
    executables: tuple[str, ...] = ()  # keys of the site file's executables
    python_modules: tuple[str, ...] = ()  # modules the worker interpreter must import
    version_command: tuple[str, ...] = ()  # preflight runs this to read the version


class ConformerSettings(BaseModel):  # threads: the site's engines.crest.execution
    model_config = _FROZEN
    nci: bool = False
    quick: bool = True
    ewin_kcal: float = 6.0


class ConformerEnsemble(BaseModel):
    model_config = _FROZEN
    kind: Literal["conformers"] = "conformers"
    members: tuple[tuple[Geometry, float | None], ...]  # a missing energy stays None
    topology_stops: tuple[Geometry, ...]


class DiscoverySettings(BaseModel):
    model_config = _FROZEN
    max_scf_iterations: int = 300
    electronic_temperature_K: float = 300.0
    scc_retry_temperature_K: float = 1000.0
    imag_cutoff_cm1: float = 50.0
    timeout_s: float = 600.0


class DiscoveryResult(BaseModel):
    model_config = _FROZEN
    kind: Literal["discovery_result"] = "discovery_result"
    outcome: Literal["product", "negative"]
    reason: str | None
    product: Geometry | None
    ts: Geometry | None
    ts_imag_cm1: float | None
    dE_act_kcal: float | None  # TS - source; kept for a negative that has a TS
    dE_rxn_kcal: float | None
    irc_connected_to_source: bool
    electronic_temperature_K: float
    job_key: str


@runtime_checkable
class Engine(Protocol):
    name: ClassVar[str]

    def requirements(self) -> Requirements:
        """Implemented as a classmethod: engines.requirements_for reads it without an engine."""
        ...

    def supports(self, method: MethodSpec) -> bool: ...


class EngineFactory(Protocol):
    """How engines.create builds an engine (concrete classes are factories themselves)."""

    def __call__(self, *, jobs: JobRunner, site: EngineSite) -> Engine: ...


@runtime_checkable
class QMEngine(Engine, Protocol):
    def energy(
        self, mol: Molecule, method: MethodSpec, *, deadline: Deadline | None = None
    ) -> Evidence | Failure: ...

    def optimize(
        self,
        mol: Molecule,
        method: MethodSpec,
        *,
        # freq Evidence at mol or within 0.5 Å per atom of it (same atom order and frame)
        init_hessian: Evidence | None = None,
        deadline: Deadline | None = None,
    ) -> Evidence | Failure: ...

    def frequencies(
        self,
        mol: Molecule,
        method: MethodSpec,
        *,
        # the opt or saddle at mol whose converged SCF starts this one (same electronic state)
        scf_guess: Evidence | None = None,
        deadline: Deadline | None = None,
    ) -> Evidence | Failure: ...


@runtime_checkable
class PathEngine(Engine, Protocol):
    def find_path(
        self,
        start: Molecule,
        end: Molecule,
        method: MethodSpec,
        *,
        images: int,
        initial_path: FileRef,  # multi-frame xyz from start to end (the relaxed path's start)
        deadline: Deadline | None = None,
    ) -> PathProfile | Failure: ...


@runtime_checkable
class SaddleRefiner(Engine, Protocol):
    def refine(
        self,
        seed: Molecule,
        method: MethodSpec,
        *,
        # freq Evidence (any Level) at seed or within 0.5 Å per atom of it (same atoms, frame)
        hessian: Evidence,
        mode: Sequence[float],  # reaction direction (3N): the only negative initial curvature
        deadline: Deadline | None = None,
    ) -> Evidence | Failure: ...  # maxiter: a Failure whose ``final`` is the last frame


@runtime_checkable
class ConformerEngine(Engine, Protocol):
    def search(
        self,
        mol: Molecule,
        method: MethodSpec,
        settings: ConformerSettings,
        *,
        deadline: Deadline | None = None,
    ) -> ConformerEnsemble | Failure: ...


@runtime_checkable
class DiscoveryEngine(Engine, Protocol):
    def explore(
        self,
        source: Molecule,
        trial: ReactionTrial,
        method: MethodSpec,
        settings: DiscoverySettings,
        *,
        deadline: Deadline | None = None,
    ) -> DiscoveryResult | Failure: ...

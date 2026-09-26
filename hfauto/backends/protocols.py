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
from hfauto.core.method import Deadline, EngineSite, MethodSpec, ThermoSettings
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
    THERMO = "thermo"


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
    notopo_atoms: tuple[int, ...] = ()


class ConformerEnsemble(BaseModel):
    model_config = _FROZEN
    kind: Literal["conformers"] = "conformers"
    members: tuple[tuple[Geometry, float | None], ...]  # a missing energy stays None
    topology_removed: int
    topology_stops: tuple[Geometry, ...]
    version: str
    job_key: str


class DiscoverySettings(BaseModel):
    model_config = _FROZEN
    max_scf_iterations: int = 300
    electronic_temperature_K: float = 300.0
    scc_retry_temperature_K: float = 1000.0
    afir_gamma_kj_mol: float = 125.0
    afir_gamma_retry_kj_mol: float = 300.0
    nt_total_force_norm: float = 0.1
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
    barrier_kj_mol: float | None
    reaction_kj_mol: float | None
    irc_connected_to_source: bool
    electronic_temperature_K: float
    job_key: str


class ThermoResult(BaseModel):
    model_config = _FROZEN
    kind: Literal["thermo_result"] = "thermo_result"
    settings_sha: str
    T_K: float
    G_hartree: float
    H_hartree: float
    E_hartree: float
    zpe_hartree: float
    n_real: int
    S_rot: float
    notes: tuple[str, ...]
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
        self, mol: Molecule, method: MethodSpec, *, deadline: Deadline | None = None
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
        initial_path: FileRef | None = None,
        refine_ts: bool = False,
        deadline: Deadline | None = None,
    ) -> PathProfile | Failure: ...


@runtime_checkable
class SaddleRefiner(Engine, Protocol):
    def refine(
        self,
        seed: Molecule,
        method: MethodSpec,
        *,
        hessian: Evidence,  # freq Evidence whose final geometry is seed (any Level)
        mode_index: int | None = None,  # into imaginary_modes (ascending); engine maps it
        deadline: Deadline | None = None,
    ) -> Evidence | Failure: ...


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


@runtime_checkable
class ThermoEngine(Engine, Protocol):
    def thermo(
        self,
        freq: Evidence,
        settings: Sequence[ThermoSettings],
        *,
        temperatures_K: Sequence[float],
        saddle: bool = False,
        deadline: Deadline | None = None,
    ) -> list[ThermoResult] | Failure:
        """One result per (settings, temperature); modes as chemistry.thermo_frequencies."""
        ...

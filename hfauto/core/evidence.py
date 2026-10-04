"""Typed observations of external jobs (design §5.1).

An ``Evidence`` exists only for a job that finished normally; backends return a
``Failure`` otherwise. Pass/fail decisions live in ``hfauto.chemistry.gates``.
"""

from __future__ import annotations

import hashlib
import json
from enum import StrEnum
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, field_validator, model_validator

_FROZEN = ConfigDict(frozen=True, extra="forbid")

_LEVEL_STRINGS = ("program", "version", "method", "basis", "dispersion", "grid")


def _short_sha(data: dict[str, Any]) -> str:
    text = json.dumps(data, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


class FileRef(BaseModel):
    model_config = _FROZEN
    path: str  # POSIX path relative to the run directory
    sha256: str


class Level(BaseModel):
    """Level of theory observed in the output (not the requested one); strings are lower-case."""

    model_config = _FROZEN
    program: str
    version: str
    method: str
    basis: str | None = None  # cartesian basis sets carry a "/cart" suffix
    dispersion: str | None = None
    charge: int
    multiplicity: int
    grid: str | None = None
    scf_tol: float | None = None
    electronic_temperature_K: float | None = None

    @field_validator(*_LEVEL_STRINGS, mode="before")
    @classmethod
    def _lower(cls, value: object) -> object:
        return value.lower() if isinstance(value, str) else value

    def full_key(self) -> str:
        return _short_sha(self.model_dump(mode="json"))


class Geometry(BaseModel):
    model_config = _FROZEN
    file: FileRef
    fingerprint: str  # chemistry.xyz.geometry_fingerprint; never computed in core
    symbols: tuple[str, ...]


Task = Literal["sp", "opt", "freq", "saddle"]


class Evidence(BaseModel):
    """Observed facts of one normally terminated external job. Its existence guarantees:
    (1) normal termination and SCF convergence  (2) geometry convergence for opt / saddle
    (3) 3N - n_external frequencies and the Hessian for freq  (4) preserved atom order
    (5) preserved input frame (echoed start geometry within 1e-4 A of the input)
    (6) observed Level matches the requested MethodSpec and the site's version pin
    (7) a gradient, when given, is the one at final.
    Backends return a Failure when any of these does not hold."""

    model_config = _FROZEN
    kind: Literal["calculation"] = "calculation"
    engine: str
    task: Task
    level: Level
    start: Geometry
    final: Geometry  # equals start for sp / freq
    energy_hartree: float  # electronic energy at final
    frequencies_cm1: tuple[float, ...] | None = None  # freq only; projected, imaginary < 0
    n_external: Literal[3, 5, 6] | None = None  # freq only; 5 if linear, 3 for an atom (no mode)
    imaginary_modes: tuple[tuple[float, ...], ...] = ()  # normalized cartesian, input frame
    hessian: FileRef | None = None  # freq only; canonical .npy (3N, 3N) in Eh/bohr^2
    gradient: tuple[float, ...] | None = None  # opt / saddle; 3N Eh/bohr at final, input frame
    s2: float | None = None  # observed <S^2> (open shell only)
    output: FileRef
    job_key: str

    @model_validator(mode="after")
    def _freq_modes(self) -> Evidence:
        """Guarantee (3), checked once here: 3N - n_external frequencies, one imaginary mode
        per negative frequency, and the Hessian they come from."""
        if self.task != "freq":
            return self
        if self.frequencies_cm1 is None or self.n_external is None or self.hessian is None:
            raise ValueError("a freq Evidence needs frequencies_cm1, n_external and hessian")
        expected = 3 * len(self.final.symbols) - self.n_external
        if len(self.frequencies_cm1) != expected:
            raise ValueError(f"frequency_count:{len(self.frequencies_cm1)}!={expected}")
        negative = sum(f < 0 for f in self.frequencies_cm1)
        if len(self.imaginary_modes) != negative:
            raise ValueError(f"imaginary_modes:{len(self.imaginary_modes)}!={negative}")
        return self


class PathProfile(BaseModel):
    model_config = _FROZEN
    kind: Literal["path"] = "path"
    engine: str
    level: Level
    images: FileRef  # multi-frame xyz
    energies_hartree: tuple[float, ...] = ()  # bead energies (nwchem_string); none from the NEB
    ts: Geometry | None = None  # TS optimized from the climbing image (pysis_neb)


class FailureKind(StrEnum):
    EXECUTABLE_MISSING = "executable_missing"
    INPUT_INVALID = "input_invalid"
    TIMEOUT = "timeout"
    NONZERO_EXIT = "nonzero_exit"
    SCF_NOT_CONVERGED = "scf_not_converged"
    GEOMETRY_MAXITER = "geometry_maxiter"
    INCOMPLETE_OUTPUT = "incomplete_output"
    METHOD_MISMATCH = "method_mismatch"
    GATE_REJECTED = "gate_rejected"


class Failure(BaseModel):
    model_config = _FROZEN
    kind: FailureKind
    reason: str
    final: Geometry | None = None  # last frame of a saddle search stopped at maxiter
    energy_hartree: float | None = None  # electronic energy at that frame
    job_key: str | None = None

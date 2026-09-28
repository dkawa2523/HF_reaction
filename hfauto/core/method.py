"""Requested methods, site execution settings and deadlines (design §5.2)."""

from __future__ import annotations

import math
import time
from dataclasses import dataclass
from typing import Literal

from pydantic import BaseModel, ConfigDict

from hfauto.core.evidence import Level

_DEFAULT_ELECTRONIC_TEMPERATURE_K = 300.0
_SCF_TOL_REL = 1e-9


class MethodSpec(BaseModel):
    """Requested level of theory (configs/methods/*.yaml), gas phase only; part of the JobStore
    key."""

    model_config = ConfigDict(frozen=True, extra="forbid")
    id: str
    kind: Literal["dft", "xtb", "wft"]
    functional: str | None = None  # dft
    wft_method: Literal["ccsd(t)"] | None = None  # wft (single points only)
    basis: str | None = None  # dft / wft
    dispersion: Literal["d3zero", "d3bj"] | None = None
    gfn: Literal[1, 2] | None = None  # xtb
    electronic_temperature_K: float | None = None  # xtb (None means 300 K)
    grid: Literal["coarse", "medium", "fine", "xfine"] | None = None
    scf_energy_tol: float | None = None

    def signature(self) -> dict[str, object]:
        """Every field except id (cache key)."""
        return self.model_dump(mode="json", exclude={"id"})


def _lower(value: str | None) -> str | None:
    return value.lower() if value is not None else None


def _diff(name: str, requested: object, observed: object) -> list[str]:
    if requested == observed:
        return []
    return [f"{name}: requested {requested!r}, observed {observed!r}"]


def _diff_ci(name: str, requested: str | None, observed: str | None) -> list[str]:
    return _diff(name, _lower(requested), _lower(observed))


def _scf_tol_diff(requested: float | None, observed: float | None) -> list[str]:
    if requested is None:
        return []
    if observed is not None and math.isclose(requested, observed, rel_tol=_SCF_TOL_REL):
        return []
    return [f"scf_tol: requested {requested!r}, observed {observed!r}"]


def _dft_mismatches(req: MethodSpec, obs: Level) -> list[str]:
    out = _diff_ci("functional", req.functional, obs.method)
    out += _diff_ci("basis", req.basis, obs.basis)
    out += _diff("dispersion", req.dispersion, obs.dispersion)
    if req.grid is not None:
        out += _diff("grid", req.grid, obs.grid)
    return out + _scf_tol_diff(req.scf_energy_tol, obs.scf_tol)


def _xtb_mismatches(req: MethodSpec, obs: Level) -> list[str]:
    out = _diff_ci("method", f"gfn{req.gfn}", obs.method)
    return out + _diff(
        "electronic_temperature_K",
        _electronic_temperature(req.electronic_temperature_K),
        _electronic_temperature(obs.electronic_temperature_K),
    )


def _electronic_temperature(value: float | None) -> float:
    return _DEFAULT_ELECTRONIC_TEMPERATURE_K if value is None else float(value)


def _wft_mismatches(req: MethodSpec, obs: Level) -> list[str]:
    out = _diff_ci("method", req.wft_method, obs.method)
    return out + _diff_ci("basis", req.basis, obs.basis)


_BY_KIND = {"dft": _dft_mismatches, "xtb": _xtb_mismatches, "wft": _wft_mismatches}


def level_mismatches(requested: MethodSpec, observed: Level, *, version_pin: str) -> list[str]:
    """Differences between the requested method and the observed level; empty means a match.

    Basis names compare case-insensitively, but a "/cart" suffix on one side is a mismatch.
    """
    out = _BY_KIND[requested.kind](requested, observed)
    return out + _diff_ci("version", version_pin, observed.version)


class ExecutionSpec(BaseModel):
    """Site execution settings, which do not change a result; not part of the JobStore key."""

    model_config = ConfigDict(frozen=True, extra="forbid")
    ranks: int = 1
    threads: int = 1
    memory_mb_per_rank: int = 1200
    timeout_s: float = 14_400
    env: dict[str, str] = {}


class EngineSite(BaseModel):
    """Per-engine site settings passed to engines.create."""

    model_config = ConfigDict(frozen=True, extra="forbid")
    version: str  # version pin: part of the JobStore key, checked against Level.version
    executables: dict[str, str] = {}
    execution: ExecutionSpec = ExecutionSpec()
    python: str | None = None  # worker interpreter (None means sys.executable)
    scratch_dir: str | None = None  # absolute path on ext4


class ThermoSettings(BaseModel):
    """Typed GoodVibes settings; temperatures come from ThermoConfig.temperatures_K."""

    model_config = ConfigDict(frozen=True, extra="forbid")
    qs: Literal["grimme", "truhlar"] = "grimme"
    cutoff_cm1: float = 100.0
    vib_scale: float = 1.0  # one factor for the frequencies and the ZPE
    sensitivity: bool = True  # qs x cutoff{50, 100, 150} band


@dataclass
class Deadline:
    end: float  # time.monotonic() reference

    @classmethod
    def after(cls, seconds: float) -> Deadline:
        return cls(end=time.monotonic() + seconds)

    def remaining(self) -> float:
        return max(0.0, self.end - time.monotonic())

    def expired(self) -> bool:
        return time.monotonic() >= self.end

"""MinimumDriver: opt → separate freq → mode-follow, and the minimum Registry (design §7.2).

Imports are limited to hfauto.core, hfauto.chemistry and hfauto.backends.protocols.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import numpy as np

from hfauto.chemistry.gates import Policy, imaginary_tier, is_minimum, spin_ok
from hfauto.chemistry.identity import BASIN_DE_HARTREE, assign, same_as_labelled
from hfauto.chemistry.modes import capped, classify_mode_follow, off_saddle
from hfauto.chemistry.topology import state_label
from hfauto.chemistry.vibrations import newton_step, stationarity_gap
from hfauto.chemistry.xyz import XYZ, Molecule, composition_key
from hfauto.core.evidence import Evidence, Failure, FailureKind, FileRef, Geometry
from hfauto.core.method import Deadline, MethodSpec
from hfauto.core.records import MinimumRecord, SpeciesRecord

if TYPE_CHECKING:
    from hfauto.backends.protocols import QMEngine

Status = Literal["minimum", "saddle", "known", "failed"]
LoadXYZ = Callable[[Geometry], XYZ]
Resolve = Callable[[FileRef], Path]
_GATES = Policy()


@dataclass(frozen=True)
class MinimumOutcome:
    status: Status
    opt: Evidence | None
    freq: Evidence | None
    history: tuple[str, ...]
    known_basin: str | None = None  # existing basin id when status == "known"
    failure: Failure | None = None
    # ± reached two minima distinct as labelled: a free TS candidate, with each side's outcome
    ts_candidate: tuple[MinimumOutcome, MinimumOutcome] | None = None
    notes: tuple[str, ...] = ()  # MinimumRecord.notes: *_imaginary_mode, spin_contaminated


def calc_id(ev: Evidence) -> str:
    """Artifact id of the calculation artifact that carries ``ev``."""
    return f"calc_{ev.job_key[:16]}"


def newton_push(hessian: np.ndarray, gradient: Sequence[float], xyz: XYZ, *,
                signed: bool = False) -> np.ndarray | None:
    """The Newton step (vibrations.newton_step, capped at the upper bound) from a point that is
    not stationary: ΔE_N (vibrations.stationarity_gap, its final gradient and freq Hessian) >
    BASIN_DE_HARTREE, the energy within which two minima are one basin. None at a stationary
    point. ``signed`` (a saddle): to the quadratic model's stationary point."""
    if stationarity_gap(hessian, gradient, xyz.symbols, xyz.coords) <= BASIN_DE_HARTREE:
        return None
    return capped(newton_step(hessian, gradient, xyz.coords, signed=signed))


@dataclass(frozen=True)
class _Point:
    """An optimized structure with its separate freq job; ``newton`` is the Newton step of a
    point that is not stationary."""

    opt: Evidence
    freq: Evidence
    xyz: XYZ
    tier: str
    notes: tuple[str, ...]
    newton: np.ndarray | None = None


@dataclass(frozen=True)
class _Ctx:
    template: Molecule  # charge and multiplicity of every structure
    method: MethodSpec
    qm: QMEngine
    load: LoadXYZ
    max_mode_follow: int
    gates: Policy
    deadline: Deadline | None
    resolve: Resolve | None

    def molecule(self, xyz: XYZ) -> Molecule:
        return Molecule(xyz=xyz, charge=self.template.charge,
                        multiplicity=self.template.multiplicity)

    def frequencies(self, opt: Evidence) -> _Point | Failure:
        """freq → is_minimum. A point below saddle order that is not stationary (``newton``) is
        soft whatever its frequencies, noted soft_imaginary_mode: noise_cm1 is numerical noise
        only at a stationary point."""
        xyz = self.load(opt.final)
        freq = self.qm.frequencies(self.molecule(xyz), self.method, scf_guess=opt,
                                   deadline=self.deadline)
        if isinstance(freq, Failure):
            return freq
        gate = is_minimum(freq, opt=opt, policy=self.gates)
        hard = [reason for reason in gate.reasons if reason != "imaginary_mode"]
        if hard:
            return Failure(kind=FailureKind.GATE_REJECTED, reason=",".join(hard),
                           job_key=freq.job_key)
        tier, notes = imaginary_tier(freq.frequencies_cm1 or (), self.gates), gate.notes
        newton = None if tier == "saddle" else self.newton(opt, freq, xyz)
        if newton is not None:
            tier, notes = "soft", ("soft_imaginary_mode",)
        return _Point(opt, freq, xyz, tier, notes + spin_ok(freq, self.gates).reasons, newton)

    def newton(self, opt: Evidence, freq: Evidence, xyz: XYZ) -> np.ndarray | None:
        """newton_push −H₊⁺g from the opt's final structure; not judged without a gradient (an
        opt stored before Evidence carried one)."""
        if opt.gradient is None or freq.hessian is None or self.resolve is None:
            return None
        return newton_push(np.load(self.resolve(freq.hessian)), opt.gradient, xyz)

    def relax(self, coords: np.ndarray, source: _Point) -> _Point | None:
        """opt → freq from coordinates displaced from ``source``, whose freq Hessian starts the
        opt (the engine writes it as its positive-definite model unless it is a higher-order
        saddle's; trust 0.3); None when either job fails."""
        xyz = XYZ(symbols=list(self.template.xyz.symbols), coords=coords)
        opt = self.qm.optimize(self.molecule(xyz), self.method, init_hessian=source.freq,
                               deadline=self.deadline)
        if isinstance(opt, Failure):
            return None
        point = self.frequencies(opt)
        return point if isinstance(point, _Point) else None


def _labels(source: _Point, sides: list[_Point | None]) -> list[str | None]:
    """Identity labels (source, plus, minus): equal labels mean one structure as labelled, so
    the two structures of a degenerate rearrangement (the NH3 inversion) get two labels."""
    seen: list[tuple[_Point, str]] = [(source, "source")]
    labels: list[str | None] = ["source"]
    for name, side in zip(("plus", "minus"), sides, strict=True):
        if side is None:
            labels.append(None)
            continue
        label = next((known for p, known in seen if same_as_labelled(
            side.xyz.coords, p.xyz.coords, side.opt.energy_hartree, p.opt.energy_hartree)), name)
        seen.append((side, label))
        labels.append(label)
    return labels


def _follow(ctx: _Ctx, point: _Point, history: list[str]
            ) -> tuple[_Point, tuple[MinimumOutcome, MinimumOutcome] | None]:
    """Up to max_mode_follow cycles from a saddle, along its modes below -saddle_cm1
    (modes.off_saddle): ± along the one of a first-order saddle, its QRC step (a TS candidate
    when both sides are minima, each side settled as its own outcome); one side along all of a
    higher-order one (± would mostly stop at first-order saddles)."""
    for cycle in range(1, ctx.max_mode_follow + 1):
        if point.tier != "saddle" or not point.freq.imaginary_modes:
            break
        step = off_saddle(point.freq, point.xyz.symbols, below_cm1=ctx.gates.saddle_cm1)
        order = sum(nu < -ctx.gates.saddle_cm1 for nu in point.freq.frequencies_cm1 or ())
        sides = [ctx.relax(point.xyz.coords + step, point),
                 ctx.relax(point.xyz.coords - step, point) if order == 1 else None]
        labels = _labels(point, sides)
        verdict = classify_mode_follow(*labels)
        history.append(f"follow{cycle}:{verdict}")
        if verdict == "same_as_source":
            break
        reached = [s for s, label in zip(sides, labels[1:], strict=True)
                   if s is not None and label != "source"]
        if verdict == "ts_candidate" and all(s.tier != "saddle" for s in reached):
            a, b = (_outcome(ctx, s, [f"opt:follow{cycle}", f"freq:{s.tier}"]) for s in reached)
            return point, (a, b)
        point = min(reached, key=lambda s: (s.tier == "saddle", s.opt.energy_hartree))
    return point, None


def _soften(ctx: _Ctx, point: _Point, history: list[str]) -> _Point:
    """One push from a soft point, then relax: its Newton step when it is not stationary
    (whether or not it has imaginary modes), else along its soft imaginary modes
    (modes.off_saddle). The relaxed point when it is no longer soft (soft:resolved); else the
    point itself (soft:persisted, noted soft_imaginary_mode)."""
    step = point.newton if point.newton is not None else off_saddle(
        point.freq, point.xyz.symbols, below_cm1=ctx.gates.noise_cm1)
    side = ctx.relax(point.xyz.coords + step, point)
    if side is not None and side.tier in ("none", "noise"):
        history.append("soft:resolved")
        return side
    history.append("soft:persisted")
    return point


def _outcome(ctx: _Ctx, point: _Point, history: list[str],
             pair: tuple[MinimumOutcome, MinimumOutcome] | None = None) -> MinimumOutcome:
    """A saddle as it is; else a minimum, after one push when it is soft."""
    status: Status = "saddle" if point.tier == "saddle" else "minimum"
    if point.tier == "soft":
        point = _soften(ctx, point, history)
    return MinimumOutcome(status, point.opt, point.freq, tuple(history), ts_candidate=pair,
                          notes=point.notes)


def _failed(history: list[str], failure: Failure, opt: Evidence | None = None) -> MinimumOutcome:
    history.append(f"failed:{failure.kind.value}")
    return MinimumOutcome("failed", opt, None, tuple(history), failure=failure)


def relax_to_minimum(
    mol: Molecule,
    method: MethodSpec,
    qm: QMEngine,
    *,
    known: Registry | None = None,
    init_hessian: Evidence | None = None,
    opt: Evidence | None = None,
    max_mode_follow: int = 2,
    gates: Policy = _GATES,
    deadline: Deadline | None = None,
    load_xyz: LoadXYZ | None = None,
    resolve: Resolve | None = None,
) -> MinimumOutcome:
    """opt (with init_hessian) → known check → freq → mode-follow (saddle) / one push (soft).

    A converged ``opt`` is taken as it is (``mol`` then gives only charge and multiplicity).
    ``load_xyz`` (default ``known.load_xyz``) reads the optimized geometry for the freq job;
    ``resolve`` opens a freq's Hessian, which judges stationarity (not judged without it).
    """
    load = load_xyz or (known.load_xyz if known is not None else None)
    if load is None:
        raise ValueError("relax_to_minimum needs load_xyz or a known Registry")
    ctx = _Ctx(mol, method, qm, load, max_mode_follow, gates, deadline, resolve)
    history = ["opt" if opt is None else "opt:reused"]
    done = opt if opt is not None else qm.optimize(
        mol, method, init_hessian=init_hessian, deadline=deadline)
    if isinstance(done, Failure):
        return _failed(history, done)
    if known is not None:
        basin = known.find(done)
        if basin is not None:
            history.append(f"known:{basin}")
            return MinimumOutcome("known", done, None, tuple(history), known_basin=basin)
    point = ctx.frequencies(done)
    if isinstance(point, Failure):
        return _failed(history, point, done)
    history.append(f"freq:{point.tier}")
    point, pair = _follow(ctx, point, history)
    return _outcome(ctx, point, history, pair)


@dataclass
class _Entry:
    record: MinimumRecord
    geometry: Geometry  # the representative's optimized structure
    xyz: XYZ  # the same, loaded


class Registry:
    """The store of the minima of a tier: per composition × level_key, one identity criterion,
    identity.assign (mirror images are one basin). ``minima`` pairs each existing record with
    the optimized structure of its representative (a MinimumRecord carries no geometry); a
    record grows members as species join."""

    def __init__(self, minima: Iterable[tuple[MinimumRecord, Geometry]],
                 load_xyz: LoadXYZ) -> None:
        self.load_xyz = load_xyz
        self._basins = {r.basin_id: _Entry(r, g, load_xyz(g)) for r, g in minima}

    @property
    def minima(self) -> Mapping[str, tuple[MinimumRecord, Geometry]]:
        """Minimum id → (current record, its representative's optimized structure)."""
        return {e.record.minimum_id: (e.record, e.geometry) for e in self._basins.values()}

    def basin(self, basin_id: str) -> MinimumRecord:
        """The current record of ``basin_id``."""
        return self._basins[basin_id].record

    def find(self, opt: Evidence, coords: np.ndarray | None = None) -> str | None:
        """Basin id that identity.assign uniquely matches with the optimized structure of
        ``opt`` (or ``coords`` in its atom order: an image of it), among the basins of its
        element list, composition (charge and multiplicity of its Level) and level_key; else
        None."""
        symbols, level = list(opt.final.symbols), opt.level
        x = self.load_xyz(opt.final).coords if coords is None else coords
        key = (composition_key(symbols, level.charge, level.multiplicity), level.full_key())
        candidates = {basin: (entry.xyz.coords, entry.record.energy_hartree)
                      for basin, entry in self._basins.items()
                      if list(entry.xyz.symbols) == symbols
                      and (entry.record.composition_id, entry.record.level_key) == key}
        return assign(symbols, x, opt.energy_hartree, candidates)

    def add(self, outcome: MinimumOutcome, species: SpeciesRecord, *,
            tier: Literal["screen", "dft"]) -> MinimumRecord:
        """Join the known basin, or the basin that ``find`` matches; else register a new
        basin."""
        if outcome.status == "known" and outcome.known_basin is not None:
            return self.join(outcome.known_basin, species.species_id)
        opt, freq = outcome.opt, outcome.freq
        if outcome.status != "minimum" or opt is None or freq is None:
            raise ValueError(f"cannot register a {outcome.status!r} outcome")
        basin = self.find(opt)
        if basin is not None:
            return self.join(basin, species.species_id)
        xyz = self.load_xyz(opt.final)
        record = _new_record(opt, freq, outcome.notes, species, tier, xyz)
        self._basins[record.basin_id] = _Entry(record, opt.final, xyz)
        return record

    def join(self, basin_id: str, species_id: str, *notes: str) -> MinimumRecord:
        """Add ``species_id`` to the members of ``basin_id`` and ``notes`` to its record."""
        entry = self._basins[basin_id]
        members = tuple(dict.fromkeys((*entry.record.members, species_id)))
        notes = tuple(dict.fromkeys((*entry.record.notes, *notes)))
        entry.record = entry.record.model_copy(update={"members": members, "notes": notes})
        return entry.record


def _new_record(opt: Evidence, freq: Evidence, notes: tuple[str, ...], species: SpeciesRecord,
                tier: Literal["screen", "dft"], xyz: XYZ) -> MinimumRecord:
    level_key = opt.level.full_key()
    stem = f"{species.species_id}_{level_key[:8]}"
    return MinimumRecord(
        minimum_id=f"min_{stem}", basin_id=f"basin_{stem}",
        composition_id=species.composition_id, species_id=species.species_id, tier=tier,
        level_key=level_key, opt_calc=calc_id(opt), freq_calc=calc_id(freq),
        energy_hartree=opt.energy_hartree, state_label=state_label(xyz.symbols, xyz.coords),
        members=(species.species_id,), notes=notes,
    )

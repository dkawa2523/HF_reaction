"""MinimumDriver: opt → separate freq → mode-follow, and the minimum Registry (design §7.2).

Imports are limited to hfauto.core, hfauto.chemistry and hfauto.backends.protocols.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np

from hfauto.chemistry.gates import Policy, imaginary_tier, is_minimum, spin_ok
from hfauto.chemistry.identity import assign, is_chiral, same_as_labelled
from hfauto.chemistry.modes import BOUNDS_A, amplitude, classify_mode_follow
from hfauto.chemistry.topology import state_label
from hfauto.chemistry.xyz import XYZ, Molecule, composition_key
from hfauto.core.evidence import Evidence, Failure, FailureKind, Geometry
from hfauto.core.method import Deadline, MethodSpec
from hfauto.core.records import MinimumRecord, SpeciesRecord

if TYPE_CHECKING:
    from hfauto.backends.protocols import QMEngine

Status = Literal["minimum", "soft_minimum", "saddle", "known", "failed"]
LoadXYZ = Callable[[Geometry], XYZ]
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


@dataclass(frozen=True)
class _Point:
    """An optimized structure with its separate freq job."""

    opt: Evidence
    freq: Evidence
    xyz: XYZ
    tier: str
    notes: tuple[str, ...]

    def step(self, below_cm1: float) -> tuple[np.ndarray, int]:
        """Sum of the imaginary modes below -below_cm1, each by its energy-target amplitude
        (modes.amplitude), with the largest atomic displacement capped at the upper bound; and
        the number of those modes."""
        symbols, freq = self.xyz.symbols, self.freq
        pairs = [(nu, np.asarray(mode, dtype=float).reshape(-1, 3)) for nu, mode in
                 zip(freq.frequencies_cm1 or (), freq.imaginary_modes, strict=False)
                 if nu < -below_cm1]
        step = np.sum([amplitude(nu, mode, symbols) / np.linalg.norm(mode, axis=1).max() * mode
                       for nu, mode in pairs], axis=0)
        largest = float(np.linalg.norm(step, axis=1).max())
        return step * min(1.0, BOUNDS_A[1] / largest), len(pairs)


@dataclass(frozen=True)
class _Ctx:
    template: Molecule  # charge and multiplicity of every structure
    method: MethodSpec
    qm: QMEngine
    load: LoadXYZ
    max_mode_follow: int
    gates: Policy
    deadline: Deadline | None

    def molecule(self, xyz: XYZ) -> Molecule:
        return Molecule(xyz=xyz, charge=self.template.charge,
                        multiplicity=self.template.multiplicity)

    def frequencies(self, opt: Evidence) -> _Point | Failure:
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
        tier = imaginary_tier(freq.frequencies_cm1 or (), self.gates)
        notes = gate.notes + spin_ok(freq, self.gates).reasons
        return _Point(opt, freq, xyz, tier, notes)

    def relax(self, coords: np.ndarray, source: _Point) -> _Point | None:
        """opt → freq from coordinates displaced from ``source``, whose freq Hessian starts the
        opt (the engine writes a first-order saddle's as its positive-definite model,
        trust 0.3); None when either job fails."""
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
    """Up to max_mode_follow cycles from a saddle: ± along the imaginary mode of a first-order
    saddle (a TS candidate when both sides are minima, each side settled as its own outcome);
    one side along all the modes below -saddle_cm1 of a higher-order one (± would mostly stop
    at first-order saddles)."""
    for cycle in range(1, ctx.max_mode_follow + 1):
        if point.tier != "saddle" or not point.freq.imaginary_modes:
            break
        step, order = point.step(ctx.gates.saddle_cm1)
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


def _soften(ctx: _Ctx, point: _Point, history: list[str]) -> tuple[_Point, Status]:
    """One displacement along the soft imaginary modes; soft_minimum when one persists."""
    if point.freq.imaginary_modes:
        side = ctx.relax(point.xyz.coords + point.step(ctx.gates.noise_cm1)[0], point)
        if side is not None and side.tier in ("none", "noise"):
            history.append("soft:resolved")
            return side, "minimum"
    history.append("soft:persisted")
    return point, "soft_minimum"


def _outcome(ctx: _Ctx, point: _Point, history: list[str],
             pair: tuple[MinimumOutcome, MinimumOutcome] | None = None) -> MinimumOutcome:
    """A saddle as it is; else a minimum, after one push when a soft mode remains."""
    status: Status = "saddle" if point.tier == "saddle" else "minimum"
    if point.tier == "soft":
        point, status = _soften(ctx, point, history)
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
) -> MinimumOutcome:
    """opt (with init_hessian) → known check → freq → mode-follow (saddle) / one push (soft).

    A converged ``opt`` is taken as it is (``mol`` then gives only charge and multiplicity).
    ``load_xyz`` (default ``known.load_xyz``) reads the optimized geometry for the freq job.
    """
    load = load_xyz or (known.load_xyz if known is not None else None)
    if load is None:
        raise ValueError("relax_to_minimum needs load_xyz or a known Registry")
    ctx = _Ctx(mol, method, qm, load, max_mode_follow, gates, deadline)
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
    xyz: XYZ


class Registry:
    """Minima per composition × level_key; one identity criterion, identity.assign (mirror
    images are one basin).

    ``minima`` pairs each existing record with the geometry of its representative
    (MinimumRecord itself carries no geometry).
    """

    def __init__(self, minima: Iterable[tuple[MinimumRecord, Geometry]],
                 load_xyz: LoadXYZ) -> None:
        self.load_xyz = load_xyz
        self._basins: dict[str, _Entry] = {
            record.basin_id: _Entry(record, load_xyz(geometry)) for record, geometry in minima
        }

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
        if outcome.status not in ("minimum", "soft_minimum") or opt is None or freq is None:
            raise ValueError(f"cannot register a {outcome.status!r} outcome")
        basin = self.find(opt)
        if basin is not None:
            return self.join(basin, species.species_id)
        xyz = self.load_xyz(opt.final)
        record = _new_record(opt, freq, outcome.notes, species, tier, xyz)
        self._basins[record.basin_id] = _Entry(record, xyz)
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
        members=(species.species_id,), chiral=is_chiral(xyz.symbols, xyz.coords), notes=notes,
    )

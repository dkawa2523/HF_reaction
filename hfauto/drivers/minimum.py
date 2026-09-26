"""MinimumDriver: opt → separate freq → mode-follow, and the minimum Registry (design §7.2).

Imports are limited to hfauto.core, hfauto.chemistry and hfauto.backends.protocols.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal

import numpy as np

from hfauto.chemistry.gates import Policy, imaginary_tier, is_minimum, spin_ok
from hfauto.chemistry.identity import assign, same_minimum
from hfauto.chemistry.modes import classify_mode_follow, displace
from hfauto.chemistry.topology import state_label
from hfauto.chemistry.xyz import XYZ, Molecule
from hfauto.core.evidence import Evidence, Failure, FailureKind, Geometry
from hfauto.core.method import Deadline, MethodSpec
from hfauto.core.records import MinimumRecord, SpeciesRecord

if TYPE_CHECKING:
    from hfauto.backends.protocols import QMEngine

Status = Literal["minimum", "soft_minimum", "saddle", "known", "failed"]
LoadXYZ = Callable[[Geometry], XYZ]


@dataclass(frozen=True)
class MinimumPolicy:
    max_mode_follow: int = 2
    amplitude_A: float = 0.1
    gates: Policy = field(default_factory=Policy)


@dataclass(frozen=True)
class MinimumOutcome:
    status: Status
    opt: Evidence | None
    freq: Evidence | None
    history: tuple[str, ...]
    known_basin: str | None = None  # existing basin id when status == "known"
    failure: Failure | None = None
    # ± displacements reached two distinct minima (their freq Evidence): a free TS candidate
    ts_candidate: tuple[Evidence, Evidence] | None = None
    notes: tuple[str, ...] = ()  # MinimumRecord.notes: *_imaginary_mode, spin_contaminated


_DEFAULT = MinimumPolicy()


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


@dataclass(frozen=True)
class _Ctx:
    template: Molecule  # charge and multiplicity of every structure
    method: MethodSpec
    qm: QMEngine
    load: LoadXYZ
    policy: MinimumPolicy
    deadline: Deadline | None

    def molecule(self, xyz: XYZ) -> Molecule:
        return Molecule(xyz=xyz, charge=self.template.charge,
                        multiplicity=self.template.multiplicity)

    def frequencies(self, opt: Evidence) -> _Point | Failure:
        xyz = self.load(opt.final)
        freq = self.qm.frequencies(self.molecule(xyz), self.method, deadline=self.deadline)
        if isinstance(freq, Failure):
            return freq
        gate = is_minimum(freq, opt=opt, policy=self.policy.gates)
        hard = [reason for reason in gate.reasons if reason != "imaginary_mode"]
        if hard:
            return Failure(kind=FailureKind.GATE_REJECTED, reason=",".join(hard),
                           job_key=freq.job_key)
        tier = imaginary_tier(freq.frequencies_cm1 or (), self.policy.gates)
        notes = gate.notes + spin_ok(freq, self.policy.gates).reasons
        return _Point(opt, freq, xyz, tier, notes)

    def relax(self, coords: np.ndarray) -> _Point | None:
        """opt → freq from displaced coordinates; None when either job fails."""
        xyz = XYZ(symbols=list(self.template.xyz.symbols), coords=coords)
        opt = self.qm.optimize(self.molecule(xyz), self.method, deadline=self.deadline)
        if isinstance(opt, Failure):
            return None
        point = self.frequencies(opt)
        return point if isinstance(point, _Point) else None


def _labels(source: _Point, sides: list[_Point | None]) -> list[str | None]:
    """Identity labels (source, plus, minus): equal labels mean the same structure."""

    def same(a: _Point, b: _Point) -> bool:
        return same_minimum(a.xyz.symbols, a.xyz.coords, b.xyz.coords,
                            a.opt.energy_hartree, b.opt.energy_hartree)

    seen: list[tuple[_Point, str]] = [(source, "source")]
    labels: list[str | None] = ["source"]
    for name, side in zip(("plus", "minus"), sides, strict=True):
        if side is None:
            labels.append(None)
            continue
        label = next((known for point, known in seen if same(side, point)), name)
        seen.append((side, label))
        labels.append(label)
    return labels


def _follow(ctx: _Ctx, point: _Point, history: list[str]
            ) -> tuple[_Point, tuple[Evidence, Evidence] | None]:
    """Up to max_mode_follow cycles of ± displacement along the lowest imaginary mode."""
    for cycle in range(1, ctx.policy.max_mode_follow + 1):
        if point.tier != "saddle" or not point.freq.imaginary_modes:
            break
        plus, minus = displace(point.xyz.coords, np.asarray(point.freq.imaginary_modes[0]),
                               ctx.policy.amplitude_A)
        sides = [ctx.relax(plus), ctx.relax(minus)]
        labels = _labels(point, sides)
        verdict = classify_mode_follow(*labels)
        history.append(f"follow{cycle}:{verdict}")
        if verdict == "same_as_source":
            break
        reached = [s for s, label in zip(sides, labels[1:], strict=True)
                   if s is not None and label != "source"]
        if verdict == "ts_candidate" and all(s.tier != "saddle" for s in reached):
            return point, (reached[0].freq, reached[1].freq)
        point = min(reached, key=lambda s: (s.tier == "saddle", s.opt.energy_hartree))
    return point, None


def _soften(ctx: _Ctx, point: _Point, history: list[str]) -> tuple[_Point, Status]:
    """One displacement along a soft imaginary mode; soft_minimum when it persists."""
    if point.freq.imaginary_modes:
        plus, _ = displace(point.xyz.coords, np.asarray(point.freq.imaginary_modes[0]),
                           ctx.policy.amplitude_A)
        side = ctx.relax(plus)
        if side is not None and side.tier in ("none", "noise"):
            history.append("soft:resolved")
            return side, "minimum"
    history.append("soft:persisted")
    return point, "soft_minimum"


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
    policy: MinimumPolicy = _DEFAULT,
    deadline: Deadline | None = None,
    load_xyz: LoadXYZ | None = None,
) -> MinimumOutcome:
    """opt (with init_hessian) → known check → freq → mode-follow (saddle) / one push (soft).

    ``load_xyz`` defaults to ``known.load_xyz``; one of them is required because the
    driver re-reads the optimized geometry before the separate freq job.
    """
    load = load_xyz or (known.load_xyz if known is not None else None)
    if load is None:
        raise ValueError("relax_to_minimum needs load_xyz or a known Registry")
    ctx = _Ctx(mol, method, qm, load, policy, deadline)
    history = ["opt"]
    opt = qm.optimize(mol, method, init_hessian=init_hessian, deadline=deadline)
    if isinstance(opt, Failure):
        return _failed(history, opt)
    if known is not None:
        final = load(opt.final)
        basin = known.find(final.symbols, final.coords, opt.energy_hartree)
        if basin is not None:
            history.append(f"known:{basin}")
            return MinimumOutcome("known", opt, None, tuple(history), known_basin=basin)
    point = ctx.frequencies(opt)
    if isinstance(point, Failure):
        return _failed(history, point, opt)
    history.append(f"freq:{point.tier}")
    point, pair = _follow(ctx, point, history)
    status: Status = "saddle" if point.tier == "saddle" else "minimum"
    if point.tier == "soft":
        point, status = _soften(ctx, point, history)
    return MinimumOutcome(status, point.opt, point.freq, tuple(history),
                          ts_candidate=pair, notes=point.notes)


@dataclass
class _Entry:
    record: MinimumRecord
    xyz: XYZ


class Registry:
    """Minima per composition × level_key; one identity criterion, identity.assign.

    ``minima`` pairs each existing record with the geometry of its representative
    (MinimumRecord itself carries no geometry).
    """

    def __init__(self, minima: Iterable[tuple[MinimumRecord, Geometry]],
                 load_xyz: LoadXYZ) -> None:
        self.load_xyz = load_xyz
        self._basins: dict[str, _Entry] = {
            record.basin_id: _Entry(record, load_xyz(geometry)) for record, geometry in minima
        }

    def find(self, symbols: Sequence[str], coords: np.ndarray, energy: float, *,
             composition_id: str | None = None, level_key: str | None = None) -> str | None:
        """Basin id uniquely matching the structure (identity.assign) among the basins with
        the same element list (and the given composition_id / level_key), else None."""
        candidates = {basin: (entry.xyz.coords, entry.record.energy_hartree)
                      for basin, entry in self._basins.items()
                      if list(entry.xyz.symbols) == list(symbols)
                      and composition_id in (None, entry.record.composition_id)
                      and level_key in (None, entry.record.level_key)}
        return assign(symbols, coords, energy, candidates)

    def members(self, basin_id: str) -> tuple[str, ...]:
        return self._basins[basin_id].record.members

    def add(self, outcome: MinimumOutcome, species: SpeciesRecord, *,
            tier: Literal["screen", "dft"]) -> MinimumRecord:
        """Join the known basin, or the basin of the same composition_id and level_key that
        identity.assign uniquely matches; else register a new basin."""
        if outcome.status == "known" and outcome.known_basin is not None:
            return self._join(outcome.known_basin, species.species_id)
        opt, freq = outcome.opt, outcome.freq
        if outcome.status not in ("minimum", "soft_minimum") or opt is None or freq is None:
            raise ValueError(f"cannot register a {outcome.status!r} outcome")
        xyz = self.load_xyz(opt.final)
        basin = self.find(xyz.symbols, xyz.coords, opt.energy_hartree,
                          composition_id=species.composition_id, level_key=opt.level.full_key())
        if basin is not None:
            return self._join(basin, species.species_id)
        record = _new_record(opt, freq, outcome.notes, species, tier, xyz)
        self._basins[record.basin_id] = _Entry(record, xyz)
        return record

    def _join(self, basin_id: str, species_id: str) -> MinimumRecord:
        entry = self._basins[basin_id]
        if species_id not in entry.record.members:
            members = (*entry.record.members, species_id)
            entry.record = entry.record.model_copy(update={"members": members})
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

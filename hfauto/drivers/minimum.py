"""MinimumDriver: opt → separate freq → certification → mode-follow, the one way down from a
saddle (``Relaxer.descend``) and the minimum Registry (design §6.1).

Imports are limited to hfauto.core, hfauto.chemistry and hfauto.backends.protocols.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import numpy as np

from hfauto.chemistry.gates import Policy, imaginary_tier, is_minimum, spin_ok
from hfauto.chemistry.identity import (
    BASIN_A,
    BASIN_DE_HARTREE,
    assign,
    carry,
    is_image,
    same_as_labelled,
)
from hfauto.chemistry.modes import classify_mode_follow, off_saddle
from hfauto.chemistry.topology import state_label
from hfauto.chemistry.vibrations import trust_region_step
from hfauto.chemistry.xyz import XYZ, Molecule, composition_key
from hfauto.core.evidence import Evidence, Failure, FailureKind, FileRef, Geometry
from hfauto.core.method import MethodSpec
from hfauto.core.records import MinimumRecord, SpeciesRecord

if TYPE_CHECKING:
    from hfauto.backends.protocols import QMEngine

Status = Literal["minimum", "saddle", "known", "failed"]
LoadXYZ = Callable[[Geometry], XYZ]
Resolve = Callable[[FileRef], Path]
Find = Callable[[Evidence, np.ndarray], str | None]  # (opt, its coordinates) -> a known basin
Map = Callable[[Callable[[Any], Any], Sequence[Any]], list[Any]]
_GATES = Policy()


def _serial(fn: Callable[[Any], Any], items: Sequence[Any]) -> list[Any]:
    return [fn(x) for x in items]


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
    # MinimumRecord.notes: *_imaginary_mode, spin_contaminated, not_stationary
    notes: tuple[str, ...] = ()


def calc_id(ev: Evidence) -> str:
    """Artifact id of the calculation artifact that carries ``ev``."""
    return f"calc_{ev.job_key[:16]}"


def not_stationary(freq: Evidence, gradient: Sequence[float] | None, coords: np.ndarray,
                   resolve: Resolve | None) -> np.ndarray | None:
    """The trust-region step (N, 3; Å) from a point whose certification fails (design §6.1):
    ΔE_TR of vibrations.trust_region_step, from its final gradient and the Hessian of ``freq``
    at ``coords``, within one basin's radius BASIN_A·√N (the identity RMSD), exceeds
    BASIN_DE_HARTREE, the energy of one basin. Minima and TSs alike, |H| over all modes: without
    the reaction mode a TS off its ridge would pass. None at a stationary point; not judged
    without a gradient, a Hessian or ``resolve`` (an older stored result)."""
    if gradient is None or freq.hessian is None or resolve is None:
        return None
    x = np.reshape(coords, (-1, 3))
    gap, step = trust_region_step(np.load(resolve(freq.hessian)), gradient, x,
                                  BASIN_A * math.sqrt(len(x)))
    return step if gap > BASIN_DE_HARTREE else None


@dataclass(frozen=True, eq=False)
class Settled:
    """An optimized structure at ``x`` (an image side: the other side's optimum carried onto
    its own start), settled: in a known ``basin`` (no freq job), or with its freq, imaginary
    tier and notes."""

    opt: Evidence
    x: np.ndarray
    history: tuple[str, ...]
    basin: str | None = None
    freq: Evidence | None = None
    tier: str = "none"
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class Relaxer:
    """What settling and descending need: ``template`` gives charge and multiplicity, ``known``
    finds a known basin (Registry.find; a case adds its own fallback), ``map`` runs independent
    optimizations at once."""

    template: Molecule
    method: MethodSpec
    qm: QMEngine
    load: LoadXYZ
    gates: Policy = _GATES
    resolve: Resolve | None = None
    known: Find | None = None
    map: Map = field(default=_serial)
    max_mode_follow: int = 2

    def molecule(self, coords: np.ndarray) -> Molecule:
        xyz = XYZ(symbols=list(self.template.xyz.symbols), coords=np.reshape(coords, (-1, 3)))
        return Molecule(xyz=xyz, charge=self.template.charge,
                        multiplicity=self.template.multiplicity)

    def settle(self, opt: Evidence, history: Sequence[str] = (), relaxed: bool = False
               ) -> Settled | Failure:
        """A known basin needs no freq; else freq → is_minimum. A point below saddle order is
        certified (``not_stationary``); one that fails relaxes once from x + s with its freq as
        the positive model, settled again, and is noted not_stationary when it still fails
        (thermo blocks it). noise_cm1 is the numerical noise of a stationary point only."""
        x = np.asarray(self.load(opt.final).coords, dtype=float)
        basin = None if self.known is None else self.known(opt, x)
        if basin is not None:
            return Settled(opt, x, (*history, f"known:{basin}"), basin=basin)
        freq = self.qm.frequencies(self.molecule(x), self.method, scf_guess=opt)
        if isinstance(freq, Failure):
            return freq
        gate = is_minimum(freq, opt=opt, policy=self.gates)
        if hard := [reason for reason in gate.reasons if reason != "imaginary_mode"]:
            return Failure(kind=FailureKind.GATE_REJECTED, reason=",".join(hard),
                           job_key=freq.job_key)
        tier = imaginary_tier(freq.frequencies_cm1 or (), self.gates)
        point = Settled(opt, x, (*history, f"freq:{tier}"), freq=freq, tier=tier,
                        notes=(*gate.notes, *spin_ok(freq, self.gates).reasons))
        step = None if tier == "saddle" else not_stationary(freq, opt.gradient, x, self.resolve)
        if step is None:
            return point
        uncertified = replace(point, notes=(*point.notes, "not_stationary"))
        if relaxed:
            return uncertified
        moved = self.qm.optimize(self.molecule(x + step), self.method, init_hessian=freq)
        again = moved if isinstance(moved, Failure) else self.settle(
            moved, (*point.history, "tr_relax"), relaxed=True)
        return uncertified if isinstance(again, Failure) else again

    def descend(self, x: np.ndarray, freq: Evidence, step: np.ndarray, *, both: bool,
                name: str, carried: bool = False, past_saddles: bool = False
                ) -> list[Settled | Failure]:
        """The one way down from a saddle at ``x`` with its ``freq``: ± ``step`` (one side
        unless ``both``), each optimized from that freq (its positive model from a first-order
        saddle, as is from a higher-order one, whose side must stay free to leave its other
        saddle directions), at once. ``carried``: a minus start that is an exact image of the
        plus one (identity.is_image) is not optimized; the plus optimum carried onto it stands
        for it. Each side is then settled, a known basin with no freq and a side in the basin of
        the side before it (an image) with that side's. ``past_saddles``: a side that stops on a
        saddle goes on down once, away from ``x`` along its imaginary modes."""
        starts = [x + step, x - step][:2 if both else 1]
        image = carried and both and is_image(self.template.xyz.symbols, starts[1], starts[0])
        order = sum(nu < -self.gates.saddle_cm1 for nu in freq.frequencies_cm1 or ())
        runs = self.map(lambda y: self.qm.optimize(
            self.molecule(y), self.method, init_hessian=freq,
            hessian_model="positive" if order <= 1 else "as_is"), starts[:1] if image else starts)
        sides: list[Settled | Failure] = []
        for opt in runs:
            side = opt if isinstance(opt, Failure) else self._side(opt, sides, f"opt:{name}")
            if past_saddles and isinstance(side, Settled) and side.tier == "saddle":
                side = self._past(side, x, name)
            sides.append(side)
        if image:
            plus = sides[0]
            sides.append(plus if isinstance(plus, Failure) else replace(plus, x=carry(
                self.template.xyz.symbols, starts[1], starts[0], plus.x)[1]))
        return sides

    def _side(self, opt: Evidence, earlier: Sequence[Settled | Failure], step: str
              ) -> Settled | Failure:
        """A side in the basin of an earlier one shares its settlement; else its own."""
        x, symbols = np.asarray(self.load(opt.final).coords), self.template.xyz.symbols
        settled = [s for s in earlier if isinstance(s, Settled)]
        found = assign(symbols, x, opt.energy_hartree,
                       {str(k): (s.x, s.opt.energy_hartree) for k, s in enumerate(settled)})
        if found is None:
            return self.settle(opt, (step,))
        return replace(settled[int(found)], opt=opt, x=x, history=(step, "same_basin"))

    def _past(self, side: Settled, x: np.ndarray, name: str) -> Settled | Failure:
        """``side`` (a saddle) settled once more, one side down along its modes below
        -saddle_cm1, the sign leading away from ``x``."""
        assert side.freq is not None
        step = off_saddle(side.freq, self.template.xyz.symbols, below_cm1=self.gates.saddle_cm1)
        sign = 1.0 if float(np.vdot(step, side.x - np.reshape(x, side.x.shape))) >= 0 else -1.0
        (down,) = self.descend(side.x, side.freq, sign * step, both=False, name=f"{name}_past")
        return down


def _labels(source: Settled, sides: list[Settled | Failure]) -> list[str | None]:
    """Identity labels (source, plus, minus): equal labels mean one structure as labelled, so
    the two structures of a degenerate rearrangement (the NH3 inversion) get two labels."""
    seen: list[tuple[Settled, str]] = [(source, "source")]
    labels: list[str | None] = ["source"]
    for name, side in zip(("plus", "minus"), sides, strict=False):
        if isinstance(side, Failure):
            labels.append(None)
            continue
        label = next((known for p, known in seen if same_as_labelled(
            side.x, p.x, side.opt.energy_hartree, p.opt.energy_hartree)), name)
        seen.append((side, label))
        labels.append(label)
    return labels + [None] * (3 - len(labels))


def _follow(r: Relaxer, point: Settled, history: list[str]
            ) -> tuple[Settled, tuple[Settled, Settled] | None]:
    """Up to max_mode_follow cycles from a saddle, along its modes below -saddle_cm1
    (modes.off_saddle, Relaxer.descend): ± along the one of a first-order saddle, its QRC step (a
    TS candidate when both sides are minima distinct as labelled); one side along all of a
    higher-order one (± would mostly stop at first-order saddles)."""
    for cycle in range(1, r.max_mode_follow + 1):
        if point.tier != "saddle" or point.freq is None or not point.freq.imaginary_modes:
            break
        order = sum(nu < -r.gates.saddle_cm1 for nu in point.freq.frequencies_cm1 or ())
        step = off_saddle(point.freq, r.template.xyz.symbols, below_cm1=r.gates.saddle_cm1)
        sides = r.descend(point.x, point.freq, step, both=order == 1, name=f"follow{cycle}")
        labels = _labels(point, sides)
        verdict = classify_mode_follow(*labels)
        history.append(f"follow{cycle}:{verdict}")
        if verdict == "same_as_source":
            break
        reached = [s for s, label in zip(sides, labels[1:], strict=False)
                   if isinstance(s, Settled) and label != "source"]
        if verdict == "ts_candidate" and all(s.tier != "saddle" for s in reached):
            return point, (reached[0], reached[1])
        point = min(reached, key=lambda s: (s.tier == "saddle", s.opt.energy_hartree))
    return point, None


def _outcome(point: Settled, history: Sequence[str],
             pair: tuple[Settled, Settled] | None = None) -> MinimumOutcome:
    """A known basin, a saddle as it is, or a minimum; a TS candidate with its sides."""
    sides = None if pair is None else (_outcome(pair[0], pair[0].history),
                                       _outcome(pair[1], pair[1].history))
    if point.basin is not None:
        return MinimumOutcome("known", point.opt, None, tuple(history), known_basin=point.basin,
                              ts_candidate=sides)
    status: Status = "saddle" if point.tier == "saddle" else "minimum"
    return MinimumOutcome(status, point.opt, point.freq, tuple(history), ts_candidate=sides,
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
    load_xyz: LoadXYZ | None = None,
    resolve: Resolve | None = None,
) -> MinimumOutcome:
    """opt (with init_hessian) → Relaxer.settle (a known basin: no freq) → mode-follow (saddle).

    A converged ``opt`` is taken as it is (``mol`` then gives only charge and multiplicity).
    ``load_xyz`` (default ``known.load_xyz``) reads the optimized geometry for the freq job;
    ``resolve`` opens a freq's Hessian for the certification (not judged without it).
    """
    load = load_xyz or (known.load_xyz if known is not None else None)
    if load is None:
        raise ValueError("relax_to_minimum needs load_xyz or a known Registry")
    find = None if known is None else (lambda ev, x: known.find(ev, coords=x))
    r = Relaxer(mol, method, qm, load, gates, resolve, find, max_mode_follow=max_mode_follow)
    history = ["opt" if opt is None else "opt:reused"]
    done = opt if opt is not None else qm.optimize(mol, method, init_hessian=init_hessian)
    if isinstance(done, Failure):
        return _failed(history, done)
    point = r.settle(done)
    if isinstance(point, Failure):
        return _failed(history, point, done)
    history += point.history
    point, pair = _follow(r, point, history)
    return _outcome(point, history, pair)


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

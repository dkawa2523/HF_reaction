"""Reaction-case context and the saddle and TS actions (design §7.3; the path actions are in
``paths``, the connection and intermediate actions in ``connection``). Jobs run through the
capability Protocols and so the JobStore (idempotent); an action returns the next CaseState and
keeps the evidence behind it (saddle, TS freq, claims, new basins) in ``Work``."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import TYPE_CHECKING, Literal, NamedTuple

import numpy as np

from hfauto.chemistry import profile, topology
from hfauto.chemistry.gates import barrier_verdict, is_first_order_saddle, qrc_drop, spin_ok
from hfauto.chemistry.geometry import declared_coordinate_gradient, most_changed_dihedral
from hfauto.chemistry.interpolation import align_mapped
from hfauto.chemistry.modes import TARGET_HARTREE, amplitude, displace, overlap
from hfauto.chemistry.xyz import (
    XYZ,
    Molecule,
    read_xyz_trajectory,
    write_xyz_trajectory,
    written_geometry,
)
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Evidence, Failure, FileRef, Geometry
from hfauto.core.method import Deadline
from hfauto.core.records import (
    BarrierVerdict,
    ConnectionClaim,
    MinimumRecord,
    ReactionRecord,
    SaddleClaim,
    SpeciesRecord,
)
from hfauto.drivers.minimum import calc_id
from hfauto.drivers.reaction_case.state import CaseRules, CaseState, Decision, Seed

if TYPE_CHECKING:
    from hfauto.drivers.reaction_case.driver import CaseRuntime

_MIN_OVERLAP = 0.3  # an xTB Hessian is used when a negative mode lies along the direction
_MODE_SEEDS = frozenset({"screen_ts", "discovery_ts", "higher_order_retry"})  # tangent: TS mode


class Profile(NamedTuple):
    """A DFT profile between the DFT minima: its frames and energies, the minima's at its ends
    (SCREEN's NEB or IDPP, or a string)."""

    frames: list[np.ndarray]
    energies: tuple[float, ...]


@dataclass
class Work:
    """Evidence behind the CaseState and the artifacts the case adds."""

    path: Profile | None = None  # the latest DFT profile (CaseState.screen), FIND_PATH's start
    saddle: Evidence | None = None
    ts_freq: Evidence | None = None
    connection: ConnectionClaim | None = None
    intermediate: tuple[MinimumRecord, SpeciesRecord] | None = None  # basin, structure reached
    split_ts: tuple[int, str] | None = None  # (split child, saddle calc) it validates (GEN-05)
    calcs: dict[str, Evidence] = field(default_factory=dict)
    species: dict[str, SpeciesRecord] = field(default_factory=dict)
    minima: dict[str, MinimumRecord] = field(default_factory=dict)


@dataclass
class Ctx:
    case: ReactionRecord
    rt: CaseRuntime
    rules: CaseRules
    deadline: Deadline
    folder: Path  # cases/<reaction_id>
    symbols: list[str]
    charge: int
    multiplicity: int
    raw: tuple[np.ndarray, np.ndarray]  # DFT minima in their own frames (JobStore reuse)
    ends: tuple[np.ndarray, np.ndarray]  # the same with b aligned onto a (path endpoints)
    energies: tuple[float, float]  # DFT minimum energies
    log: Callable[[dict[str, object]], None]
    work: Work = field(default_factory=Work)

    @property
    def resolution(self) -> float:
        return self.rules.gates.resolution_kcal / HARTREE_TO_KCAL_MOL

    def mol(self, coords: np.ndarray) -> Molecule:
        xyz = XYZ(list(self.symbols), np.asarray(coords, dtype=float).reshape(-1, 3))
        return Molecule(xyz, self.charge, self.multiplicity)

    def coords(self, geometry: Geometry) -> np.ndarray:
        return np.asarray(self.rt.load_xyz(geometry).coords, dtype=float)

    def note(self, text: str) -> None:
        self.log({"note": text})

    def keep(self, ev: Evidence) -> Evidence:
        self.work.calcs[calc_id(ev)] = ev
        return ev

    def sp(self, coords: np.ndarray) -> float | None:
        return self.sps([coords])[0]

    def sps(self, frames: Sequence[np.ndarray]) -> list[float | None]:
        """DFT SPs at ``frames``, independent jobs run at once (CaseRuntime.map); kept and
        noted in order."""
        rt = self.rt
        return [self._energy(ev) for ev in rt.map(
            lambda x: rt.qm.energy(self.mol(x), rt.method, deadline=self.deadline), frames)]

    def _energy(self, ev: Evidence | Failure) -> float | None:
        if isinstance(ev, Failure):
            self.note(f"sp:{ev.kind.value}")
            return None
        return self.keep(ev).energy_hartree

    def geometry(self, name: str, coords: np.ndarray) -> Geometry:
        return written_geometry(self.mol(coords).write(self.folder / f"{name}.xyz"),
                                self.rt.file_ref)

    def path_file(self, name: str, frames: Sequence[np.ndarray]) -> FileRef:
        images, path = [self.mol(f).xyz for f in frames], self.folder / f"{name}.xyz"
        return self.rt.file_ref(write_xyz_trajectory(images, path))

    def frames(self, ref: FileRef) -> list[np.ndarray]:
        images = read_xyz_trajectory(self.rt.resolve(ref))
        return [np.asarray(i.coords, dtype=float) for i in images]

    def verdict(self, frames: list[np.ndarray], inner: Sequence[float],
                source: Literal["screen", "string"]) -> BarrierVerdict:
        """Class of the sequence-aligned path ``frames`` between the DFT minima (``inner``: its
        interior DFT energies), kept as the latest profile; its maximum bounds the saddle from
        above. A barrierless class is accepted only after densifying: DFT SPs at the midpoints
        of the two segments beside the highest interior node, where a barrier the nodes step
        over would rise (nodes lie ~0.2 A apart, so the midpoints stay near the path)."""
        energies = [self.energies[0], *inner, self.energies[1]]
        self.work.path = Profile(frames, tuple(energies))
        verdict = barrier_verdict(energies, source=source, policy=self.rules.gates)
        if verdict.verdict != "barrierless":
            return verdict
        k = 1 + int(np.argmax(inner))
        mids = [0.5 * (frames[i] + frames[i + 1]) for i in (k - 1, k)]
        low, high = self.sps(mids)
        if low is None or high is None:
            return BarrierVerdict(verdict="unavailable", source=source,
                                  reasons=("midpoint_single_point",))
        frames = [*frames[:k], mids[0], frames[k], mids[1], *frames[k + 1:]]
        energies = [*energies[:k], low, energies[k], high, *energies[k + 1:]]
        self.work.path = Profile(frames, tuple(energies))
        verdict = barrier_verdict(energies, source=source, policy=self.rules.gates)
        self.note(f"{source}_midpoints:{verdict.verdict}")
        return verdict

    def direction(self, x: np.ndarray, seed: Seed | None = None) -> np.ndarray:
        """Reaction direction at a seed (design §6): a TS's own imaginary mode, else the path
        tangent (the endpoint chord off a path) on the reaction centre; with no bond change,
        the gradient of the declared coordinate, else of the most changed dihedral."""
        if seed is not None and seed.tangent is not None and seed.source in _MODE_SEEDS:
            return np.asarray(seed.tangent)
        a, b = (align_mapped(x, end) for end in self.ends)
        centre, bonded = topology.reaction_centre(self.symbols, a, b)
        if not centre and (terms := self.case.coordinate or most_changed_dihedral(bonded, a, b)):
            return declared_coordinate_gradient(terms, x)
        t = np.reshape(b - a if seed is None or seed.tangent is None else seed.tangent, (-1, 3))
        rows, v = sorted(centre) or list(range(len(t))), np.zeros_like(t)
        v[rows] = t[rows]
        return v.ravel() / np.linalg.norm(v)

    def record(self, basin_id: str) -> MinimumRecord:
        return next(r for r, _ in self.rt.minima.values() if r.basin_id == basin_id)


def peak_seed(ctx: Ctx, name: str, source: Literal["screen_hei", "path_hei"]) -> Seed | None:
    """The highest peak detected on the latest profile, refined by a parabola."""
    path = ctx.work.path
    peaks = () if path is None else profile.interior_maxima(path.energies, ctx.resolution)
    if path is None or not peaks:
        return None
    k, _, x = profile.hei(path.frames, path.energies, max(peaks, key=path.energies.__getitem__))
    tangent = tuple(profile.tangent(path.frames, round(k)).ravel())
    return Seed(ctx.geometry(name, x), source, tangent)


def xtb_modes(ctx: Ctx, coords: np.ndarray, below_cm1: float
              ) -> tuple[Evidence | None, list[np.ndarray]]:
    """The xTB freq at ``coords`` (None when unavailable) and its modes below -below_cm1."""
    rt = ctx.rt
    if rt.screen_qm is None or rt.screen_method is None:
        return None, []
    freq = rt.screen_qm.frequencies(ctx.mol(coords), rt.screen_method, deadline=ctx.deadline)
    if isinstance(freq, Failure):
        return None, []
    nus = sorted(freq.frequencies_cm1 or ())
    return freq, [np.asarray(m) for m, nu in zip(freq.imaginary_modes, nus, strict=False)
                  if nu < -below_cm1]


def _seed_hessian(ctx: Ctx, seed: Seed, x: np.ndarray, direction: np.ndarray
                  ) -> Evidence | Failure:
    """The seed's own TS freq, else xTB when one of its negative modes lies along the
    direction (chem 20), else DFT at the seed."""
    if seed.hessian is not None:
        ctx.note("saddle_hessian:ts_freq")
        return seed.hessian
    xtb, modes = xtb_modes(ctx, x, ctx.rules.gates.noise_cm1)
    score = max((overlap(m, direction) for m in modes), default=0.0)
    use_xtb = xtb is not None and score >= _MIN_OVERLAP
    ctx.note(f"saddle_hessian:{'xtb' if use_xtb else 'dft'}:overlap:{score:.2f}")
    return xtb if xtb is not None and use_xtb else ctx.rt.qm.frequencies(
        ctx.mol(x), ctx.rt.method, deadline=ctx.deadline)


def refine_saddle(ctx: Ctx, state: CaseState, decision: Decision) -> CaseState:
    """The front seed's reaction direction and initial Hessian → saddle.refine with only that
    direction negative. A search stalled at maxiter restarts once from its last frame with a
    fresh Hessian: a new front seed that is not counted (the attempts count per case, a split
    child starts from 0; only the walltime deadline is shared). A new search drops the claim
    and connection verdict of an earlier TS whose sides joined one basin (row 8)."""
    rt, seed = ctx.rt, state.seeds[0]
    counted = int(seed.source != "saddle_restart")
    state = replace(state, seeds=state.seeds[1:], saddle_attempts=state.saddle_attempts + counted,
                    last_saddle="failed", ts_check=None, claim=None, connection=None)
    x = ctx.coords(seed.geometry)
    direction = ctx.direction(x, seed)
    hessian = _seed_hessian(ctx, seed, x, direction)
    if isinstance(hessian, Failure):
        ctx.note(f"saddle_hessian:{hessian.kind.value}")
        return state
    result = rt.saddle.refine(ctx.mol(x), rt.method, hessian=hessian, mode=tuple(direction),
                              deadline=ctx.deadline)
    if isinstance(result, Failure):
        ctx.note(f"saddle:{result.kind.value}:{result.reason}")
        if result.final is None or seed.source == "saddle_restart":
            return state
        restart = Seed(result.final, "saddle_restart", tuple(direction))
        return replace(state, seeds=(restart, *state.seeds))
    ctx.work.saddle = result
    return replace(state, last_saddle="converged")


def mode_amplitude(ctx: Ctx, freq: Evidence, index: int) -> float:
    """modes.amplitude of imaginary mode ``index`` for max(3 x QRC drop, TARGET_HARTREE)."""
    target = max(3.0 * qrc_drop(freq.level), TARGET_HARTREE)
    return amplitude((freq.frequencies_cm1 or ())[index], np.asarray(freq.imaginary_modes[index]),
                     ctx.symbols, target_hartree=target)


def _pushed(ctx: Ctx, freq: Evidence, x: np.ndarray, name: str) -> Seed:
    """A higher-order saddle pushed once along its most negative mode other than the reaction
    mode (the one along the direction); its verified freq is the new seed's Hessian."""
    direction = ctx.direction(x)
    r = max(range(len(freq.imaginary_modes)),
            key=lambda i: overlap(np.asarray(freq.imaginary_modes[i]), direction))
    other = 1 if r == 0 else 0  # imaginary modes are ordered by frequency
    step = mode_amplitude(ctx, freq, other)
    pushed, _ = displace(x, np.asarray(freq.imaginary_modes[other]), step)
    return Seed(ctx.geometry(name, pushed), "higher_order_retry", freq.imaginary_modes[r], freq)


def validate_ts(ctx: Ctx, state: CaseState, decision: Decision) -> CaseState:
    """Separate DFT freq on the saddle → is_first_order_saddle and spin_ok."""
    rt, saddle, gates = ctx.rt, ctx.work.saddle, ctx.rules.gates
    if saddle is None:
        return replace(state, last_saddle="failed")
    x = ctx.coords(saddle.final)
    freq = rt.qm.frequencies(ctx.mol(x), rt.method, scf_guess=saddle, deadline=ctx.deadline)
    if isinstance(freq, Failure):
        ctx.note(f"ts_freq:{freq.kind.value}")
        return replace(state, last_saddle="failed")
    gate = is_first_order_saddle(freq, saddle=saddle, policy=gates)
    if gate:
        ctx.work.ts_freq = ctx.keep(freq)
        notes = (*gate.notes, *spin_ok(freq, gates).reasons)
        claim = SaddleClaim(saddle_calc=calc_id(ctx.keep(saddle)), freq_calc=calc_id(freq),
                            imag_cm1=min(freq.frequencies_cm1 or (0.0,)),
                            energy_hartree=freq.energy_hartree, notes=notes)
        return replace(state, ts_check="ok", claim=claim, connection_attempts=0)
    ctx.note(f"ts_rejected:{','.join(gate.reasons)}")
    if "higher_order" in gate.reasons:
        seed = _pushed(ctx, freq, x, f"retry{state.saddle_attempts}")
        return replace(state, last_saddle="failed", seeds=(seed, *state.seeds))
    if "no_imaginary_mode" in gate.reasons:
        return replace(state, ts_check="collapsed")
    return replace(state, last_saddle="failed")

"""Reaction-case actions (design §7.3). Jobs run through the capability Protocols and so the
JobStore (idempotent); an action returns the next CaseState and keeps the evidence behind it
(saddle, TS freq, claims, new basins) in ``Work``."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import TYPE_CHECKING, Literal, NamedTuple

import numpy as np

from hfauto.chemistry import profile, topology
from hfauto.chemistry import xyz_trajectory as trajectory
from hfauto.chemistry.gates import (
    barrier_verdict,
    connection,
    is_first_order_saddle,
    qrc_drop,
    spin_ok,
)
from hfauto.chemistry.geometry import (
    declared_coordinate,
    declared_coordinate_gradient,
    most_changed_dihedral,
)
from hfauto.chemistry.identity import mapped_equivalent, periodic_nearest
from hfauto.chemistry.interpolation import align_mapped, align_sequential, idpp
from hfauto.chemistry.modes import amplitude, displace, overlap
from hfauto.chemistry.xyz import XYZ, Molecule, composition_key, geometry_fingerprint
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Evidence, Failure, FileRef, Geometry
from hfauto.core.ids import species_id
from hfauto.core.method import Deadline
from hfauto.core.records import (
    BarrierVerdict,
    ConnectionClaim,
    MinimumRecord,
    ReactionRecord,
    SaddleClaim,
    SpeciesRecord,
)
from hfauto.drivers.minimum import calc_id, relax_to_minimum
from hfauto.drivers.reaction_case import state as case_state
from hfauto.drivers.reaction_case.state import Action, CasePolicy, CaseState, Decision, Seed

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
    intermediate: tuple[MinimumRecord, SpeciesRecord] | None = None
    calcs: dict[str, Evidence] = field(default_factory=dict)
    species: dict[str, SpeciesRecord] = field(default_factory=dict)
    minima: dict[str, MinimumRecord] = field(default_factory=dict)


@dataclass
class Ctx:
    case: ReactionRecord
    rt: CaseRuntime
    policy: CasePolicy
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
        return self.policy.gates.resolution_kcal / HARTREE_TO_KCAL_MOL

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
        ev = self.rt.qm.energy(self.mol(coords), self.rt.method, deadline=self.deadline)
        if isinstance(ev, Failure):
            self.note(f"sp:{ev.kind.value}")
            return None
        return self.keep(ev).energy_hartree

    def geometry(self, name: str, coords: np.ndarray) -> Geometry:
        path = self.mol(coords).write(self.folder / f"{name}.xyz")
        xyz = trajectory.read_xyz_trajectory(path)[0]  # fingerprint of the file as written
        return Geometry(file=self.rt.file_ref(path), symbols=tuple(xyz.symbols),
                        fingerprint=geometry_fingerprint(xyz.symbols, xyz.coords))

    def path_file(self, name: str, frames: Sequence[np.ndarray]) -> FileRef:
        images, path = [self.mol(f).xyz for f in frames], self.folder / f"{name}.xyz"
        return self.rt.file_ref(trajectory.write_xyz_trajectory(images, path))

    def frames(self, ref: FileRef) -> list[np.ndarray]:
        images = trajectory.read_xyz_trajectory(self.rt.resolve(ref))
        return [np.asarray(i.coords, dtype=float) for i in images]

    def verdict(self, frames: list[np.ndarray], inner: Sequence[float],
                source: Literal["screen", "string"]) -> BarrierVerdict:
        """Class of the path ``frames`` between the DFT minima (``inner``: its interior DFT
        energies), kept as the latest profile; its maximum bounds the saddle from above."""
        energies = (self.energies[0], *inner, self.energies[1])
        self.work.path = Profile(frames, energies)
        return barrier_verdict(energies, source=source, policy=self.policy.gates,
                               max_node_spacing_A=profile.max_node_spacing(frames))

    def direction(self, x: np.ndarray, seed: Seed | None = None) -> np.ndarray:
        """Reaction direction at a seed (design §6): a TS's own imaginary mode, else the path
        tangent (the endpoint chord off a path) on the reaction centre; with no bond change,
        the gradient of the declared coordinate, else of the most changed dihedral."""
        if seed is not None and seed.tangent is not None and seed.source in _MODE_SEEDS:
            return np.asarray(seed.tangent)
        a, b = (align_mapped(x, end) for end in self.ends)
        centre, bonded = _reaction_centre(self.symbols, a, b)
        if not centre and (terms := self.case.coordinate or most_changed_dihedral(bonded, a, b)):
            return declared_coordinate_gradient(terms, x)
        t = np.reshape(b - a if seed is None or seed.tangent is None else seed.tangent, (-1, 3))
        rows, v = sorted(centre) or list(range(len(t))), np.zeros_like(t)
        v[rows] = t[rows]
        return v.ravel() / np.linalg.norm(v)

    def record(self, basin_id: str) -> MinimumRecord:
        return next(r for r, _ in self.rt.minima.values() if r.basin_id == basin_id)


def _reaction_centre(symbols: list[str], a: np.ndarray, b: np.ndarray
                     ) -> tuple[set[int], frozenset[tuple[int, int]]]:
    """Atoms whose bonds change from a to b and their neighbours; the bonds of a."""
    formed, broken = topology.bond_changes(symbols, a, b)
    changed, bonded = {i for pair in formed | broken for i in pair}, topology.bonds(symbols, a)
    return changed | {k for pair in bonded if changed & set(pair) for k in pair}, bonded


def _unavailable(reason: str) -> BarrierVerdict:
    return BarrierVerdict(verdict="unavailable", source="screen", reasons=(reason,))


def _peak_seed(ctx: Ctx, name: str, source: Literal["screen_hei", "path_hei"]) -> Seed | None:
    """The highest peak detected on the latest profile, refined by a parabola."""
    path = ctx.work.path
    peaks = () if path is None else profile.interior_maxima(path.energies, ctx.resolution)
    if path is None or not peaks:
        return None
    k, _, x = profile.hei(path.frames, path.energies, max(peaks, key=path.energies.__getitem__))
    tangent = tuple(profile.tangent(path.frames, round(k)).ravel())
    return Seed(ctx.geometry(name, x), source, tangent)


def _xtb_modes(ctx: Ctx, coords: np.ndarray, below_cm1: float
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


def _xtb_ts_mode(ctx: Ctx, coords: np.ndarray) -> tuple[float, ...] | None:
    """A low-level TS's imaginary mode: its xTB freq has exactly one mode below -saddle_cm1."""
    _, modes = _xtb_modes(ctx, coords, ctx.policy.gates.saddle_cm1)
    return tuple(modes[0]) if len(modes) == 1 else None


def _shortcut(ctx: Ctx, geometry: Geometry) -> tuple[BarrierVerdict, Seed | None] | None:
    """Step 1: a low-level TS already known (discovery / mode-follow): three points."""
    x = ctx.coords(geometry)
    mode = _xtb_ts_mode(ctx, x)
    if mode is None:
        ctx.note("low_level_ts_rejected")
        return None
    e_ts = ctx.sp(x)
    if e_ts is None:
        return None
    verdict = barrier_verdict((ctx.energies[0], e_ts, ctx.energies[1]), source="screen",
                              policy=ctx.policy.gates)
    return verdict, Seed(geometry, "discovery_ts", mode)


def _neb(ctx: Ctx, frames: list[np.ndarray]) -> tuple[list[np.ndarray], Geometry | None]:
    """The xTB CI-NEB from ``frames`` (fixed ends) and its TS; ``frames`` without a NEB."""
    engine, method = ctx.rt.screen_path, ctx.rt.screen_method
    if engine is None or method is None:
        return frames, None
    neb = engine.find_path(ctx.mol(frames[0]), ctx.mol(frames[-1]), method, images=len(frames),
                           initial_path=ctx.path_file("screen_idpp", frames),
                           deadline=ctx.deadline)
    if isinstance(neb, Failure):
        ctx.note(f"screen_neb:{neb.kind.value}")
        return frames, None
    return align_sequential(ctx.frames(neb.images)), neb.ts


def _screen_path(ctx: Ctx) -> tuple[BarrierVerdict, Seed | None]:
    """Steps 2-3: the NEB from the IDPP between the DFT minima, DFT SPs on its interior; a single
    peak seeds at the NEB's TS when its xTB freq confirms it, else at the DFT peak."""
    try:
        frames, ts = _neb(ctx, idpp(ctx.symbols, *ctx.ends, ctx.policy.screen_images))
    except ValueError as exc:
        return _unavailable(f"idpp:{exc}"), None
    inner = [ctx.sp(f) for f in frames[1:-1]]
    if any(e is None for e in inner):
        return _unavailable("screen_single_point"), None
    verdict = ctx.verdict(frames, [e for e in inner if e is not None], "screen")
    if verdict.verdict != "single":
        return verdict, None
    mode = None if ts is None else _xtb_ts_mode(ctx, ctx.coords(ts))
    if ts is not None and mode is not None:
        return verdict, Seed(ts, "screen_ts", mode)
    return verdict, _peak_seed(ctx, "screen_hei", "screen_hei")


def screen(ctx: Ctx, state: CaseState, decision: Decision) -> CaseState:
    """Barrier pre-check: shortcut from a low-level TS, else the low-level path (chem 16)."""
    found = None if ctx.case.low_level_ts is None else _shortcut(ctx, ctx.case.low_level_ts)
    if found is None or found[0].verdict != "single":
        found = _screen_path(ctx)
    verdict, seed = found
    ctx.note(f"screen:{verdict.verdict}:{','.join(verdict.reasons)}")
    return case_state.record_profile(state, verdict, seed)


def _seed_hessian(ctx: Ctx, seed: Seed, x: np.ndarray, direction: np.ndarray
                  ) -> Evidence | Failure:
    """The seed's own TS freq, else xTB when one of its negative modes lies along the
    direction (chem 20), else DFT at the seed."""
    if seed.hessian is not None:
        ctx.note("saddle_hessian:ts_freq")
        return seed.hessian
    xtb, modes = _xtb_modes(ctx, x, ctx.policy.gates.noise_cm1)
    score = max((overlap(m, direction) for m in modes), default=0.0)
    use_xtb = xtb is not None and score >= _MIN_OVERLAP
    ctx.note(f"saddle_hessian:{'xtb' if use_xtb else 'dft'}:overlap:{score:.2f}")
    return xtb if xtb is not None and use_xtb else ctx.rt.qm.frequencies(
        ctx.mol(x), ctx.rt.method, deadline=ctx.deadline)


def refine_saddle(ctx: Ctx, state: CaseState, decision: Decision) -> CaseState:
    """The front seed's reaction direction and initial Hessian → saddle.refine with only that
    direction negative. A search stalled at maxiter restarts once from its last frame with a
    fresh Hessian (a new front seed, same attempt budget)."""
    rt, seed = ctx.rt, state.seeds[0]
    state = replace(state, seeds=state.seeds[1:], saddle_attempts=state.saddle_attempts + 1,
                    last_saddle="failed", ts_check=None)
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


def _amplitude(ctx: Ctx, freq: Evidence, index: int) -> float:
    """modes.amplitude of imaginary mode ``index`` for max(3 x QRC drop, qrc_target_hartree)."""
    p = ctx.policy
    target = max(3.0 * qrc_drop(freq.level, p.gates), p.qrc_target_hartree)
    return amplitude((freq.frequencies_cm1 or ())[index], np.asarray(freq.imaginary_modes[index]),
                     ctx.symbols, target_hartree=target, bounds_A=p.qrc_bounds_A)


def _pushed(ctx: Ctx, freq: Evidence, x: np.ndarray, name: str) -> Seed:
    """A higher-order saddle pushed once along its most negative mode other than the reaction
    mode (the one along the direction); its verified freq is the new seed's Hessian."""
    direction = ctx.direction(x)
    r = max(range(len(freq.imaginary_modes)),
            key=lambda i: overlap(np.asarray(freq.imaginary_modes[i]), direction))
    other = 1 if r == 0 else 0  # imaginary modes are ordered by frequency
    pushed, _ = displace(x, np.asarray(freq.imaginary_modes[other]), _amplitude(ctx, freq, other))
    return Seed(ctx.geometry(name, pushed), "higher_order_retry", freq.imaginary_modes[r], freq)


def validate_ts(ctx: Ctx, state: CaseState, decision: Decision) -> CaseState:
    """Separate DFT freq on the saddle → is_first_order_saddle and spin_ok."""
    rt, saddle, gates = ctx.rt, ctx.work.saddle, ctx.policy.gates
    if saddle is None:
        return replace(state, last_saddle="failed")
    x = ctx.coords(saddle.final)
    freq = rt.qm.frequencies(ctx.mol(x), rt.method, deadline=ctx.deadline)
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


def _register(ctx: Ctx, coords: np.ndarray, name: str,
              source: Literal["connection", "intermediate"]) -> MinimumRecord | None:
    """relax_to_minimum → Registry: the known basin, a new one, or None (no minimum)."""
    rt = ctx.rt
    out = relax_to_minimum(ctx.mol(coords), rt.method, rt.qm, known=rt.registry,
                           deadline=ctx.deadline, gates=ctx.policy.gates)
    if out.status == "known" and out.known_basin is not None:
        return ctx.record(out.known_basin)
    if out.status not in ("minimum", "soft_minimum") or out.opt is None or out.freq is None:
        ctx.note(f"{name}:{out.status}:{out.failure.reason if out.failure else ''}")
        return None
    x, opt = ctx.coords(out.opt.final), out.opt
    species = SpeciesRecord(
        species_id=species_id(ctx.case.reaction_id, name),
        composition_id=composition_key(ctx.symbols, ctx.charge, ctx.multiplicity),
        charge=ctx.charge, multiplicity=ctx.multiplicity, geometry=opt.final, source=source,
        state_label=topology.state_label(ctx.symbols, x), energy_hartree=opt.energy_hartree,
        level_key=opt.level.full_key(),
    )
    record = rt.registry.add(out, species, tier="dft")
    ctx.keep(opt)
    ctx.keep(out.freq)
    ctx.work.species[species.species_id] = rt.species[species.species_id] = species
    ctx.work.minima[record.minimum_id] = record
    geometry = rt.minima[record.minimum_id][1] if record.minimum_id in rt.minima else opt.final
    rt.minima[record.minimum_id] = (record, geometry)
    return record


def _assign(ctx: Ctx, side: Evidence, x: np.ndarray, name: str) -> str | None:
    """Registry.find; a torsional case falls back to the nearest declared dihedral (CH-04)."""
    basin = ctx.rt.registry.find(side)
    if basin is not None:
        return ctx.record(basin).minimum_id
    terms = ctx.case.coordinate
    if ctx.case.torsional and terms and all(t.kind == "dihedral" for t in terms):
        values = [declared_coordinate(terms, end) for end in ctx.raw]
        return ctx.case.minima[periodic_nearest(declared_coordinate(terms, x), values)]
    record = _register(ctx, x, name, "connection")
    return None if record is None else record.minimum_id


def connect(ctx: Ctx, state: CaseState, decision: Decision) -> CaseState:
    """QRC: displace ± along the TS mode by an energy target, optimize, assign, gate."""
    rt, p, freq = ctx.rt, ctx.policy, ctx.work.ts_freq
    attempt = state.connection_attempts + 1
    state = replace(state, connection_attempts=attempt)
    if freq is None or not freq.imaginary_modes:
        return replace(state, connection="failed")
    step = min(_amplitude(ctx, freq, 0) * p.qrc_retry_factor ** (attempt - 1), p.qrc_bounds_A[1])
    mode = np.asarray(freq.imaginary_modes[0])
    plus, minus = (rt.qm.optimize(ctx.mol(y), rt.method, init_hessian=freq, deadline=ctx.deadline)
                   for y in displace(ctx.coords(freq.final), mode, step))  # the TS Hessian
    if isinstance(plus, Failure) or isinstance(minus, Failure):
        ctx.note(f"qrc{attempt}:side_failed")
        return replace(state, connection="failed")
    finals = [ctx.coords(plus.final), ctx.coords(minus.final)]
    first, second = (_assign(ctx, s, x, f"qrc{attempt}_{i}")
                     for i, (s, x) in enumerate(zip((plus, minus), finals, strict=True)))
    distinct = mapped_equivalent(ctx.symbols, *finals) if ctx.case.degenerate else True
    bonds = tuple(topology.bonds(ctx.symbols, x) for x in (*ctx.ends, finals[1], finals[0]))
    gate, label = connection(freq, (plus, minus), (first, second), frozenset(ctx.case.minima),
                             degenerate=ctx.case.degenerate, sides_distinct=distinct,
                             bond_sets=bonds, policy=p.gates)
    if label == "failed" or first is None or second is None:
        ctx.note(f"qrc{attempt}:{','.join(gate.reasons)}")
        retry = "sides_same_basin" in gate.reasons  # the only failure worth a wider displacement
        return replace(state, connection="same_basin" if retry else "failed")
    sides = (calc_id(ctx.keep(plus)), calc_id(ctx.keep(minus)))
    ctx.work.connection = ConnectionClaim(side_calcs=sides, minima=(first, second),
                                          amplitude_A=step)
    return replace(state, connection=label)


def _past_the_well(ctx: Ctx, state: CaseState, path: Profile, name: str) -> CaseState:
    """The well relaxed into an endpoint: barrierless when no interior point rises a resolution
    above the higher end, else the highest peak seeds the saddle search."""
    v, e = state.screen, path.energies
    if v is not None and max(e[1:-1]) - max(e[0], e[-1]) < ctx.resolution:
        return replace(state, screen=v.model_copy(update={"verdict": "barrierless",
                                                          "reasons": ("well_is_endpoint",)}))
    source = "screen_hei" if v is not None and v.source == "screen" else "path_hei"
    seed = _peak_seed(ctx, f"{name}_hei", source)
    return replace(state, seeds=(*state.seeds, seed)) if seed else state


def validate_intermediate(ctx: Ctx, state: CaseState, decision: Decision) -> CaseState:
    """relax_to_minimum → Registry on the collapsed saddle (or soft TS whose QRC failed: its
    claim is withdrawn) or on the latest profile's lowest well; consumes its trigger."""
    path = ctx.work.path if decision.reason == "path_intermediate" else None
    state = replace(state, ts_check=None, last_saddle=None, claim=None, connection=None)
    if path is not None:
        wells = profile.interior_maxima([-e for e in path.energies], ctx.resolution)
        x = path.frames[min(wells, key=path.energies.__getitem__)]
    elif ctx.work.saddle is not None:
        x = ctx.coords(ctx.work.saddle.final)
    else:
        return replace(state, intermediate="same_as_endpoint")
    name = f"int{state.saddle_attempts}_{state.path_runs}"
    record = _register(ctx, x, name, "intermediate")
    ends = {ctx.rt.minima[m][0].basin_id for m in ctx.case.minima}
    if record is not None and record.basin_id not in ends:
        ctx.work.intermediate = (record, ctx.rt.species[record.species_id])
        return replace(state, intermediate="distinct")
    state = replace(state, intermediate="same_as_endpoint")
    return state if path is None else _past_the_well(ctx, state, path, name)


def _initial_path(ctx: Ctx, beads: int, name: str) -> FileRef | None:
    """The latest DFT profile's path (SCREEN's or the last chunk's), else a new IDPP, resampled
    between the DFT minima and aligned in sequence (freezeN moves ends with rigid jumps)."""
    frames = None if ctx.work.path is None else ctx.work.path.frames
    if frames is None:
        try:
            frames = idpp(ctx.symbols, *ctx.ends, beads)
        except ValueError as exc:
            ctx.note(f"idpp:{exc}")
            return None
    images = trajectory.resample_xyz_trajectory([ctx.mol(f).xyz for f in frames], beads)
    inner = [np.asarray(i.coords, dtype=float) for i in images[1:-1]]
    return ctx.path_file(name, align_sequential([ctx.ends[0], *inner, ctx.ends[1]]))


def find_path(ctx: Ctx, state: CaseState, decision: Decision) -> CaseState:
    """One DFT string chunk, classified at once (unconverged too: its maximum bounds the saddle
    and its peak is a seed); a single peak seeds the saddle search. The path gets the DFT
    minima back at its ends, which NWChem's freezeN moves."""
    rt, p, run_name = ctx.rt, ctx.policy, f"string{state.path_runs}"
    initial = _initial_path(ctx, p.string_beads, f"{run_name}_initial")
    if initial is None:
        return replace(state, path_runs=state.path_runs + 1)
    run = rt.path.find_path(ctx.mol(ctx.ends[0]), ctx.mol(ctx.ends[1]), rt.method,
                            images=p.string_beads, initial_path=initial, deadline=ctx.deadline)
    if isinstance(run, Failure):
        ctx.note(f"{run_name}:{run.kind.value}")
        return replace(state, path_runs=state.path_runs + 1)
    frames = align_sequential([ctx.ends[0], *ctx.frames(run.images)[1:-1], ctx.ends[1]])
    verdict = ctx.verdict(frames, run.energies_hartree[1:-1], "string")
    seed = _peak_seed(ctx, f"{run_name}_hei", "path_hei") if verdict.verdict == "single" else None
    return case_state.record_profile(state, verdict, seed)


Handler = Callable[[Ctx, CaseState, Decision], CaseState]
HANDLERS: dict[Action, Handler] = {
    Action.SCREEN: screen,
    Action.REFINE_SADDLE: refine_saddle,
    Action.VALIDATE_TS: validate_ts,
    Action.FIND_PATH: find_path,
    Action.CONNECT: connect,
    Action.VALIDATE_INTERMEDIATE: validate_intermediate,
}

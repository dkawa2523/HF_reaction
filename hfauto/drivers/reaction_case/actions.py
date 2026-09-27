"""Reaction-case actions (design §7.3). Jobs run through the capability Protocols and so the
JobStore (idempotent); an action returns the next CaseState and keeps the evidence behind it
(saddle, TS freq, claims, new basins) in ``Work``."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import TYPE_CHECKING, Literal, NamedTuple

import numpy as np

from hfauto.chemistry import profile
from hfauto.chemistry import xyz_trajectory as trajectory
from hfauto.chemistry.gates import (
    barrier_verdict,
    connection,
    is_first_order_saddle,
    qrc_drop,
    spin_ok,
)
from hfauto.chemistry.identity import mapped_equivalent, mapped_rmsd, periodic_nearest
from hfauto.chemistry.interpolation import align_mapped, align_sequential, idpp
from hfauto.chemistry.modes import displace, overlap, qrc_amplitude
from hfauto.chemistry.topology import declared_coordinate, declared_coordinate_gradient, state_label
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

_COLLAPSE_A = 0.05  # xTB endpoints this close (identity mapping) collapsed into one structure
_MIN_OVERLAP = 0.3  # mode vs tangent: use the xTB Hessian / below it only a note
_RETRY_A = 0.1  # step down the second imaginary mode of a higher-order saddle
_SETTLED = 0.1  # a string ends once its bead energies move < 0.1 x resolution in 3 iterations


class Profile(NamedTuple):
    """A DFT profile with the DFT minima energies at its ends. ``exact`` when its end frames
    are the DFT minima themselves (IDPP, string), not the xTB minima a GS path starts from."""

    frames: list[np.ndarray]
    energies: tuple[float, ...]
    exact: bool


@dataclass
class Work:
    """Evidence behind the CaseState and the artifacts the case adds."""

    initial: list[np.ndarray] | None = None  # screen GS path or screen IDPP (row 15)
    path: Profile | None = None  # the latest DFT profile, behind CaseState.screen
    saddle: Evidence | None = None
    ts_freq: Evidence | None = None
    saddle_notes: tuple[str, ...] = ()
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
                source: Literal["screen", "string"], *, exact: bool = True) -> BarrierVerdict:
        """Class of the profile along ``frames`` (``inner`` between the DFT minima energies),
        kept as the latest path. Only a path between the DFT minima themselves bounds the
        saddle: a barrierless GS path (xTB ends) leaves the question to FIND_PATH."""
        energies = (self.energies[0], *inner, self.energies[1])
        self.work.path = Profile(frames, energies, exact)
        verdict = barrier_verdict(energies, source=source, policy=self.policy.gates,
                                  max_node_spacing_A=profile.max_node_spacing(frames))
        if verdict.verdict != "barrierless" or exact:
            return verdict
        return verdict.model_copy(update={"verdict": "unavailable", "reasons": ("low_level_ends",)})

    def chord(self, coords: np.ndarray) -> tuple[float, ...]:
        """Endpoint difference aligned onto ``coords``: the tangent of a seed off a path."""
        a, b = (align_mapped(coords, end) for end in self.ends)
        return tuple((b - a).ravel())

    def direction(self, coords: np.ndarray, tangent: tuple[float, ...] | None
                  ) -> np.ndarray | None:
        """Gradient of the declared coordinate at ``coords``, else the seed tangent."""
        if self.case.coordinate:
            return declared_coordinate_gradient(self.case.coordinate, coords)
        return None if tangent is None else np.asarray(tangent)

    def record(self, basin_id: str) -> MinimumRecord:
        return next(r for r, _ in self.rt.minima.values() if r.basin_id == basin_id)


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


def _xtb_saddle_ok(ctx: Ctx, coords: np.ndarray) -> bool:
    """Low-level freq: exactly one mode below -saddle_cm1, and it is the lowest."""
    rt = ctx.rt
    if rt.screen_qm is None or rt.screen_method is None:
        return False
    freq = rt.screen_qm.frequencies(ctx.mol(coords), rt.screen_method, deadline=ctx.deadline)
    if isinstance(freq, Failure):
        return False
    return sum(nu < -ctx.policy.gates.saddle_cm1 for nu in freq.frequencies_cm1 or ()) == 1


def _shortcut(ctx: Ctx, geometry: Geometry) -> tuple[BarrierVerdict, Seed | None] | None:
    """Step 1: a low-level TS already known (discovery / mode-follow): three points."""
    x = ctx.coords(geometry)
    if not _xtb_saddle_ok(ctx, x):
        ctx.note("low_level_ts_rejected")
        return None
    e_ts = ctx.sp(x)
    if e_ts is None:
        return None
    verdict = barrier_verdict((ctx.energies[0], e_ts, ctx.energies[1]), source="screen",
                              policy=ctx.policy.gates)
    return verdict, Seed(geometry, "discovery_ts", ctx.chord(x))


def _screened(ctx: Ctx, frames: list[np.ndarray], *, exact: bool, ts: Geometry | None = None
              ) -> tuple[BarrierVerdict, Seed | None]:
    """DFT single points on the interior nodes; a single peak seeds at the GS's TSOpt structure
    when its xTB freq confirms it, else at the DFT peak."""
    inner = [ctx.sp(f) for f in frames[1:-1]]
    if any(e is None for e in inner):
        return _unavailable("screen_single_point"), None
    verdict = ctx.verdict(frames, [e for e in inner if e is not None], "screen", exact=exact)
    if verdict.verdict != "single":
        return verdict, None
    if ts is not None and _xtb_saddle_ok(ctx, ctx.coords(ts)):
        return verdict, Seed(ts, "screen_ts", ctx.chord(ctx.coords(ts)))
    return verdict, _peak_seed(ctx, "screen_hei", "screen_hei")


def _screen_idpp(ctx: Ctx) -> tuple[BarrierVerdict, Seed | None]:
    """Step 2 when xTB collapses the endpoints: the DFT IDPP between the DFT minima."""
    try:
        frames = idpp(ctx.symbols, *ctx.ends, ctx.policy.screen_images)
    except ValueError as exc:
        return _unavailable(f"idpp:{exc}"), None
    ctx.work.initial = frames  # FIND_PATH (row 15) starts from it
    return _screened(ctx, frames, exact=True)


def _screen_path(ctx: Ctx) -> tuple[BarrierVerdict, Seed | None]:
    """Steps 2-4: xTB endpoints, then GS with TSOpt (the DFT IDPP when xTB collapses them)."""
    qm, path, method = ctx.rt.screen_qm, ctx.rt.screen_path, ctx.rt.screen_method
    if qm is None or path is None or method is None:
        return _unavailable("no_screen_engines"), None
    ea, eb = (qm.optimize(ctx.mol(x), method, deadline=ctx.deadline) for x in ctx.ends)
    if isinstance(ea, Failure) or isinstance(eb, Failure):
        kind = next(o.kind.value for o in (ea, eb) if isinstance(o, Failure))
        return _unavailable(f"screen_endpoint:{kind}"), None
    low_a, low_b = ctx.coords(ea.final), ctx.coords(eb.final)
    if mapped_rmsd(low_a, low_b) <= _COLLAPSE_A:
        return _screen_idpp(ctx)
    prof = path.find_path(ctx.mol(low_a), ctx.mol(align_mapped(low_a, low_b)), method,
                          images=ctx.policy.screen_images, refine_ts=True, deadline=ctx.deadline)
    if isinstance(prof, Failure):
        return _unavailable(f"screen_path:{prof.kind.value}"), None
    frames = ctx.work.initial = align_sequential(ctx.frames(prof.images))  # dlc images rotate
    return _screened(ctx, frames, exact=False, ts=prof.ts)


def screen(ctx: Ctx, state: CaseState, decision: Decision) -> CaseState:
    """Barrier pre-check: shortcut from a low-level TS, else the low-level path (chem 16)."""
    found = None if ctx.case.low_level_ts is None else _shortcut(ctx, ctx.case.low_level_ts)
    if found is None or found[0].verdict != "single":
        found = _screen_path(ctx)
    verdict, seed = found
    ctx.note(f"screen:{verdict.verdict}:{','.join(verdict.reasons)}")
    return case_state.record_profile(state, verdict, seed)


def _mode_index(freq: Evidence, direction: np.ndarray | None, noise_cm1: float
                ) -> tuple[int | None, float | None]:
    """Imaginary mode (ascending order) most parallel to ``direction`` and that overlap."""
    freqs = freq.frequencies_cm1 or ()
    modes = [m for m, nu in zip(freq.imaginary_modes, sorted(freqs), strict=False)
             if nu < -noise_cm1]
    if direction is None or not modes:
        return None, None
    scores = [overlap(np.asarray(m), direction) for m in modes]
    best = int(np.argmax(scores))
    return best, scores[best]


def _first_hessian(ctx: Ctx, mol: Molecule, direction: np.ndarray | None) -> Evidence | None:
    """xTB Hessian at the seed if it has one negative mode along the tangent (chem 20)."""
    rt = ctx.rt
    if rt.screen_qm is None or rt.screen_method is None:
        return None
    freq = rt.screen_qm.frequencies(mol, rt.screen_method, deadline=ctx.deadline)
    if isinstance(freq, Failure):
        return None
    negative = [nu for nu in freq.frequencies_cm1 or () if nu < -ctx.policy.gates.noise_cm1]
    _, score = _mode_index(freq, direction, ctx.policy.gates.noise_cm1)
    return freq if len(negative) == 1 and score is not None and score >= _MIN_OVERLAP else None


def refine_saddle(ctx: Ctx, state: CaseState, decision: Decision) -> CaseState:
    """Initial Hessian (xTB first when usable, else DFT) → saddle.refine from the front seed."""
    rt, seed = ctx.rt, state.seeds[0]
    state = replace(state, seeds=state.seeds[1:], saddle_attempts=state.saddle_attempts + 1,
                    last_saddle=None, ts_check=None, intermediate=None)
    x = ctx.coords(seed.geometry)
    mol, direction = ctx.mol(x), ctx.direction(x, seed.tangent)
    hessian = _first_hessian(ctx, mol, direction) if state.saddle_attempts == 1 else None
    if hessian is None:  # later attempts, or no usable xTB Hessian: DFT at the seed
        dft = rt.qm.frequencies(mol, rt.method, deadline=ctx.deadline)
        if isinstance(dft, Failure):
            ctx.note(f"saddle_hessian:{dft.kind.value}")
            return replace(state, last_saddle="failed")
        hessian = dft
    index, score = _mode_index(hessian, direction, ctx.policy.gates.noise_cm1)
    low = score is not None and score < _MIN_OVERLAP
    ctx.work.saddle_notes = ("mode_overlap_below_0.3",) if low else ()
    result = rt.saddle.refine(mol, rt.method, hessian=hessian, mode_index=index,
                              deadline=ctx.deadline)
    if isinstance(result, Failure):
        ctx.note(f"saddle:{result.kind.value}:{result.reason}")
        return replace(state, last_saddle="failed")
    ctx.work.saddle = result
    return replace(state, last_saddle="converged")


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
    gate = is_first_order_saddle(freq, saddle=saddle, endpoint_energies=ctx.energies,
                                 torsional=ctx.case.torsional, policy=gates)
    if gate:
        ctx.work.ts_freq = ctx.keep(freq)
        notes = (*gate.notes, *ctx.work.saddle_notes, *spin_ok(freq, gates).reasons)
        claim = SaddleClaim(saddle_calc=calc_id(ctx.keep(saddle)), freq_calc=calc_id(freq),
                            imag_cm1=min(freq.frequencies_cm1 or (0.0,)),
                            energy_hartree=freq.energy_hartree, notes=notes)
        return replace(state, ts_check="ok", claim=claim)
    ctx.note(f"ts_rejected:{','.join(gate.reasons)}")
    if "higher_order" in gate.reasons and len(freq.imaginary_modes) > 1:
        down, _ = displace(x, np.asarray(freq.imaginary_modes[1]), _RETRY_A)
        name = f"retry{state.saddle_attempts}"
        seed = Seed(ctx.geometry(name, down), "higher_order_retry", ctx.chord(down))
        return replace(state, ts_check="higher_order", seeds=(seed, *state.seeds))
    if "no_imaginary_mode" in gate.reasons:
        return replace(state, ts_check="collapsed")
    return replace(state, last_saddle="failed")


def _register(ctx: Ctx, coords: np.ndarray, name: str,
              source: Literal["connection", "intermediate"]) -> MinimumRecord | None:
    """relax_to_minimum → Registry: the known basin, a new one, or None (no minimum)."""
    rt = ctx.rt
    out = relax_to_minimum(ctx.mol(coords), rt.method, rt.qm, known=rt.registry,
                           deadline=ctx.deadline)
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
        state_label=state_label(ctx.symbols, x), energy_hartree=opt.energy_hartree,
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
    if freq is None or freq.hessian is None or not freq.imaginary_modes:
        return replace(state, connection="failed")
    mode, target = np.asarray(freq.imaginary_modes[0]), 3.0 * qrc_drop(freq.level, p.gates)
    first = qrc_amplitude(np.load(rt.resolve(freq.hessian)), mode, bounds_A=p.qrc_bounds_A,
                          target_hartree=max(target, p.qrc_target_hartree))
    amplitude = min(first * p.qrc_retry_factor ** (attempt - 1), p.qrc_bounds_A[1])
    runs = [rt.qm.optimize(ctx.mol(y), rt.method, init_hessian=freq, deadline=ctx.deadline)
            for y in displace(ctx.coords(freq.final), mode, amplitude)]  # the TS Hessian
    plus, minus = runs
    if isinstance(plus, Failure) or isinstance(minus, Failure):
        ctx.note(f"qrc{attempt}:side_failed")
        return replace(state, connection="failed")
    finals = [ctx.coords(plus.final), ctx.coords(minus.final)]
    first, second = (_assign(ctx, s, x, f"qrc{attempt}_{i}")
                     for i, (s, x) in enumerate(zip((plus, minus), finals, strict=True)))
    distinct = mapped_equivalent(ctx.symbols, *finals) if ctx.case.degenerate else True
    gate, label = connection(freq, (plus, minus), (first, second), frozenset(ctx.case.minima),
                             degenerate=ctx.case.degenerate, sides_distinct=distinct,
                             policy=p.gates)
    if label == "failed" or first is None or second is None:
        ctx.note(f"qrc{attempt}:{','.join(gate.reasons)}")
        return replace(state, connection="failed")
    sides = (calc_id(ctx.keep(plus)), calc_id(ctx.keep(minus)))
    ctx.work.connection = ConnectionClaim(side_calcs=sides, minima=(first, second),
                                          amplitude_A=amplitude)
    return replace(state, connection=label)


def _past_the_well(ctx: Ctx, state: CaseState, path: Profile, name: str) -> CaseState:
    """The well relaxed into an endpoint: barrierless (with exact ends) when no interior point
    rises a resolution above the higher end, else the highest peak seeds the saddle search."""
    v, e = state.screen, path.energies
    if v is not None and max(e[1:-1]) - max(e[0], e[-1]) < ctx.resolution:
        flat = "barrierless" if path.exact else "unavailable"
        return replace(state, screen=v.model_copy(update={"verdict": flat,
                                                          "reasons": ("well_is_endpoint",)}))
    source = "screen_hei" if v is not None and v.source == "screen" else "path_hei"
    seed = _peak_seed(ctx, f"{name}_hei", source)
    return replace(state, seeds=(*state.seeds, seed)) if seed else state


def validate_intermediate(ctx: Ctx, state: CaseState, decision: Decision) -> CaseState:
    """relax_to_minimum → Registry on the collapsed saddle or on the lowest well of the latest
    DFT profile."""
    path = ctx.work.path if decision.reason == "path_intermediate" else None
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
    """Screen GS path (or screen IDPP) resampled onto the DFT endpoints, else a new IDPP."""
    frames = ctx.work.initial
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
    """DFT string in chunks until its bead energies settle (at most string_chunks), classified
    as found; a single peak seeds the saddle search. NWChem freezeN moves the frozen ends:
    seams and the final path get the DFT minima back, drift in a chunk is a note."""
    rt, p, run_name = ctx.rt, ctx.policy, f"string{state.path_runs}"
    initial = _initial_path(ctx, p.string_beads, f"{run_name}_initial")
    start, end = (ctx.mol(x) for x in ctx.ends)
    last = None
    for i in range(1, p.string_chunks + 1):
        run = rt.path.find_path(start, end, rt.method, images=p.string_beads,
                                initial_path=initial, deadline=ctx.deadline)
        if isinstance(run, Failure):
            ctx.note(f"{run_name}:{run.kind.value}")
            break
        beads = ctx.frames(run.images)
        drift = max(mapped_rmsd(beads[0], ctx.ends[0]), mapped_rmsd(beads[-1], ctx.ends[1]))
        ctx.note(f"{run_name}:c{i}:end_drift:{drift:.3f}")  # monitored, never a gate
        last = align_sequential([ctx.ends[0], *beads[1:-1], ctx.ends[1]]), run.energies_hartree
        initial = ctx.path_file(f"{run_name}_c{i}", last[0])  # true minima, no rigid jumps
        if profile.energies_settled(run.energy_history, _SETTLED * ctx.resolution):
            break
    if last is None:
        return replace(state, path_runs=state.path_runs + 1)
    verdict = ctx.verdict(last[0], last[1][1:-1], "string")  # unsettled too: a bound and a seed
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

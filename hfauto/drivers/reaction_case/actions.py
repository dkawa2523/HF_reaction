"""Reaction-case actions (design §7.3). Jobs run through the capability Protocols and so the
JobStore (idempotent); an action returns the next CaseState and keeps the evidence behind it
(saddle, TS freq, claims, new basins) in ``Work``."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import TYPE_CHECKING, Literal

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
from hfauto.chemistry.vibrations import projected_frequencies
from hfauto.chemistry.xyz import XYZ, Molecule, composition_key, geometry_fingerprint
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Evidence, Failure, FileRef, Geometry, PathProfile
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


@dataclass
class Work:
    """Evidence behind the CaseState and the artifacts the case adds."""

    initial: list[np.ndarray] | None = None  # screen GS path or screen IDPP (row 16)
    path: tuple[list[np.ndarray], tuple[float, ...]] | None = None  # last DFT string
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
        return self.policy.gates.barrier_proceed_kcal / HARTREE_TO_KCAL_MOL

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

    def verdict(self, dft: Sequence[float], *, frames: Sequence[np.ndarray] | None = None,
                low: Sequence[float] | None = None, low_ends: tuple[float, float] | None = None,
                seed: Geometry | None = None) -> BarrierVerdict:
        """barrier_verdict; ``frames`` add the node spacing and the tangent mode (below_zpe)."""
        return barrier_verdict(
            dft, dft_endpoints=self.energies, low_profile=low, low_endpoints=low_ends, seed=seed,
            policy=self.policy.gates,
            max_node_spacing_A=None if frames is None else profile.max_node_spacing(frames),
            tangent_mode_cm1=None if frames is None else _tangent_mode_cm1(self, frames),
        )

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
    return BarrierVerdict(verdict="unavailable", reasons=(reason,))


def _single_points(ctx: Ctx, frames: Sequence[np.ndarray]) -> list[float] | None:
    energies = [ctx.sp(f) for f in frames]
    return None if any(e is None for e in energies) else [e for e in energies if e is not None]


def _tangent_mode_cm1(ctx: Ctx, frames: Sequence[np.ndarray]) -> float | None:
    """Frequency of the higher DFT endpoint's mode most parallel to the path there."""
    high = int(ctx.energies[1] > ctx.energies[0])
    x = ctx.raw[high]
    step = align_mapped(x, frames[-2] if high else frames[1]) - x
    freq = ctx.rt.qm.frequencies(ctx.mol(x), ctx.rt.method, deadline=ctx.deadline)
    if isinstance(freq, Failure) or freq.hessian is None:
        return None
    freqs, modes, _ = projected_frequencies(np.load(ctx.rt.resolve(freq.hessian)), ctx.symbols, x)
    return float(freqs[int(np.argmax([overlap(m, step) for m in modes]))])


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
    verdict = ctx.verdict([ctx.energies[0], e_ts, ctx.energies[1]], seed=geometry)
    return verdict, Seed(geometry, "discovery_ts", ctx.chord(x))


def _screen_idpp(ctx: Ctx) -> BarrierVerdict:
    """Step 2 when xTB collapses the endpoints: DFT single points along the DFT IDPP."""
    try:
        frames = idpp(ctx.symbols, *ctx.ends, ctx.policy.screen_images)
    except ValueError as exc:
        return _unavailable(f"idpp:{exc}")
    ctx.work.initial = frames  # FIND_PATH (row 16) starts from it; no seed (§7.3)
    dft = _single_points(ctx, frames)
    if dft is None:
        return _unavailable("screen_single_point")
    return ctx.verdict(dft, frames=frames)


def _screen_seed(ctx: Ctx, prof: PathProfile, frames: list[np.ndarray], dft: list[float]
                 ) -> tuple[Seed, list[float], list[float]]:
    """(seed, DFT profile, low profile): the GS's TSOpt structure when its xTB freq confirms
    it (its energies follow the DFT HEI in both profiles), else the DFT HEI node."""
    low, k = list(prof.energies_hartree), 1 + int(np.argmax(dft[1:-1]))
    ts = prof.ts if prof.ts is not None and _xtb_saddle_ok(ctx, ctx.coords(prof.ts)) else None
    e_ts = None if ts is None else ctx.sp(ctx.coords(ts))
    if ts is None or e_ts is None:
        tangent = tuple(profile.tangent(frames, k).ravel())
        return Seed(ctx.geometry("screen_hei", frames[k]), "screen_hei", tangent), dft, low
    if prof.ts_energy_hartree is not None:
        low.insert(k + 1, prof.ts_energy_hartree)
    return Seed(ts, "screen_ts", ctx.chord(ctx.coords(ts))), [*dft[:k + 1], e_ts, *dft[k + 1:]], low


def _screen_path(ctx: Ctx) -> tuple[BarrierVerdict, Seed | None]:
    """Steps 2-4: xTB endpoints, then GS with TSOpt, DFT single points on every node."""
    qm, path, method = ctx.rt.screen_qm, ctx.rt.screen_path, ctx.rt.screen_method
    if qm is None or path is None or method is None:
        return _unavailable("no_screen_engines"), None
    ea, eb = (qm.optimize(ctx.mol(x), method, deadline=ctx.deadline) for x in ctx.ends)
    if isinstance(ea, Failure) or isinstance(eb, Failure):
        kind = next(o.kind.value for o in (ea, eb) if isinstance(o, Failure))
        return _unavailable(f"screen_endpoint:{kind}"), None
    low_a, low_b = ctx.coords(ea.final), ctx.coords(eb.final)
    if mapped_rmsd(low_a, low_b) <= _COLLAPSE_A:
        return _screen_idpp(ctx), None
    prof = path.find_path(ctx.mol(low_a), ctx.mol(align_mapped(low_a, low_b)), method,
                          images=ctx.policy.screen_images, refine_ts=True, deadline=ctx.deadline)
    if isinstance(prof, Failure):
        return _unavailable(f"screen_path:{prof.kind.value}"), None
    frames = ctx.work.initial = align_sequential(ctx.frames(prof.images))  # dlc images rotate
    dft = _single_points(ctx, frames)
    if dft is None:
        return _unavailable("screen_single_point"), None
    seed, dft, low = _screen_seed(ctx, prof, frames, dft)
    low_ends = (ea.energy_hartree, eb.energy_hartree)
    return ctx.verdict(dft, frames=frames, low=low, low_ends=low_ends, seed=seed.geometry), seed


def screen(ctx: Ctx, state: CaseState, decision: Decision) -> CaseState:
    """Barrier pre-check: shortcut from a low-level TS, else the low-level path (chem 16)."""
    found = None if ctx.case.low_level_ts is None else _shortcut(ctx, ctx.case.low_level_ts)
    if found is None or found[0].verdict != "proceed":
        found = _screen_path(ctx)
    verdict, seed = found
    ctx.note(f"screen:{verdict.verdict}:{','.join(verdict.reasons)}")
    seeds = (seed, *state.seeds) if seed and verdict.verdict == "proceed" else state.seeds
    return replace(state, screen=verdict, seeds=seeds)


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
    basin = ctx.rt.registry.find(ctx.symbols, x, side.energy_hartree)
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


def validate_intermediate(ctx: Ctx, state: CaseState, decision: Decision) -> CaseState:
    """relax_to_minimum → Registry on the collapsed saddle or on the lowest node between the
    maxima of a multi_max DFT path (then its highest maximum seeds when it is no intermediate)."""
    work, peaks = ctx.work, ()
    frames, energies = work.path if work.path is not None else ([], ())
    if decision.reason == "path_multi_max" and work.path is not None:
        peaks = profile.interior_maxima(energies, ctx.resolution)
        x = frames[min(range(peaks[0], peaks[-1] + 1), key=lambda i: energies[i])]
    elif work.saddle is not None:
        x = ctx.coords(work.saddle.final)
    else:
        return replace(state, intermediate="same_as_endpoint")
    name = f"int{state.saddle_attempts}_{len(state.path_runs)}"
    record = _register(ctx, x, name, "intermediate")
    ends = {ctx.rt.minima[m][0].basin_id for m in ctx.case.minima}
    if record is not None and record.basin_id not in ends:
        work.intermediate = (record, ctx.rt.species[record.species_id])
        return replace(state, intermediate="distinct")
    seeds = state.seeds
    if peaks:  # multi_max whose valley is an endpoint: try the highest maximum
        top = max(peaks, key=lambda i: energies[i])
        seed = Seed(ctx.geometry(f"{name}_hei", frames[top]), "path_hei",
                    tuple(profile.tangent(frames, top).ravel()))
        seeds = (*seeds, seed)
    return replace(state, intermediate="same_as_endpoint", seeds=seeds)


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
    inner = [align_mapped(ctx.ends[0], np.asarray(i.coords, dtype=float)) for i in images[1:-1]]
    return ctx.path_file(name, [ctx.ends[0], *inner, ctx.ends[1]])


def find_path(ctx: Ctx, state: CaseState, decision: Decision) -> CaseState:
    """DFT string in chunks until its bead energies settle (at most string_chunks); the shape
    is recorded as found and only a single_max profile yields a seed. NWChem freezeN moves the
    frozen ends: seams and the final path get the DFT minima back, drift in a chunk is a note."""
    rt, p, run_name = ctx.rt, ctx.policy, f"string{len(state.path_runs)}"
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
        return case_state.record_path(state, "failed")
    frames, energies = ctx.work.path = last
    shape = profile.shape(energies, ctx.resolution)  # a path maximum bounds the saddle
    if shape != "single_max":
        return case_state.record_path(state, shape)
    k, _, x = profile.hei(frames, energies)  # unsettled too: the HEI is only a seed
    seed = Seed(ctx.geometry(f"{run_name}_hei", x), "path_hei",
                tuple(profile.tangent(frames, round(k)).ravel()))
    return case_state.record_path(state, shape, seed)


Handler = Callable[[Ctx, CaseState, Decision], CaseState]
HANDLERS: dict[Action, Handler] = {
    Action.SCREEN: screen,
    Action.REFINE_SADDLE: refine_saddle,
    Action.VALIDATE_TS: validate_ts,
    Action.FIND_PATH: find_path,
    Action.CONNECT: connect,
    Action.VALIDATE_INTERMEDIATE: validate_intermediate,
}

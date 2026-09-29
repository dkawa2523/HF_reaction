"""Reaction-case context and the saddle, TS, connection and intermediate actions (design §7.3;
the path actions are in ``paths``). Jobs run through the capability Protocols and so the
JobStore (idempotent); an action returns the next CaseState and keeps the evidence behind it
(saddle, TS freq, claims, new basins) in ``Work``."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import TYPE_CHECKING, Literal, NamedTuple

import numpy as np

from hfauto.chemistry import profile, topology
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
from hfauto.chemistry.identity import IMAGE_A, carry, mapped_equivalent, periodic_nearest
from hfauto.chemistry.interpolation import align_mapped
from hfauto.chemistry.modes import BOUNDS_A, TARGET_HARTREE, amplitude, displace, overlap
from hfauto.chemistry.xyz import (
    XYZ,
    Molecule,
    composition_key,
    read_xyz_trajectory,
    write_xyz_trajectory,
    written_geometry,
)
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
from hfauto.drivers.reaction_case.state import CaseRules, CaseState, Decision, Seed

if TYPE_CHECKING:
    from hfauto.drivers.reaction_case.driver import CaseRuntime

_MIN_OVERLAP = 0.3  # an xTB Hessian is used when a negative mode lies along the direction
QRC_RETRY_FACTOR = 2.0  # the second QRC amplitude, capped at BOUNDS_A[1]
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
    child starts from 0; only the walltime deadline is shared)."""
    rt, seed = ctx.rt, state.seeds[0]
    counted = int(seed.source != "saddle_restart")
    state = replace(state, seeds=state.seeds[1:], saddle_attempts=state.saddle_attempts + counted,
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
    pushed, _ = displace(x, np.asarray(freq.imaginary_modes[other]), _amplitude(ctx, freq, other))
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


def _register(ctx: Ctx, coords: np.ndarray, name: str,
              source: Literal["connection", "intermediate"], opt: Evidence | None = None
              ) -> MinimumRecord | None:
    """relax_to_minimum (from ``opt`` when it has converged already, not optimized again) →
    Registry: the known basin, a new one, or None (no minimum)."""
    rt = ctx.rt
    out = relax_to_minimum(ctx.mol(coords), rt.method, rt.qm, known=rt.registry, opt=opt,
                           deadline=ctx.deadline, gates=ctx.rules.gates)
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
    """Registry.find of the side at ``x`` (its own structure, or its image's); a torsional case
    falls back to the nearest declared dihedral (CH-04); else the converged side is registered."""
    basin = ctx.rt.registry.find(side, coords=x)
    if basin is not None:
        return ctx.record(basin).minimum_id
    terms = ctx.case.coordinate
    if ctx.case.torsional and terms and all(t.kind == "dihedral" for t in terms):
        values = [declared_coordinate(terms, end) for end in ctx.raw]
        return ctx.case.minima[periodic_nearest(declared_coordinate(terms, x), values)]
    record = _register(ctx, x, name, "connection", opt=side)
    return None if record is None else record.minimum_id


def _sides(ctx: Ctx, freq: Evidence, starts: tuple[np.ndarray, np.ndarray], attempt: int
           ) -> tuple[tuple[Evidence, Evidence], tuple[np.ndarray, np.ndarray]] | None:
    """The QRC sides optimized from the TS Hessian and their final structures; None when one
    fails. At a symmetric TS the minus start is an exact image of the plus start (identity.carry
    within IMAGE_A), so on the invariant PES the minus optimum is the plus one's image: only the
    plus side runs and stands for both, the minus structure carried by that image. Two sides
    run at once (CaseRuntime.map)."""
    rt = ctx.rt
    image = carry(ctx.symbols, starts[1], starts[0], starts[0])[0] <= IMAGE_A
    if image:
        ctx.note(f"qrc{attempt}:minus_is_image")
    runs = rt.map(lambda y: rt.qm.optimize(ctx.mol(y), rt.method, init_hessian=freq,
                                           deadline=ctx.deadline), starts[:1 if image else 2])
    plus, minus = runs[0], runs[-1]
    if isinstance(plus, Failure) or isinstance(minus, Failure):
        ctx.note(f"qrc{attempt}:side_failed")
        return None
    x = ctx.coords(plus.final)
    y = carry(ctx.symbols, starts[1], starts[0], x)[1] if image else ctx.coords(minus.final)
    return (plus, minus), (x, y)


def connect(ctx: Ctx, state: CaseState, decision: Decision) -> CaseState:
    """QRC: displace ± along the TS mode by an energy target, optimize, assign, gate."""
    freq = ctx.work.ts_freq
    attempt = state.connection_attempts + 1
    state = replace(state, connection_attempts=attempt)
    if freq is None or not freq.imaginary_modes:
        return replace(state, connection="failed")
    step = min(_amplitude(ctx, freq, 0) * QRC_RETRY_FACTOR ** (attempt - 1), BOUNDS_A[1])
    optimized = _sides(ctx, freq, displace(ctx.coords(freq.final),
                                           np.asarray(freq.imaginary_modes[0]), step), attempt)
    if optimized is None:
        return replace(state, connection="failed")
    sides, finals = optimized
    first, second = (_assign(ctx, s, x, f"qrc{attempt}_{i}")
                     for i, (s, x) in enumerate(zip(sides, finals, strict=True)))
    distinct = mapped_equivalent(ctx.symbols, *finals) if ctx.case.degenerate else True
    bonds = tuple(topology.bonds(ctx.symbols, x) for x in (*ctx.ends, finals[1], finals[0]))
    gate, label = connection(freq, sides, (first, second), frozenset(ctx.case.minima),
                             degenerate=ctx.case.degenerate, sides_distinct=distinct,
                             bond_sets=bonds)
    if label == "failed" or first is None or second is None:
        ctx.note(f"qrc{attempt}:{','.join(gate.reasons)}")
        retry = "sides_same_basin" in gate.reasons  # the only failure worth a wider displacement
        return replace(state, connection="same_basin" if retry else "failed")
    side_calcs = (calc_id(ctx.keep(sides[0])), calc_id(ctx.keep(sides[1])))
    if (two := _two_steps(ctx, state, label, (first, second))) is not None:
        return two
    ctx.work.connection = ConnectionClaim(side_calcs=side_calcs, minima=(first, second))
    return replace(state, connection=label)


def _two_steps(ctx: Ctx, state: CaseState, label: str, assigned: tuple[str, str]
               ) -> CaseState | None:
    """GEN-05: a reassigned TS that joins exactly one endpoint's basin to another DFT basin
    makes the case two steps (row 5); the split child between those two validates this TS (a
    degenerate case: child 1). None for any other connection."""
    ends = [ctx.rt.minima[m][0].basin_id for m in ctx.case.minima]
    records = [ctx.rt.minima[m][0] for m in assigned]
    inside = [ends.index(r.basin_id) for r in records if r.basin_id in ends]
    if label != "reassigned" or len(inside) != 1 or state.claim is None:
        return None
    well = next(r for r in records if r.basin_id not in ends)
    ctx.note(f"qrc{state.connection_attempts}:end{inside[0]}_to_new_basin:{well.minimum_id}")
    ctx.work.intermediate = (well, ctx.rt.species[well.species_id])
    ctx.work.split_ts = (inside[0] + 1, state.claim.saddle_calc)
    return replace(state, connection=None, claim=None, intermediate="distinct")


def _with_peak(ctx: Ctx, state: CaseState, name: str) -> CaseState:
    """The latest profile's highest peak joins the seeds."""
    v = state.screen
    source = "screen_hei" if v is not None and v.source == "screen" else "path_hei"
    seed = peak_seed(ctx, f"{name}_hei", source)
    return replace(state, seeds=(*state.seeds, seed)) if seed else state


def _past_the_well(ctx: Ctx, state: CaseState, path: Profile, name: str) -> CaseState:
    """The well relaxed into an endpoint: barrierless when no interior point rises a resolution
    above the higher end, else the highest peak seeds the saddle search."""
    v, e = state.screen, path.energies
    if v is not None and max(e[1:-1]) - max(e[0], e[-1]) < ctx.resolution:
        return replace(state, screen=v.model_copy(update={"verdict": "barrierless",
                                                          "reasons": ("well_is_endpoint",)}))
    return _with_peak(ctx, state, name)


def validate_intermediate(ctx: Ctx, state: CaseState, decision: Decision) -> CaseState:
    """relax_to_minimum → Registry on the collapsed saddle (or soft TS whose QRC failed: its
    claim is withdrawn) or on the latest profile's lowest well; consumes its trigger. A failed
    relaxation is no chemical result: a well's profile goes on from its highest peak."""
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
    if record is None:
        ctx.note("int:relax_failed")
        state = replace(state, intermediate="relax_failed")
        return state if path is None else _with_peak(ctx, state, name)
    ends = {ctx.rt.minima[m][0].basin_id for m in ctx.case.minima}
    if record.basin_id not in ends:
        ctx.work.intermediate = (record, ctx.rt.species[record.species_id])
        return replace(state, intermediate="distinct")
    state = replace(state, intermediate="same_as_endpoint")
    return state if path is None else _past_the_well(ctx, state, path, name)

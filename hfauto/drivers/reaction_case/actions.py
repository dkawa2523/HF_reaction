"""Reaction-case context, the saddle search and the TS checks (design §7.3; the path actions
and their seeds are in ``paths``, VALIDATE_AND_CONNECT and the intermediate action in
``connection``). Jobs run through the capability Protocols and so the JobStore (idempotent); an
action returns the next CaseState and keeps the evidence behind it (saddle, TS freq, claims, new
basins) in ``Work``."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from hfauto.chemistry import topology
from hfauto.chemistry.gates import (
    Change,
    is_first_order_saddle,
    reaction_mode_character,
    reaction_mode_chi,
    same_spin_state,
    spin_ok,
)
from hfauto.chemistry.geometry import declared_coordinate_gradient, most_changed_dihedral
from hfauto.chemistry.identity import carry, mapped_rmsd
from hfauto.chemistry.interpolation import align_mapped
from hfauto.chemistry.modes import off_saddle
from hfauto.chemistry.profile import Point, Profile, Sample
from hfauto.chemistry.xyz import (
    XYZ,
    Molecule,
    read_xyz_trajectory,
    write_xyz_trajectory,
    written_geometry,
)
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Evidence, Failure, FileRef, Geometry
from hfauto.core.records import (
    ConnectionClaim,
    CoordinateTerm,
    MinimumRecord,
    ReactionRecord,
    SaddleClaim,
    SpeciesRecord,
)
from hfauto.drivers.minimum import calc_id, newton_push
from hfauto.drivers.reaction_case.state import CaseRules, CaseState, Decision, Seed

if TYPE_CHECKING:
    from hfauto.drivers.reaction_case.driver import CaseRuntime

MAX_DEPTH = 2  # continuations per fresh seed: a push whose search stalls still restarts


@dataclass
class Work:
    """Evidence behind the CaseState and the artifacts the case adds."""

    path: Profile | None = None  # the latest DFT profile (CaseState.screen), FIND_PATH's start
    sample: Sample | None = None  # its energy function at new nodes; None: Ctx.points
    saddle: Evidence | None = None
    depth: int = 0  # the saddle's continuation depth (Seed.depth of the search that found it)
    connection: ConnectionClaim | None = None
    intermediate: tuple[MinimumRecord, SpeciesRecord] | None = None  # basin, structure reached
    split_ts: tuple[int, str] | None = None  # (split child, saddle calc) it validates
    calcs: dict[str, Evidence] = field(default_factory=dict)
    species: dict[str, SpeciesRecord] = field(default_factory=dict)
    minima: dict[str, MinimumRecord] = field(default_factory=dict)


@dataclass
class Ctx:
    case: ReactionRecord
    rt: CaseRuntime
    rules: CaseRules
    folder: Path  # cases/<reaction_id>
    symbols: list[str]
    charge: int
    multiplicity: int
    raw: tuple[np.ndarray, np.ndarray]  # DFT minima in their own frames (JobStore reuse)
    ends: tuple[np.ndarray, np.ndarray]  # the same with b aligned onto a (path endpoints)
    energies: tuple[float, float]  # DFT minimum energies (an association: the monomers' sum)
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

    def sps(self, frames: Sequence[np.ndarray]) -> list[Evidence | None]:
        """DFT SPs at ``frames``, independent jobs run at once (CaseRuntime.map), kept and
        noted in order; per frame the lowest of its guesses' (``_guesses``)."""
        rt = self.rt
        jobs = [(i, *job) for i, x in enumerate(frames) for job in self._guesses(x)]
        found: list[Evidence | None] = [None] * len(frames)
        for (i, _, _), ev in zip(jobs, rt.map(lambda job: rt.qm.energy(
                self.mol(job[2]), rt.method, scf_guess=job[1]), jobs),
                                  strict=True):
            if isinstance(ev, Failure):
                self.note(f"sp:{ev.kind.value}")
                continue
            self.keep(ev)
            if (best := found[i]) is None or ev.energy_hartree < best.energy_hartree:
                found[i] = ev
        return found

    def points(self, new: list[tuple[int, np.ndarray]]) -> list[Point | None]:
        """The energy function (profile.Sample) of SCREEN's and a string's profile: ``sps``."""
        return [None if ev is None else (ev.energy_hartree, ev.s2)
                for ev in self.sps([x for _, x in new])]

    def _guesses(self, x: np.ndarray) -> list[tuple[Evidence | None, np.ndarray]]:
        """(scf_guess, structure) of each SP at ``x``: a closed shell's from scratch; an open
        shell's from the opt of the nearer DFT minimum, or of both when they differ in spin
        coupling (gates.same_spin_state: the broken-symmetry CH3···O2 complex 1.71, its adduct
        0.754). From scratch the middle frames of CH3· + O2 fell onto another branch up to 28.9
        kcal/mol high; an end's vectors passed directly reproduce the frame-by-frame
        continuation within 5e-8 Eh, and the lower of the two ends' is continuous. NWChem reads
        the vectors as stored, so x, aligned onto that end, goes into the opt's atom order and
        frame as the end is carried onto it (identity.carry; ends[1] is rotated onto ends[0]);
        energy and ⟨S²⟩ do not change."""
        if self.multiplicity == 1:
            return [(None, x)]
        opts, ends = self._end_opts(), self.ends
        near = int(mapped_rmsd(x, ends[1]) < mapped_rmsd(x, ends[0]))
        apart = not same_spin_state(*opts, self.rules.gates)
        return [(opts[k], carry(self.symbols, self.coords(opts[k].final), ends[k],
                                align_mapped(ends[k], x))[1])
                for k in ((0, 1) if apart else (near,))]

    def end_records(self) -> tuple[MinimumRecord, MinimumRecord]:
        """The current Registry records of the case's DFT minima."""
        minima = self.rt.registry.minima
        return minima[self.case.minima[0]][0], minima[self.case.minima[1]][0]

    def _end_opts(self) -> tuple[Evidence, Evidence]:
        a, b = self.end_records()
        return self.rt.calcs[a.opt_calc], self.rt.calcs[b.opt_calc]

    def end_s2(self) -> tuple[float | None, float | None]:
        """The DFT minima's ⟨S²⟩ (None: a closed shell)."""
        if self.multiplicity == 1:
            return None, None
        a, b = self._end_opts()
        return a.s2, b.s2

    def change(self) -> Change:
        """(formed, broken) between the labelled case ends: the one bond change that ρ, χ, the
        retry's mode and the degenerate QRC check share. A TS lent with other atom labels is
        measured on it too (its lender's ends are not on the record)."""
        return topology.bond_changes(self.symbols, *self.ends)

    def geometry(self, name: str, coords: np.ndarray) -> Geometry:
        return written_geometry(self.mol(coords).write(self.folder / f"{name}.xyz"),
                                self.rt.file_ref)

    def path_file(self, name: str, frames: Sequence[np.ndarray]) -> FileRef:
        images, path = [self.mol(f).xyz for f in frames], self.folder / f"{name}.xyz"
        return self.rt.file_ref(write_xyz_trajectory(images, path))

    def frames(self, ref: FileRef) -> list[np.ndarray]:
        images = read_xyz_trajectory(self.rt.resolve(ref))
        return [np.asarray(i.coords, dtype=float) for i in images]

    def direction(self, x: np.ndarray, mode: Sequence[float] | None = None
                  ) -> tuple[str, np.ndarray]:
        """(kind, 3N vector) of the reaction direction at ``x`` (design §7.3): a TS seed's own
        imaginary ``mode``; else ρ = ∇(Σ_broken r − Σ_formed r) of the case's bond change; with
        none, the gradient of the declared coordinate, else of the most changed dihedral, else
        the chord between the ends aligned onto x."""
        if mode is not None:
            return "mode", np.asarray(mode)
        formed, broken = self.change()
        if formed or broken:
            rho = [CoordinateTerm(kind="distance", atoms=bond,
                                  coefficient=1.0 if bond in broken else -1.0)
                   for bond in sorted(formed | broken)]
            return "rho", declared_coordinate_gradient(rho, x)
        a, b = (align_mapped(x, end) for end in self.ends)
        if terms := self.case.coordinate or most_changed_dihedral(
                topology.bonds(self.symbols, a), a, b):
            return "coordinate", declared_coordinate_gradient(terms, x)
        return "chord", (b - a).ravel()


def xtb_freq(ctx: Ctx, coords: np.ndarray) -> Evidence | None:
    """The xTB freq at ``coords``; None without a low-level engine or when it fails."""
    rt = ctx.rt
    if rt.screen_qm is None or rt.screen_method is None:
        return None
    freq = rt.screen_qm.frequencies(ctx.mol(coords), rt.screen_method)
    return None if isinstance(freq, Failure) else freq


def _seed_hessian(ctx: Ctx, seed: Seed, x: np.ndarray) -> tuple[str, Evidence | Failure]:
    """The seed's own TS freq, else the xTB freq at the seed whenever there is one (the model
    keeps only its curvatures off the direction), else the DFT freq."""
    if seed.hessian is not None:
        return "ts_freq", seed.hessian
    if (xtb := xtb_freq(ctx, x)) is not None:
        return "xtb", xtb
    return "dft", ctx.rt.qm.frequencies(ctx.mol(x), ctx.rt.method)


def _search(ctx: Ctx, seed: Seed) -> Evidence | Failure:
    """saddle.refine from ``seed`` with only its reaction direction negative in its initial
    Hessian; a push starts its SCF from its TS freq."""
    rt, x = ctx.rt, ctx.coords(seed.geometry)
    kind, direction = ctx.direction(x, seed.mode)
    source, hessian = _seed_hessian(ctx, seed, x)
    if isinstance(hessian, Failure):
        ctx.note(f"saddle_hessian:{hessian.kind.value}")
        return hessian
    ctx.note(f"saddle_hessian:{source}:{kind}")
    result = rt.saddle.refine(ctx.mol(x), rt.method, hessian=hessian, mode=tuple(direction),
                              scf_guess=seed.hessian)
    if isinstance(result, Failure):
        ctx.note(f"saddle:{result.kind.value}:{result.reason}")
    return result


def _restart(ctx: Ctx, seed: Seed, stalled: Failure) -> Seed | None:
    """A search stalled at maxiter goes on from its last frame with a fresh Hessian, one
    continuation deeper (up to MAX_DEPTH), unless that frame lies more than a resolution above
    the latest DFT profile's maximum: a continuous path bounds the saddle from above, so the
    search has climbed past the barrier it was to find. A stalled search stored without its
    energy is not bounded."""
    path, energy = ctx.work.path, stalled.energy_hartree
    if stalled.final is None or seed.depth >= MAX_DEPTH:
        return None
    if energy is not None and path is not None and energy > max(path.energies) + ctx.resolution:
        ctx.note(f"saddle:above_path_bound:{path.source}")
        return None
    return replace(seed, geometry=stalled.final, mode=None, hessian=None, depth=seed.depth + 1)


def refine_saddle(ctx: Ctx, state: CaseState, decision: Decision) -> CaseState:
    """The front seed → a saddle search, restarted once within the attempt when it stalls
    (``_restart``). Every seed counts one attempt (per case: a split child starts from 0). A
    new search drops the claim and connection verdict of an earlier TS whose sides joined one
    basin or state (row 7)."""
    seed = state.seeds[0]
    state = replace(state, seeds=state.seeds[1:], saddle_attempts=state.saddle_attempts + 1,
                    last_saddle="failed", claim=None, connection=None)
    result = _search(ctx, seed)
    if isinstance(result, Failure) and (restart := _restart(ctx, seed, result)) is not None:
        seed, result = restart, _search(ctx, restart)
    if isinstance(result, Failure):
        return state
    ctx.work.saddle, ctx.work.depth = result, seed.depth
    return replace(state, last_saddle="converged")


def _newton_start(ctx: Ctx, saddle: Evidence, freq: Evidence, x: np.ndarray) -> np.ndarray | None:
    """The signed Newton step (minimum.newton_push) from a saddle that is not stationary by a
    minimum's criterion, when a freq there has at most one mode below -saddle_cm1: its second
    imaginary mode came from the missing stationarity (OH + CH4: ν2 −54.8i and −77.1i
    turned real one step away), noted higher_order:not_stationary. None otherwise: a saddle
    with no gradient, a stationary one, or a second mode that persists at the step."""
    rt = ctx.rt
    if saddle.gradient is None or freq.hessian is None:
        return None
    newton = newton_push(np.load(rt.resolve(freq.hessian)), saddle.gradient, ctx.mol(x).xyz,
                         signed=True)
    if newton is None:
        return None
    start = np.reshape(x, (-1, 3)) + newton
    stepped = rt.qm.frequencies(ctx.mol(start), rt.method, scf_guess=saddle)
    saddle_cm1 = ctx.rules.gates.saddle_cm1
    if isinstance(stepped, Failure) or sum(
            f < -saddle_cm1 for f in stepped.frequencies_cm1 or ()) >= 2:
        return None
    ctx.note("higher_order:not_stationary")
    return start


def _retry(ctx: Ctx, saddle: Evidence, freq: Evidence, x: np.ndarray, name: str) -> Seed:
    """The seed after a higher-order verdict, one continuation deeper than the saddle: its
    verified freq is the seed's Hessian and its reaction mode (the imaginary mode of largest χ,
    the TS gate's measure, on the case's bond change, else on its reaction direction) the
    seed's mode. It starts at the saddle's Newton step when the second imaginary mode is an
    artefact of a non-stationary point (_newton_start), else at the saddle pushed once off its
    other modes below -saddle_cm1 (modes.off_saddle)."""
    (formed, broken), (_, along), modes = ctx.change(), ctx.direction(x), freq.imaginary_modes
    r = int(np.argmax([reaction_mode_chi(ctx.symbols, m, x, formed | broken, along) or 0.0
                       for m in modes]))
    start = _newton_start(ctx, saddle, freq, x)
    if start is None:
        start = x + off_saddle(freq, ctx.symbols, below_cm1=ctx.rules.gates.saddle_cm1, keep=r)
    return Seed(ctx.geometry(name, start), "higher_order_retry", modes[r], freq,
                depth=ctx.work.depth + 1)


def validate_ts(ctx: Ctx, state: CaseState) -> CaseState:
    """Separate DFT freq on the saddle → is_first_order_saddle, reaction_mode_character (on the
    case's bond change, else the declared coordinate) and spin_ok: the claim of an accepted
    TS, its freq kept. A rejected saddle, one without an imaginary mode included, stays a
    counted attempt and gets no QRC; a higher-order one is retried (_retry) while its depth
    allows. A stationary point with -saddle_cm1 < ν < -noise_cm1 is accepted as a TS: χ and QRC
    decide whether it is a TS of this case."""
    rt, saddle, gates = ctx.rt, ctx.work.saddle, ctx.rules.gates
    state = replace(state, last_saddle="failed")
    if saddle is None:
        return state
    x = ctx.coords(saddle.final)
    freq = rt.qm.frequencies(ctx.mol(x), rt.method, scf_guess=saddle)
    if isinstance(freq, Failure):
        ctx.note(f"ts_freq:{freq.kind.value}")
        return state
    gate = is_first_order_saddle(freq, saddle=saddle, policy=gates)
    if gate:
        terms, (formed, broken) = ctx.case.coordinate, ctx.change()
        mode = reaction_mode_character(ctx.symbols, freq.imaginary_modes[0], x, formed | broken,
                                       declared_coordinate_gradient(terms, x) if terms else None)
        gate = replace(gate, ok=mode.ok, reasons=mode.reasons)
    if gate:
        freq_calc, notes = calc_id(ctx.keep(freq)), (*gate.notes, *spin_ok(freq, gates).reasons)
        claim = SaddleClaim(saddle_calc=calc_id(ctx.keep(saddle)), freq_calc=freq_calc,
                            imag_cm1=min(freq.frequencies_cm1 or (0.0,)),
                            energy_hartree=freq.energy_hartree, notes=notes)
        return replace(state, last_saddle=None, claim=claim)
    ctx.note(f"ts_rejected:{','.join(gate.reasons)}")
    if "higher_order" in gate.reasons and ctx.work.depth < MAX_DEPTH:
        seed = _retry(ctx, saddle, freq, x, f"retry{state.saddle_attempts}")
        return replace(state, seeds=(seed, *state.seeds))
    return state

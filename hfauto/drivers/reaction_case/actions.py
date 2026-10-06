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
from hfauto.chemistry.profile import Point, Profile
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
from hfauto.drivers.minimum import calc_id, not_stationary
from hfauto.drivers.reaction_case.state import CaseRules, CaseState, Seed

if TYPE_CHECKING:
    from hfauto.drivers.reaction_case.driver import CaseRuntime

# A profile's energy function at new nodes, naming the calculation of each (profile.Sample's
# Point and its calc id; None when its SP failed)
Sampler = Callable[[list[tuple[int, np.ndarray]]], list[tuple[Point, str] | None]]


@dataclass
class Work:
    """Evidence behind the CaseState and the artifacts the case adds."""

    path: Profile | None = None  # the latest DFT profile (CaseState.screen), FIND_PATH's start
    points: tuple[str, ...] = ()  # the calc id of each of its nodes ("": none of its own)
    sample: Sampler | None = None  # its energy function at new nodes; None: Ctx.points
    saddle: Evidence | None = None
    depth: int = 0  # Seed.depth of the search that found the saddle (1: continued)
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

    def points(self, new: list[tuple[int, np.ndarray]]) -> list[tuple[Point, str] | None]:
        """The energy function (``Sampler``) of SCREEN's and a string's profile: ``sps``."""
        return [None if ev is None else ((ev.energy_hartree, ev.s2), calc_id(ev))
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
    """The seed's own DFT freq (a TS seed's verified one, a continuation's), else the xTB freq at
    the seed whenever there is one (the model keeps only its curvatures off the direction), else
    the DFT freq."""
    if seed.hessian is not None:
        return "seed_freq", seed.hessian
    if (xtb := xtb_freq(ctx, x)) is not None:
        return "xtb", xtb
    return "dft", ctx.rt.qm.frequencies(ctx.mol(x), ctx.rt.method)


def _search(ctx: Ctx, seed: Seed) -> Evidence | Failure:
    """saddle.refine from ``seed`` with only its reaction direction negative in its initial
    Hessian; a seed with its own freq starts its SCF from it."""
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


def continuation(ctx: Ctx, x: np.ndarray, freq: Evidence, name: str, *, push: bool) -> Seed:
    """The one continuation of a search (design §7.3), once from a fresh seed (depth 1): the
    DFT ``freq`` at ``x`` is its Hessian and SCF guess, shaped along the imaginary mode of
    largest χ (the TS gate's measure on the case's bond change, else on its reaction direction;
    with none, the reaction direction). It starts at ``x``, whose gradient leads the search;
    ``push`` (a stationary higher-order saddle, whose gradient vanishes): pushed once off its
    other modes below -saddle_cm1 (modes.off_saddle)."""
    modes, r = freq.imaginary_modes, None
    if modes:
        (formed, broken), (_, along) = ctx.change(), ctx.direction(x)
        r = int(np.argmax([reaction_mode_chi(ctx.symbols, m, x, formed | broken, along) or 0.0
                           for m in modes]))
    start = x + off_saddle(freq, ctx.symbols, below_cm1=ctx.rules.gates.saddle_cm1,
                           keep=r) if push else x
    return Seed(ctx.geometry(name, start), "continuation", None if r is None else modes[r], freq,
                depth=1)


def _stalled(ctx: Ctx, seed: Seed, stalled: Failure, name: str) -> Seed | None:
    """A fresh seed's search stopped at maxiter continues from its last frame (``continuation``
    from a DFT freq there), unless that frame lies more than a resolution above the latest DFT
    profile's maximum: a continuous path bounds the saddle from above, so the search has
    climbed past the barrier it was to find. A stalled search stored without its energy is not
    bounded."""
    path, energy = ctx.work.path, stalled.energy_hartree
    if stalled.final is None or seed.depth:
        return None
    if energy is not None and path is not None and energy > max(path.energies) + ctx.resolution:
        ctx.note(f"saddle:above_path_bound:{path.source}")
        return None
    x = ctx.coords(stalled.final)
    freq = ctx.rt.qm.frequencies(ctx.mol(x), ctx.rt.method)
    if isinstance(freq, Failure):
        ctx.note(f"continuation_freq:{freq.kind.value}")
        return None
    return continuation(ctx, x, freq, name, push=False)


def refine_saddle(ctx: Ctx, state: CaseState) -> CaseState:
    """The front seed → a saddle search, continued within the attempt when it stalls
    (``_stalled``). Every seed is one attempt (per case: a split child starts from 0); a new
    search drops the claim and connection of an earlier TS."""
    seed = state.seeds[0]
    state = replace(state, seeds=state.seeds[1:], attempts=state.attempts + 1,
                    saddle_pending=False, claim=None, connection=None)
    result = _search(ctx, seed)
    if isinstance(result, Failure) and (
            cont := _stalled(ctx, seed, result, f"continuation{state.attempts}")) is not None:
        seed, result = cont, _search(ctx, cont)
    if isinstance(result, Failure):
        return state
    ctx.work.saddle, ctx.work.depth = result, seed.depth
    return replace(state, saddle_pending=True)


def validate_ts(ctx: Ctx, state: CaseState) -> CaseState:
    """Separate DFT freq on the saddle → is_first_order_saddle, reaction_mode_character (on the
    case's bond change, else the declared coordinate) and the certification
    (minimum.not_stationary) → the claim of an accepted TS, its freq kept. A higher-order saddle,
    or a first-order one that fails only the certification, is continued once from a fresh seed
    (``continuation``, the next seed); a continued one that still fails the certification is
    accepted, noted not_stationary (thermo blocks it). Any other rejection, one without an
    imaginary mode included, is a failure token without QRC. A stationary point with
    -saddle_cm1 < ν < -noise_cm1 is a TS like any other: χ and QRC decide whether it is this
    case's."""
    rt, saddle, gates = ctx.rt, ctx.work.saddle, ctx.rules.gates
    state = replace(state, saddle_pending=False)
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
    step = not_stationary(freq, saddle.gradient, x, rt.resolve)
    fresh = not ctx.work.depth
    if gate and (step is None or not fresh):
        freq_calc = calc_id(ctx.keep(freq))
        notes = (*gate.notes, *spin_ok(freq, gates).reasons,
                 *(() if step is None else ("not_stationary",)))
        claim = SaddleClaim(saddle_calc=calc_id(ctx.keep(saddle)), freq_calc=freq_calc,
                            imag_cm1=min(freq.frequencies_cm1 or (0.0,)),
                            energy_hartree=freq.energy_hartree, notes=notes)
        return replace(state, claim=claim)
    ctx.note(f"ts_rejected:{','.join(gate.reasons or ('not_stationary',))}")
    higher = "higher_order" in gate.reasons
    if fresh and (higher or gate):
        seed = continuation(ctx, x, freq, f"continuation{state.attempts}",
                            push=higher and step is None)
        return replace(state, seeds=(seed, *state.seeds))
    return state

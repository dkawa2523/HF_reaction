"""Reaction-case path actions (design §7.3): SCREEN's shortcut from the low-level TSs, its
low-level barrier pre-check or an association's relaxed scan, and FIND_PATH's DFT string chunk.
Each DFT profile between the ends (an association: from its separated monomers, spin-projected
for a low-spin-coupled pair) is judged by ``profile.judge`` and may seed the saddle search at
its peak; its verdict names the calculation of each node (``BarrierVerdict.points``)."""

from __future__ import annotations

from collections.abc import Collection, Sequence
from dataclasses import replace

import numpy as np

from hfauto.chemistry import profile, topology
from hfauto.chemistry.electronic_state import low_spin_coupled
from hfauto.chemistry.gates import spin_ok
from hfauto.chemistry.interpolation import align_mapped, align_sequential, idpp, resample
from hfauto.chemistry.profile import Profile
from hfauto.core.evidence import Evidence, Failure, FileRef, Geometry
from hfauto.core.records import BarrierVerdict, ProfileSource
from hfauto.drivers.minimum import calc_id
from hfauto.drivers.reaction_case import state as case_state
from hfauto.drivers.reaction_case.actions import Ctx, Sampler, xtb_freq
from hfauto.drivers.reaction_case.state import CaseState, Seed

SCREEN_IMAGES = 11  # SCREEN's NEB images including both ends
STRING_BEADS = 9  # DFT string beads including both ends
# An association's scan: constrained points SCAN_STEP_A apart along the formed bond, the first
# SCAN_POINTS - 1 of them inwards to the adduct (the 8th node), at most SCAN_MAX_POINTS in all
# (5.1 A beyond the adduct's bond)
SCAN_STEP_A = 1.5 / 7
SCAN_POINTS = 8
SCAN_MAX_POINTS = 24


def peak_seeds(ctx: Ctx, name: str, path: Profile | None = None) -> tuple[Seed, ...]:
    """The highest peak detected on ``path`` (default: the latest profile), refined by a
    parabola (none without): a screen_hei seed on SCREEN's profile, else a path_hei one."""
    path = path or ctx.work.path
    peaks = () if path is None else profile.interior_maxima(path.energies, ctx.resolution)
    if path is None or not peaks:
        return ()
    _, _, x = profile.hei(path.frames, path.energies, max(peaks, key=path.energies.__getitem__))
    return (Seed(ctx.geometry(name, x), "screen_hei" if path.source == "screen" else "path_hei"),)


def judge(ctx: Ctx, path: Profile, ids: Sequence[str], sample: Sampler
          ) -> tuple[BarrierVerdict, Profile, tuple[str, ...]]:
    """``profile.judge`` of ``path``, whose nodes are the calculations ``ids`` ("": a node of no
    calculation of its own, as the separated monomers' sum or a string bead), with ``sample``
    at its new nodes: the verdict naming the calculation of each node it was judged on
    (``BarrierVerdict.points``), that profile and those ids."""
    named: list[str] = []

    def energies(new: list[tuple[int, np.ndarray]]) -> list[profile.Point | None]:
        found = sample(new)
        named.extend(p[1] for p in found if p is not None)
        return [None if p is None else p[0] for p in found]

    verdict, judged = profile.judge(path, ctx.rules.gates, energies)
    ids, points, k = (*ids, *[""] * (len(path.frames) - len(ids))), [], 0
    for frame in judged.frames:  # the old nodes in order, a new one (a midpoint) between them
        old = k < len(path.frames) and frame is path.frames[k]
        points.append(ids[k] if old else named.pop(0))
        k += old
    return verdict.model_copy(update={"points": tuple(points)}), judged, tuple(points)


def judged(ctx: Ctx, path: Profile, ids: Sequence[str], sample: Sampler) -> BarrierVerdict:
    """``judge`` of ``path``; the profile it was judged on (with its new nodes, noted), their
    calculations and its energy function become the latest ones."""
    verdict, ctx.work.path, ctx.work.points = judge(ctx, path, ids, sample)
    ctx.work.sample = sample
    if len(ctx.work.path.energies) > len(path.energies):
        ctx.note(f"{path.source}_midpoints:{verdict.verdict}")
    return verdict


def _ends(ctx: Ctx) -> tuple[str, str]:
    """The calc ids of the DFT minima at a profile's ends (their opts)."""
    a, b = ctx.end_records()
    return a.opt_calc, b.opt_calc


def _unavailable(reason: str, source: ProfileSource = "screen") -> BarrierVerdict:
    return BarrierVerdict(verdict="unavailable", source=source, reasons=(reason,))


def _xtb_ts_mode(ctx: Ctx, coords: np.ndarray) -> tuple[float, ...] | None:
    """A low-level TS's imaginary mode: its xTB freq has exactly one mode below -saddle_cm1."""
    freq, below = xtb_freq(ctx, coords), -ctx.rules.gates.saddle_cm1
    if freq is None or sum(nu < below for nu in freq.frequencies_cm1 or ()) != 1:
        return None
    return tuple(freq.imaginary_modes[0])  # ordered by frequency


def _shortcut(ctx: Ctx) -> tuple[Seed, ...]:
    """Step 1: the low-level TSs already known for this pair of states (discovery /
    mode-follow), in order. Each whose xTB freq confirms its mode gets a DFT SP; those lying at
    least a resolution above both DFT minima seed the saddle search in that order. The three
    points only pick the seeds: they are no profile and give no verdict."""
    seeds = []
    for geometry in ctx.case.low_level_ts:
        if (mode := _xtb_ts_mode(ctx, ctx.coords(geometry))) is None:
            ctx.note("low_level_ts_rejected")
        else:
            seeds.append(Seed(geometry, "discovery_ts", mode))
    found = ctx.sps([ctx.coords(s.geometry) for s in seeds]) if seeds else []
    top = max(ctx.energies) + ctx.resolution
    return tuple(seed for seed, ev in zip(seeds, found, strict=True)
                 if ev is not None and ev.energy_hartree >= top)


def _neb(ctx: Ctx, frames: list[np.ndarray]) -> tuple[list[np.ndarray], Geometry | None]:
    """The xTB CI-NEB from ``frames`` (fixed ends) and its TS; ``frames`` without a NEB."""
    engine, method = ctx.rt.screen_path, ctx.rt.screen_method
    if engine is None or method is None:
        return frames, None
    neb = engine.find_path(ctx.mol(frames[0]), ctx.mol(frames[-1]), method, images=len(frames),
                           initial_path=ctx.path_file("screen_idpp", frames))
    if isinstance(neb, Failure):
        ctx.note(f"screen_neb:{neb.kind.value}")
        return frames, None
    return align_sequential(ctx.frames(neb.images)), neb.ts


def _screen_path(ctx: Ctx) -> tuple[BarrierVerdict, tuple[Seed, ...]]:
    """Steps 2-3: the NEB from the IDPP between the DFT minima, DFT SPs on its interior; a single
    peak seeds at the NEB's TS when its xTB freq confirms it, else at the DFT peak."""
    try:
        frames, ts = _neb(ctx, idpp(ctx.symbols, *ctx.ends, SCREEN_IMAGES))
    except ValueError as exc:
        return _unavailable(f"idpp:{exc}"), ()
    inner = [ev for ev in ctx.sps(frames[1:-1]) if ev is not None]
    if len(inner) < len(frames) - 2:
        return _unavailable("screen_single_point"), ()
    (first, last), e, (a, b) = ctx.end_s2(), ctx.energies, _ends(ctx)
    verdict = judged(ctx, Profile(frames, (e[0], *(ev.energy_hartree for ev in inner), e[1]),
                                  "screen", (first, *(ev.s2 for ev in inner), last)),
                     (a, *map(calc_id, inner), b), ctx.points)
    if verdict.verdict != "single":
        return verdict, ()
    mode = None if ts is None else _xtb_ts_mode(ctx, ctx.coords(ts))
    if ts is not None and mode is not None:
        return verdict, (Seed(ts, "screen_ts", mode),)
    return verdict, peak_seeds(ctx, "screen_hei")


def _moved(x: np.ndarray, bond: tuple[int, int], fragment: Collection[int], r: float
           ) -> np.ndarray:
    """``x`` with ``fragment`` (holding j) translated rigidly along i→j until r_ij = r."""
    (i, j), y = bond, np.array(x, dtype=float)
    axis = y[j] - y[i]
    y[list(fragment)] += (r / float(np.linalg.norm(axis)) - 1.0) * axis
    return y


def _scan(ctx: Ctx) -> tuple[BarrierVerdict, tuple[Seed, ...]]:
    """An association: the separated monomers → the adduct along a relaxed scan of the formed
    bond (i, j), its point k at r_P + k·SCAN_STEP_A (r_P: the adduct's). Points SCAN_POINTS - 1
    down to 1 run inwards: the first is the adduct with the fragment holding j moved out along
    i→j (its SCF from scratch: a broken-symmetry pair forms apart), each next the previous
    optimum moved in, its SCF from the previous point's (the pair stays on its continuous
    branch). The scan then reaches outwards a point at a time from its longest one (moved out,
    its SCF from it) until that point lies within a resolution of the separated monomers'
    energy sum, its neighbour on the profile: through any precursor complex, no hill or well of
    a resolution is left beyond the scan. Each point is optimized with r_ij fixed. A failed
    point, or the sum not reached within SCAN_MAX_POINTS points, leaves no profile
    (unavailable). The profile is [ΣE(monomers), scan..., E(adduct)] (the points' energies by
    ``_projected``), the longest point standing as the monomers' frame; a new node is the SP at
    the linear midpoint of two neighbours (a point of the continuous path), started from the
    SCF of the one before it, and projected alike."""
    rt, adduct, asymptote = ctx.rt, ctx.ends[1], ctx.energies[0]
    [(i, j)] = ctx.change()[0]  # the one formed bond
    fragment = next(f for f in topology.fragments(ctx.symbols, ctx.ends[0]) if j in f)
    r_p = float(np.linalg.norm(adduct[j] - adduct[i]))
    points: list[Evidence] = []  # the longest first

    def point(k: int, x: np.ndarray, guess: Evidence | None) -> Evidence | None:
        r = r_p + k * SCAN_STEP_A
        opt = rt.qm.optimize(ctx.mol(_moved(x, (i, j), fragment, r)), rt.method,
                             fixed_bond=(i, j, r), scf_guess=guess)
        if isinstance(opt, Failure):
            ctx.note(f"scan{k}:{opt.kind.value}")
            return None
        return ctx.keep(opt)

    x, guess = adduct, None
    for k in range(SCAN_POINTS - 1, 0, -1):
        if (opt := point(k, x, guess)) is None:
            return _unavailable("scan_point", "scan"), ()
        points.append(opt)
        x, guess = ctx.coords(opt.final), opt
    frames = [ctx.coords(p.final) for p in points]
    energies = _projected(ctx, frames, points, [f"scan{k}" for k in range(len(points), 0, -1)])
    k = len(points)
    while abs(energies[0] - asymptote) >= ctx.resolution:
        k += 1
        if k > SCAN_MAX_POINTS:
            return _unavailable("scan_reach", "scan"), ()
        if (opt := point(k, frames[0], points[0])) is None:
            return _unavailable("scan_point", "scan"), ()
        points.insert(0, opt)
        frames.insert(0, ctx.coords(opt.final))
        energies.insert(0, _projected(ctx, frames[:1], [opt], [f"scan{k}"])[0])
    ctx.note(f"scan:{i}-{j}:{r_p:.3f}:{len(points)}_points")

    def sample(new: list[tuple[int, np.ndarray]]) -> list[tuple[profile.Point, str] | None]:
        # node n >= 1 is scan point n - 1; its SCF is read in its own frame
        sps = rt.map(lambda n: rt.qm.energy(ctx.mol(align_mapped(frames[n[0] - 1], n[1])),
                                            rt.method, scf_guess=points[n[0] - 1]), new)
        if failed := [ev for ev in sps if isinstance(ev, Failure)]:
            ctx.note(f"scan_mid:{failed[0].kind.value}")
            return [None] * len(new)
        mids = [ctx.keep(ev) for ev in sps]
        names = [f"scan_mid{n}" for n, _ in new]
        return [((e, ev.s2), calc_id(ev)) for e, ev in zip(
            _projected(ctx, [x for _, x in new], mids, names), mids, strict=True)]

    (_, s2_adduct), (_, end) = ctx.end_s2(), _ends(ctx)
    nodes = align_sequential([*frames, adduct])
    path = Profile([nodes[0], *nodes], (asymptote, *energies, ctx.energies[1]), "scan",
                   (None, *(ev.s2 for ev in points), s2_adduct))
    verdict = judged(ctx, path, ("", *map(calc_id, points), end), sample)
    return verdict, peak_seeds(ctx, "scan_hei") if verdict.verdict == "single" else ()


def projected(e_bs: float, s2_bs: float, e_hs: float, s2_hs: float, spin: float) -> float:
    """Yamaguchi's approximate spin projection (CPL 149, 537 (1988)): the energy of the pure
    spin-``spin`` state from a broken-symmetry solution and the high-spin one at its structure,
    E_LS = E_BS + α(E_BS − E_HS), α = (⟨S²⟩_BS − S(S+1)) / (⟨S²⟩_HS − ⟨S²⟩_BS)."""
    alpha = (s2_bs - spin * (spin + 1)) / (s2_hs - s2_bs)
    return e_bs + alpha * (e_bs - e_hs)


def _projected(ctx: Ctx, frames: list[np.ndarray], points: list[Evidence],
               names: Sequence[str]) -> list[float]:
    """The energies of scan structures ``frames`` from their ``points``: those of a
    low-spin-coupled pair (CH3· + O2 as a doublet, electronic_state.low_spin_coupled),
    broken-symmetry where the monomers separate, are spin-projected (``projected``) with an SP
    of the high-spin coupling Σ(m_i − 1) + 1 at each structure. One whose high-spin SP fails or
    is itself contaminated (spin_ok) keeps its broken-symmetry energy, noted
    ``<its name>:ap_skipped``. The monomers' sum and the adduct keep theirs; spin_ok and the
    gates of stationary points and thermochemistry are untouched."""
    rt, energies = ctx.rt, [ev.energy_hartree for ev in points]
    spins = [rt.species[rt.registry.minima[m][0].species_id].multiplicity
             for m in ctx.case.monomers]
    if not low_spin_coupled(spins, ctx.multiplicity):
        return energies
    high, spin = sum(m - 1 for m in spins) + 1, (ctx.multiplicity - 1) / 2
    sps = rt.map(lambda x: rt.qm.energy(replace(ctx.mol(x), multiplicity=high), rt.method),
                 frames)
    for k, (name, bs, hs) in enumerate(zip(names, points, sps, strict=True)):
        if isinstance(hs, Failure):
            ctx.note(f"{name}:ap_skipped:{hs.kind.value}")
            continue
        ctx.keep(hs)
        if hs.s2 is None or bs.s2 is None or not spin_ok(hs, ctx.rules.gates):
            ctx.note(f"{name}:ap_skipped:s2={hs.s2}")
        else:
            energies[k] = projected(bs.energy_hartree, bs.s2, hs.energy_hartree, hs.s2, spin)
    return energies


def screen(ctx: Ctx, state: CaseState) -> CaseState:
    """Barrier pre-check (design §7.3): an association's scan; else first, once, the shortcut's
    seeds from the low-level TSs, and the low-level path when it gives none or once they have
    failed (a second SCREEN)."""
    if not ctx.case.monomers and not state.shortcut_done:
        state = replace(state, shortcut_done=True)
        if seeds := _shortcut(ctx):
            return replace(state, seeds=(*state.seeds, *seeds))
    verdict, seeds = (_scan if ctx.case.monomers else _screen_path)(ctx)
    ctx.note(f"{verdict.source}:{verdict.verdict}:{','.join(verdict.reasons)}")
    return case_state.record_profile(replace(state, neb_done=True), verdict, seeds)


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
    inner = resample(frames, beads)[1:-1]
    return ctx.path_file(name, align_sequential([ctx.ends[0], *inner, ctx.ends[1]]))


def find_path(ctx: Ctx, state: CaseState) -> CaseState:
    """One DFT string chunk, classified at once (unconverged too: its maximum bounds the saddle
    and its peak is a seed); a single peak seeds the saddle search. The path gets the DFT
    minima back at its ends, which NWChem's freezeN moves. A chunk without an answer (failed,
    or unavailable) is one failure token."""
    rt, run_name = ctx.rt, f"string{state.path_runs}"
    failed = replace(state, path_runs=state.path_runs + 1, attempts=state.attempts + 1)
    initial = _initial_path(ctx, STRING_BEADS, f"{run_name}_initial")
    if initial is None:
        return failed
    run = rt.path.find_path(ctx.mol(ctx.ends[0]), ctx.mol(ctx.ends[1]), rt.method,
                            images=STRING_BEADS, initial_path=initial)
    if isinstance(run, Failure):
        ctx.note(f"{run_name}:{run.kind.value}")
        return failed
    frames = align_sequential([ctx.ends[0], *ctx.frames(run.images)[1:-1], ctx.ends[1]])
    e, (a, b) = ctx.energies, _ends(ctx)
    verdict = judged(ctx, Profile(frames, (e[0], *run.energies_hartree[1:-1], e[1]), "string"),
                     (a, *[""] * (len(frames) - 2), b), ctx.points)
    seeds = peak_seeds(ctx, f"{run_name}_hei") if verdict.verdict == "single" else ()
    state = case_state.record_profile(state, verdict, seeds)
    return replace(state, attempts=state.attempts + int(verdict.verdict == "unavailable"))

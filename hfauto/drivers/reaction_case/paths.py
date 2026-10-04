"""Reaction-case path actions (design §7.3): SCREEN's low-level barrier pre-check, or an
association's relaxed scan, and FIND_PATH's DFT string chunk. Each classifies a DFT profile
between the DFT minima (an association: from its separated monomers, spin-projected for a
low-spin-coupled pair) and may seed the saddle search at its peak."""

from __future__ import annotations

from collections.abc import Collection
from dataclasses import replace

import numpy as np

from hfauto.chemistry import profile, topology
from hfauto.chemistry.electronic_state import low_spin_coupled
from hfauto.chemistry.gates import barrier_verdict, spin_ok
from hfauto.chemistry.interpolation import align_sequential, idpp, resample
from hfauto.core.evidence import Evidence, Failure, FileRef, Geometry
from hfauto.core.records import BarrierVerdict
from hfauto.drivers.reaction_case import state as case_state
from hfauto.drivers.reaction_case.actions import Ctx, Profile, xtb_freq
from hfauto.drivers.reaction_case.state import CaseState, Decision, Seed

SCREEN_IMAGES = 11  # SCREEN's NEB images including both ends
STRING_BEADS = 9  # DFT string beads including both ends
# An association's scan: its longest point lies SCAN_REACH_A beyond the adduct's
# bond, outside the Coulson-Fischer region (~2.2 A for C-O), and the SCAN_POINTS - 1 constrained
# points step evenly to the adduct minimum, the last point.
SCAN_REACH_A = 1.5
SCAN_POINTS = 8


def peak_seeds(ctx: Ctx, name: str) -> tuple[Seed, ...]:
    """The highest peak detected on the latest profile, refined by a parabola (none without):
    a screen_hei seed on SCREEN's profile, else a path_hei one."""
    path = ctx.work.path
    peaks = () if path is None else profile.interior_maxima(path.energies, ctx.resolution)
    if path is None or not peaks:
        return ()
    _, _, x = profile.hei(path.frames, path.energies, max(peaks, key=path.energies.__getitem__))
    return (Seed(ctx.geometry(name, x), "screen_hei" if path.source == "screen" else "path_hei"),)


def _unavailable(reason: str) -> BarrierVerdict:
    return BarrierVerdict(verdict="unavailable", source="screen", reasons=(reason,))


def _xtb_ts_mode(ctx: Ctx, coords: np.ndarray) -> tuple[float, ...] | None:
    """A low-level TS's imaginary mode: its xTB freq has exactly one mode below -saddle_cm1."""
    freq, below = xtb_freq(ctx, coords), -ctx.rules.gates.saddle_cm1
    if freq is None or sum(nu < below for nu in freq.frequencies_cm1 or ()) != 1:
        return None
    return tuple(freq.imaginary_modes[0])  # ordered by frequency


def _shortcut(ctx: Ctx) -> tuple[BarrierVerdict, tuple[Seed, ...]] | None:
    """Step 1: the low-level TSs already known for this pair of states (discovery /
    mode-follow), in order. Each whose xTB freq confirms its mode gets a DFT SP and a
    three-point verdict [minimum, TS SP, minimum] (Ctx.classify: an open shell's with its
    ⟨S²⟩); the single-step ones seed the saddle search in that order, under the first one's
    verdict. None without a single-step one."""
    seeds = []
    for geometry in ctx.case.low_level_ts:
        if (mode := _xtb_ts_mode(ctx, ctx.coords(geometry))) is None:
            ctx.note("low_level_ts_rejected")
        else:
            seeds.append(Seed(geometry, "discovery_ts", mode))
    found = ctx.sps([ctx.coords(s.geometry) for s in seeds]) if seeds else []
    verdicts = [(ctx.classify((ctx.energies[0], ev.energy_hartree, ctx.energies[1]), (ev.s2,),
                              "screen"), seed)
                for seed, ev in zip(seeds, found, strict=True) if ev is not None]
    for verdict, _ in verdicts:
        if "scf_branch_jump" in verdict.reasons:
            ctx.note("shortcut:scf_branch_jump")
    singles = [(verdict, seed) for verdict, seed in verdicts if verdict.verdict == "single"]
    return (singles[0][0], tuple(seed for _, seed in singles)) if singles else None


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
    verdict = ctx.verdict(frames, [ev.energy_hartree for ev in inner], "screen",
                          [ev.s2 for ev in inner])
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
    """An association: the separated monomers → the adduct along a relaxed scan of
    the formed bond (i, j), from r_P + SCAN_REACH_A inwards (r_P: the adduct's). The first point
    is the adduct with the fragment holding j moved out along i→j, each next one the previous
    optimum moved in; each is optimized with r_ij fixed, its SCF started from the previous
    point's (a broken-symmetry pair stays on its continuous branch). The profile is
    [ΣE(monomers), scan..., E(adduct)] (the scan points' energies by ``_projected``); the
    longest point stands as the monomers' frame. A failed point leaves no profile: unavailable,
    no point dropped."""
    rt, adduct = ctx.rt, ctx.ends[1]
    [(i, j)] = ctx.change()[0]  # the one formed bond
    fragment = next(f for f in topology.fragments(ctx.symbols, ctx.ends[0]) if j in f)
    r_p = float(np.linalg.norm(adduct[j] - adduct[i]))
    ctx.note(f"scan:{i}-{j}:{r_p:.3f}+{SCAN_REACH_A}A:{SCAN_POINTS}_points")
    x, frames, points = adduct, [], []
    for k in range(SCAN_POINTS - 1):
        r = r_p + SCAN_REACH_A * (1.0 - k / (SCAN_POINTS - 1))
        opt = rt.qm.optimize(ctx.mol(_moved(x, (i, j), fragment, r)), rt.method,
                             fixed_bond=(i, j, r), scf_guess=points[-1] if points else None)
        if isinstance(opt, Failure):
            ctx.note(f"scan{k}:{opt.kind.value}")
            return BarrierVerdict(verdict="unavailable", source="scan",
                                  reasons=("scan_point",)), ()
        points.append(ctx.keep(opt))
        x = ctx.coords(opt.final)
        frames.append(x)
    energies = _projected(ctx, frames, points)
    frames = align_sequential([*frames, adduct])
    ctx.work.path = Profile([frames[0], *frames], (ctx.energies[0], *energies, ctx.energies[1]),
                            "scan")
    verdict = barrier_verdict(ctx.work.path.energies, source="scan", policy=ctx.rules.gates)
    return verdict, peak_seeds(ctx, "scan_hei") if verdict.verdict == "single" else ()


def projected(e_bs: float, s2_bs: float, e_hs: float, s2_hs: float, spin: float) -> float:
    """Yamaguchi's approximate spin projection (CPL 149, 537 (1988)): the energy of the pure
    spin-``spin`` state from a broken-symmetry solution and the high-spin one at its structure,
    E_LS = E_BS + α(E_BS − E_HS), α = (⟨S²⟩_BS − S(S+1)) / (⟨S²⟩_HS − ⟨S²⟩_BS)."""
    alpha = (s2_bs - spin * (spin + 1)) / (s2_hs - s2_bs)
    return e_bs + alpha * (e_bs - e_hs)


def _projected(ctx: Ctx, frames: list[np.ndarray], points: list[Evidence]) -> list[float]:
    """The scan points' energies: those of a low-spin-coupled pair (CH3· + O2 as a
    doublet, electronic_state.low_spin_coupled), broken-symmetry where the monomers separate,
    are spin-projected (``projected``) with an SP of the high-spin coupling Σ(m_i − 1) + 1 at
    each point's structure. A point whose high-spin SP fails or is itself contaminated (spin_ok)
    keeps its broken-symmetry energy, noted ap_skipped. The monomers' sum and the adduct keep
    theirs; spin_ok and the gates of stationary points and thermochemistry are untouched."""
    rt, energies = ctx.rt, [ev.energy_hartree for ev in points]
    spins = [rt.species[rt.registry.minima[m][0].species_id].multiplicity
             for m in ctx.case.monomers]
    if not low_spin_coupled(spins, ctx.multiplicity):
        return energies
    high, spin = sum(m - 1 for m in spins) + 1, (ctx.multiplicity - 1) / 2
    ctx.note(f"scan:ap:{high}")
    sps = rt.map(lambda x: rt.qm.energy(replace(ctx.mol(x), multiplicity=high), rt.method),
                 frames)
    for k, (bs, hs) in enumerate(zip(points, sps, strict=True)):
        if isinstance(hs, Failure):
            ctx.note(f"scan{k}:ap_skipped:{hs.kind.value}")
            continue
        ctx.keep(hs)
        if hs.s2 is None or bs.s2 is None or not spin_ok(hs, ctx.rules.gates):
            ctx.note(f"scan{k}:ap_skipped:s2={hs.s2}")
        else:
            energies[k] = projected(bs.energy_hartree, bs.s2, hs.energy_hartree, hs.s2, spin)
    return energies


def screen(ctx: Ctx, state: CaseState, decision: Decision) -> CaseState:
    """Barrier pre-check (design §7.3): an association's scan; else first the shortcut from the
    low-level TSs, the low-level path when none gives a single step or once their seeds have
    failed (a second SCREEN)."""
    found = None if ctx.case.monomers or state.screen is not None else _shortcut(ctx)
    if found is None:
        run = _scan if ctx.case.monomers else _screen_path
        found, state = run(ctx), replace(state, neb_done=True)
    verdict, seeds = found
    ctx.note(f"{verdict.source}:{verdict.verdict}:{','.join(verdict.reasons)}")
    return case_state.record_profile(state, verdict, seeds)


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


def find_path(ctx: Ctx, state: CaseState, decision: Decision) -> CaseState:
    """One DFT string chunk, classified at once (unconverged too: its maximum bounds the saddle
    and its peak is a seed); a single peak seeds the saddle search. The path gets the DFT
    minima back at its ends, which NWChem's freezeN moves."""
    rt, run_name = ctx.rt, f"string{state.path_runs}"
    initial = _initial_path(ctx, STRING_BEADS, f"{run_name}_initial")
    if initial is None:
        return replace(state, path_runs=state.path_runs + 1)
    run = rt.path.find_path(ctx.mol(ctx.ends[0]), ctx.mol(ctx.ends[1]), rt.method,
                            images=STRING_BEADS, initial_path=initial)
    if isinstance(run, Failure):
        ctx.note(f"{run_name}:{run.kind.value}")
        return replace(state, path_runs=state.path_runs + 1)
    frames = align_sequential([ctx.ends[0], *ctx.frames(run.images)[1:-1], ctx.ends[1]])
    verdict = ctx.verdict(frames, run.energies_hartree[1:-1], "string")
    seeds = peak_seeds(ctx, f"{run_name}_hei") if verdict.verdict == "single" else ()
    return case_state.record_profile(state, verdict, seeds)

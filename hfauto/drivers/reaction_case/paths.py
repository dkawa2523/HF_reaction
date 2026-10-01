"""Reaction-case path actions (design §7.3): SCREEN's low-level barrier pre-check and FIND_PATH's
DFT string chunk. Each classifies a DFT profile between the DFT minima and may seed the saddle
search at its peak."""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from hfauto.chemistry.gates import barrier_verdict
from hfauto.chemistry.interpolation import align_sequential, idpp, resample
from hfauto.core.evidence import Failure, FileRef, Geometry
from hfauto.core.records import BarrierVerdict
from hfauto.drivers.reaction_case import state as case_state
from hfauto.drivers.reaction_case.actions import Ctx, peak_seed, xtb_freq
from hfauto.drivers.reaction_case.state import CaseState, Decision, Seed

SCREEN_IMAGES = 11  # SCREEN's NEB images including both ends
STRING_BEADS = 9  # DFT string beads including both ends


def _unavailable(reason: str) -> BarrierVerdict:
    return BarrierVerdict(verdict="unavailable", source="screen", reasons=(reason,))


def _xtb_ts_mode(ctx: Ctx, coords: np.ndarray) -> tuple[float, ...] | None:
    """A low-level TS's imaginary mode: its xTB freq has exactly one mode below -saddle_cm1."""
    freq, below = xtb_freq(ctx, coords), -ctx.rules.gates.saddle_cm1
    if freq is None or sum(nu < below for nu in freq.frequencies_cm1 or ()) != 1:
        return None
    return tuple(freq.imaginary_modes[0])  # ordered by frequency


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
                              policy=ctx.rules.gates)
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
        frames, ts = _neb(ctx, idpp(ctx.symbols, *ctx.ends, SCREEN_IMAGES))
    except ValueError as exc:
        return _unavailable(f"idpp:{exc}"), None
    inner = ctx.sps(frames[1:-1])
    if any(e is None for e in inner):
        return _unavailable("screen_single_point"), None
    verdict = ctx.verdict(frames, [e for e in inner if e is not None], "screen")
    if verdict.verdict != "single":
        return verdict, None
    mode = None if ts is None else _xtb_ts_mode(ctx, ctx.coords(ts))
    if ts is not None and mode is not None:
        return verdict, Seed(ts, "screen_ts", mode)
    return verdict, peak_seed(ctx, "screen_hei", "screen_hei")


def screen(ctx: Ctx, state: CaseState, decision: Decision) -> CaseState:
    """Barrier pre-check (chem 16): first the shortcut from a low-level TS; the low-level path
    when there is none, it is not a single step, or its seed has failed (a second SCREEN)."""
    ts = ctx.case.low_level_ts
    found = None if ts is None or state.screen is not None else _shortcut(ctx, ts)
    if found is None or found[0].verdict != "single":
        found, state = _screen_path(ctx), replace(state, neb_done=True)
    verdict, seed = found
    ctx.note(f"screen:{verdict.verdict}:{','.join(verdict.reasons)}")
    return case_state.record_profile(state, verdict, seed)


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
                            images=STRING_BEADS, initial_path=initial, deadline=ctx.deadline)
    if isinstance(run, Failure):
        ctx.note(f"{run_name}:{run.kind.value}")
        return replace(state, path_runs=state.path_runs + 1)
    frames = align_sequential([ctx.ends[0], *ctx.frames(run.images)[1:-1], ctx.ends[1]])
    verdict = ctx.verdict(frames, run.energies_hartree[1:-1], "string")
    seed = peak_seed(ctx, f"{run_name}_hei", "path_hei") if verdict.verdict == "single" else None
    return case_state.record_profile(state, verdict, seed)

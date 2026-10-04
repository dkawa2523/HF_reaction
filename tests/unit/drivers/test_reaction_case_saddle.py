"""A reaction case's saddle search, TS checks and QRC (design §7.3) on the fake surfaces and case
builders of test_reaction_case_actions: one VALIDATE_AND_CONNECT action from the TS freq to QRC
(G7-P5), a stalled search restarted within its attempt unless it stalled above the latest DFT
profile (G1-P3), pushes and restarts bounded by the continuation depth (G3-P6), and no string
after one that left no seed (G3-P2)."""

from dataclasses import replace

import fakes
import numpy as np
import pytest
from unit.drivers.test_reaction_case_actions import (
    DFT,
    EXHAUSTED,
    OFF_AXIS,
    K,
    NoOpt,
    act,
    case_ctx,
)

from hfauto.chemistry import modes
from hfauto.chemistry.identity import mapped_equivalent, mapped_rmsd
from hfauto.chemistry.modes import BOUNDS_A, displace
from hfauto.core.evidence import Failure, FailureKind
from hfauto.core.records import BarrierVerdict, CaseOutcome
from hfauto.drivers.reaction_case import connection
from hfauto.drivers.reaction_case.actions import MAX_DEPTH, Profile
from hfauto.drivers.reaction_case.state import (
    Action,
    CaseRules,
    Decision,
    Seed,
    decide,
    record_profile,
)

VALIDATE = Decision(Action.VALIDATE_AND_CONNECT, "saddle_converged")


class HigherOrderQM(fakes.FakeQM):  # every freq job also reports a second imaginary mode
    def frequencies(self, mol, method, **kw):
        ev = super().frequencies(mol, method, **kw)
        nu = sorted(ev.frequencies_cm1)
        return ev.model_copy(update={"frequencies_cm1": (nu[0], -200.0, *nu[2:]),
                                     "imaginary_modes": (*ev.imaginary_modes, OFF_AXIS)})


class OffAxisQM(fakes.FakeQM):  # every freq job's one imaginary mode moves H off the A-H-B axis
    def frequencies(self, *args, **kw):
        return super().frequencies(*args, **kw).model_copy(update={"imaginary_modes": (OFF_AXIS,)})


class ScriptedSaddle(fakes.FakeSaddle):
    """Each search takes the next ``script`` entry (else "stall"): "ok" converges as the fake
    does, "stall" stops at maxiter like NWChem, its last frame the start moved 0.01 Å and its
    last energy ``energy``."""

    def __init__(self, root, pes, script=(), energy=None):
        super().__init__(root, pes)
        self.script, self.energy, self.starts = list(script), energy, []

    def refine(self, seed, method, **kw):
        self.starts.append(np.array(seed.xyz.coords))
        if (self.script.pop(0) if self.script else "stall") == "ok":
            return super().refine(seed, method, **kw)
        self.calls.append("refine")
        last = fakes.write_geometry(self.root, f"last{len(self.calls)}.xyz", seed.xyz.symbols,
                                    seed.xyz.coords + 0.01)
        return Failure(kind=FailureKind.GEOMETRY_MAXITER, reason="maxiter", final=last,
                       energy_hartree=self.energy)


def seeded(ctx, state, depth=0, attempts=0):
    """The state with one seed at the double well's TS (no mode: ρ, the xTB Hessian)."""
    ts = ctx.geometry("ts", ctx.rt.qm.pes.points["ts"])
    return replace(state, seeds=(Seed(ts, "screen_hei", depth=depth),), saddle_attempts=attempts)


def qrc_step(ctx, freq) -> float:
    """The largest atomic displacement of the latest QRC's + side from its TS."""
    side = ctx.work.calcs[ctx.work.connection.side_calcs[0]]
    return float(np.linalg.norm(ctx.coords(side.start) - ctx.coords(freq.final), axis=1).max())


# VALIDATE_AND_CONNECT ------------------------------------------------------------------------

def test_saddle_hessians_validation_and_qrc_on_a_double_well(tmp_path, monkeypatch) -> None:
    ctx, state = case_ctx(tmp_path, fakes.double_well())
    state = act(ctx, state, Action.FIND_PATH)
    seeds, n_freq = state.seeds, ctx.rt.qm.calls.count("frequencies")
    ctx.log = (logged := []).append
    ctx.rt = replace(ctx.rt, screen_qm=(xtb := HigherOrderQM(tmp_path, ctx.rt.qm.pes)))
    state = act(ctx, state, Action.REFINE_SADDLE)  # xTB whenever there is one, 2 negative modes too
    assert state.last_saddle == "converged" and xtb.calls == ["frequencies"]
    assert ctx.rt.qm.calls.count("frequencies") == n_freq and ctx.rt.saddle.calls == ["refine"]
    ctx.rt = replace(ctx.rt, screen_qm=None)  # else DFT
    state = act(ctx, replace(state, seeds=seeds), Action.REFINE_SADDLE)
    assert ctx.rt.qm.calls.count("frequencies") == n_freq + 1 and state.saddle_attempts == 2
    notes = [r["note"] for r in logged if "saddle_hessian" in r["note"]]
    assert notes == ["saddle_hessian:xtb:rho", "saddle_hessian:dft:rho"]
    assert decide(ctx.case, state, ctx.rules) == VALIDATE
    jobs = len(ctx.rt.qm.calls)
    state = act(ctx, state, Action.VALIDATE_AND_CONNECT)  # a separate freq job, then QRC
    claim, connected = state.claim, ctx.work.connection
    assert claim.imag_cm1 < -50 and not claim.notes and claim.freq_calc != claim.saddle_calc
    assert (state.last_saddle, state.connection) == (None, "elementary")
    assert set(connected.minima) == set(ctx.case.minima)
    # one amplitude; both sides from the TS Hessian, run at once
    assert ctx.rt.qm.calls[jobs:] == ["frequencies"] + ["optimize+init_hessian"] * 2
    assert len(set(connected.side_calcs)) == 2 == ctx.rt.map[-1]
    assert "minus_is_image" not in str(logged)
    freq = ctx.work.calcs[claim.freq_calc]
    first = qrc_step(ctx, freq)
    connection.connect(ctx, state, freq, 2)  # the second amplitude: × 2
    assert qrc_step(ctx, freq) == pytest.approx(2 * first, abs=1e-6)
    monkeypatch.setattr(modes, "BOUNDS_A", (0.03, 0.05))  # the first one is capped
    amplitudes = []
    for attempt in (1, 2):  # C19: × 2, then clipped, so the 2nd amplitude equals the 1st
        amplitudes.append((connection.connect(ctx, state, freq, attempt).connection,
                           qrc_step(ctx, freq)))
    assert amplitudes == [("elementary", pytest.approx(0.05, abs=1e-6))] * 2


def test_a_symmetric_ts_optimizes_one_qrc_side_and_carries_its_image(tmp_path) -> None:
    """CA-1b: at the F-H-F TS the minus start is the plus start with F1 and F2 exchanged, so
    only the plus side is optimized; its image is the minus side's own optimum as labelled and
    makes a degenerate step with distinct sides, both claimed by the one optimization."""
    pes = fakes.symmetric_double_well()
    ctx, state = case_ctx(tmp_path, pes)
    freq = ctx.rt.qm.frequencies(pes.molecule("ts"), DFT)
    starts = displace(pes.points["ts"], np.asarray(freq.imaginary_modes[0]), 0.1)
    (plus, minus), (x_plus, x_minus) = connection._sides(ctx, freq, starts, 1)
    alone = ctx.coords(ctx.rt.qm.optimize(ctx.mol(starts[1]), DFT, init_hessian=freq).final)
    assert plus is minus and mapped_rmsd(x_minus, alone) < 1e-3
    assert mapped_equivalent(pes.symbols, x_plus, x_minus)  # H at F2, H at F1: one basin
    jobs, ctx.log = len(ctx.rt.qm.calls), (logged := []).append
    state = connection.connect(ctx, state, freq, 1)
    assert ctx.rt.qm.calls[jobs:] == ["optimize+init_hessian"] and state.connection == "degenerate"
    assert len(set(ctx.work.connection.side_calcs)) == 1 == ctx.rt.map[-1]
    assert {"note": "qrc1:minus_is_image"} in logged


def test_a_qrc_side_in_a_new_basin_is_registered_from_its_own_optimization(tmp_path) -> None:
    """K5: the converged QRC side toward the intermediate is the new minimum's opt; registering
    it adds only its freq job, no second optimization."""
    pes = fakes.triple_well()
    ctx, state = case_ctx(tmp_path, pes)
    freq, jobs = ctx.rt.qm.frequencies(pes.molecule("ts1"), DFT), len(ctx.rt.qm.calls)
    state, claim = connection.connect(ctx, state, freq, 1), ctx.work.connection
    assert ctx.rt.qm.calls[jobs:] == ["optimize+init_hessian"] * 2 + ["frequencies"]
    (well,) = ctx.work.minima.values()
    assert state.connection == "reassigned" and well.minimum_id in claim.minima
    assert well.opt_calc in claim.side_calcs


@pytest.mark.parametrize("screen,seeds,action", [
    (True, 0, Action.SCREEN), (False, 0, Action.FIND_PATH), (False, 1, Action.REFINE_SADDLE)])
def test_a_saddle_without_an_imaginary_mode_is_a_failed_attempt(tmp_path, screen, seeds, action):
    """G3-P1: a saddle search that converged onto a minimum (the reactant) is rejected like any
    other saddle: no QRC, no well relaxed from it, and SCREEN, a string or the next seed goes
    on."""
    ctx, state = case_ctx(tmp_path, fakes.double_well())
    ctx.rules, ctx.log = CaseRules(screen=screen), (logged := []).append
    ctx.work.saddle = ctx.rt.qm.optimize(ctx.mol(ctx.ends[0]), DFT)
    seed = Seed(ctx.work.saddle.final, "path_hei")
    state = replace(state, seeds=(seed,) * seeds, saddle_attempts=1, last_saddle="converged")
    assert decide(ctx.case, state, ctx.rules) == VALIDATE
    jobs = len(ctx.rt.qm.calls)
    state = act(ctx, state, Action.VALIDATE_AND_CONNECT)
    assert logged == [{"note": "ts_rejected:no_imaginary_mode"}]
    assert ctx.rt.qm.calls[jobs:] == ["frequencies"] and state.claim is None
    assert (state.last_saddle, state.saddle_attempts, state.connection) == ("failed", 1, None)
    assert decide(ctx.case, state, ctx.rules).action is action


def test_a_saddle_whose_mode_leaves_the_bond_change_gets_no_qrc(tmp_path) -> None:
    """G1-P2: chi on the labelled ends' bonds (F-H-F: one basin); H off the axis gets no QRC."""
    ctx, state = case_ctx(tmp_path, fakes.symmetric_double_well())
    state = act(ctx, act(ctx, state, Action.FIND_PATH), Action.REFINE_SADDLE)
    ctx.rt, ctx.log = replace(ctx.rt, qm=OffAxisQM(tmp_path, ctx.rt.qm.pes)), (logged := []).append
    state = act(ctx, state, Action.VALIDATE_AND_CONNECT)
    assert logged == [{"note": "ts_rejected:not_reaction_mode"}] and state.saddle_attempts == 1
    assert ctx.rt.qm.calls == ["frequencies"] and state.claim is None
    assert decide(ctx.case, state, ctx.rules) == Decision(Action.FIND_PATH, "dft_path")


# restart, push and the continuation depth ------------------------------------------------------

def test_a_higher_order_saddle_is_pushed_once_and_refined_from_its_ts_hessian(tmp_path):
    """U6-P3: the reaction mode is the negative mode along ρ; the other modes below
    -saddle_cm1 (here the second) push the saddle once by the energy target (modes.off_saddle),
    and the TS freq is the seed's Hessian. The push is a counted attempt one continuation
    deeper."""
    ctx, state = case_ctx(tmp_path, fakes.double_well())
    state = act(ctx, act(ctx, state, Action.FIND_PATH), Action.REFINE_SADDLE)
    x = ctx.coords(ctx.work.saddle.final)
    ctx.rt = replace(ctx.rt, qm=HigherOrderQM(tmp_path, ctx.rt.qm.pes))
    state = act(ctx, state, Action.VALIDATE_AND_CONNECT)
    seed = state.seeds[0]
    assert (state.last_saddle, state.claim, seed.source, seed.depth) == (
        "failed", None, "higher_order_retry", 1)
    assert seed.hessian.task == "freq" and seed.mode == seed.hessian.imaginary_modes[0]
    push = (ctx.coords(seed.geometry) - x).ravel()  # H along y: the second mode
    assert np.flatnonzero(np.abs(push) > 1e-9).tolist() == [4]
    assert BOUNDS_A[0] < abs(push[4]) < BOUNDS_A[1]
    decision = decide(ctx.case, state, ctx.rules)
    assert decision == Decision(Action.REFINE_SADDLE, "seed:higher_order_retry")
    ctx.log = (logged := []).append
    jobs = ctx.rt.qm.calls.count("frequencies"), len(ctx.rt.screen_qm.calls)
    state = act(ctx, state, Action.REFINE_SADDLE)  # within 0.5 Å of the TS freq: no new Hessian
    assert state.last_saddle == "converged" and state.saddle_attempts == 2
    assert [r["note"] for r in logged] == ["saddle_hessian:ts_freq:mode"]
    assert (ctx.rt.qm.calls.count("frequencies"), len(ctx.rt.screen_qm.calls)) == jobs
    assert ctx.work.depth == 1


def test_a_stalled_saddle_restarts_once_within_its_attempt(tmp_path) -> None:
    """U6-P6 / R3: maxiter restarts the search once, inside its attempt, from its last frame
    with a fresh Hessian, so the last attempt still restarts; a restart that stalls again ends
    the attempt, and a used budget ends the case without a string."""
    ctx, state = case_ctx(tmp_path, fakes.double_well(), saddle=ScriptedSaddle(tmp_path, None))
    xtb, ctx.log = ctx.rt.screen_qm.calls.count("frequencies"), (logged := []).append
    state = act(ctx, seeded(ctx, state, attempts=1), Action.REFINE_SADDLE)
    saddle = ctx.rt.saddle
    assert saddle.calls == ["refine"] * 2 and np.allclose(saddle.starts[1], saddle.starts[0] + 0.01)
    assert ctx.rt.screen_qm.calls.count("frequencies") == xtb + 2  # an xTB Hessian at each start
    assert [r["note"] for r in logged] == ["saddle_hessian:xtb:rho",
                                           "saddle:geometry_maxiter:maxiter"] * 2
    assert (state.seeds, state.saddle_attempts, state.last_saddle) == ((), 2, "failed")
    assert decide(ctx.case, state, ctx.rules) == EXHAUSTED


@pytest.mark.parametrize("above_kcal,profile,restarts", [
    (1.5, "string", False),  # climbed past the profile's maximum
    (0.5, "string", True),  # within a resolution of it
    (None, "string", True),  # a stored stalled search without its energy: no bound
    (1.5, None, True),  # no DFT profile yet (a low-level TS or ts_calc): no bound
])
def test_a_stalled_search_above_the_latest_profile_is_not_restarted(tmp_path, above_kcal, profile,
                                                                    restarts):
    """G1-P3: a continuous DFT path bounds the saddle from above, so a search that stalled more
    than a resolution above the latest profile's maximum has climbed past the barrier: no
    restart, noted with the profile that bounds it. VAL7 s5: frames 29-32 kcal/mol above their
    HEI seed."""
    peak = 6.0 * K  # the profile's maximum, an interior node
    energy = None if above_kcal is None else peak + above_kcal * K
    ctx, state = case_ctx(tmp_path, fakes.double_well(),
                          saddle=ScriptedSaddle(tmp_path, None, energy=energy))
    a, b = ctx.ends
    frames = [a + t * (b - a) for t in np.linspace(0.0, 1.0, 5)]
    ctx.work.path = None if profile is None else Profile(
        frames, tuple(e * K for e in (0, 3, 6, 2, 1)), profile)
    ctx.log = (logged := []).append
    state = act(ctx, seeded(ctx, state), Action.REFINE_SADDLE)
    assert len(ctx.rt.saddle.calls) == 1 + restarts and state.saddle_attempts == 1
    bounded = {"note": "saddle:above_path_bound:string"} in logged
    assert bounded is (above_kcal == 1.5 and profile is not None)


def test_a_push_whose_search_stalls_restarts_and_is_connected(tmp_path) -> None:
    """G3-P6 (W3 s6 split2_split2): a second-order saddle → its push (depth 1, the second
    attempt) → maxiter → the restart from its last frame (depth 2, the same attempt) → a TS
    that QRC connects to the ends."""
    pes = fakes.double_well()
    ctx, state = case_ctx(tmp_path, pes, saddle=ScriptedSaddle(tmp_path, pes, ["ok", "stall", "ok"]))
    qm = ctx.rt.qm
    state = act(ctx, seeded(ctx, state), Action.REFINE_SADDLE)
    ctx.rt = replace(ctx.rt, qm=HigherOrderQM(tmp_path, pes))
    state = act(ctx, state, Action.VALIDATE_AND_CONNECT)
    assert [(s.source, s.depth) for s in state.seeds] == [("higher_order_retry", 1)]
    ctx.rt = replace(ctx.rt, qm=qm)
    assert decide(ctx.case, state, ctx.rules) == Decision(Action.REFINE_SADDLE,
                                                          "seed:higher_order_retry")
    state = act(ctx, state, Action.REFINE_SADDLE)
    assert (state.last_saddle, state.saddle_attempts, ctx.work.depth) == ("converged", 2, 2)
    assert len(ctx.rt.saddle.calls) == 3
    state = act(ctx, state, decide(ctx.case, state, ctx.rules).action)
    assert decide(ctx.case, state, ctx.rules) == Decision(
        Action.COMPLETE, "connection:elementary", CaseOutcome.ELEMENTARY_STEP)


def test_the_continuation_depth_bounds_restarts_and_pushes(tmp_path) -> None:
    """G3-P6: a restart and a push each add one continuation, up to MAX_DEPTH, so restarts and
    pushes cannot alternate without end: a depth-2 push that stalls is not restarted, and a
    second-order saddle that a depth-2 search found is not pushed."""
    assert MAX_DEPTH == 2
    pes = fakes.double_well()
    ctx, state = case_ctx(tmp_path, pes, saddle=ScriptedSaddle(tmp_path, pes, ["stall", "ok"]))
    state = act(ctx, seeded(ctx, state), Action.REFINE_SADDLE)  # stalls, restarts (depth 1)
    assert ctx.work.depth == 1 and len(ctx.rt.saddle.calls) == 2
    ctx.rt = replace(ctx.rt, qm=HigherOrderQM(tmp_path, pes))
    state = act(ctx, state, Action.VALIDATE_AND_CONNECT)
    assert [(s.source, s.depth) for s in state.seeds] == [("higher_order_retry", 2)]
    state = act(ctx, state, Action.REFINE_SADDLE)  # the push stalls: no restart at depth 2
    assert len(ctx.rt.saddle.calls) == 3 and state.last_saddle == "failed"
    assert decide(ctx.case, state, ctx.rules) == EXHAUSTED
    ctx, state = case_ctx(tmp_path / "b", pes,
                          saddle=ScriptedSaddle(tmp_path / "b", pes, ["stall", "ok"]))
    state = act(ctx, seeded(ctx, state, depth=1), Action.REFINE_SADDLE)  # a push, restarted
    ctx.rt, ctx.log = replace(ctx.rt, qm=HigherOrderQM(tmp_path / "b", pes)), (logged := []).append
    state = act(ctx, state, Action.VALIDATE_AND_CONNECT)
    assert ctx.work.depth == 2 and not state.seeds and state.last_saddle == "failed"
    assert logged == [{"note": "ts_rejected:higher_order"}]


# the merged FIND_PATH row -----------------------------------------------------------------------

def test_a_string_that_left_no_seed_ends_the_case(tmp_path) -> None:
    """G3-P2: a string whose well cannot be relaxed and whose profile rises monotonically from
    it to the product leaves no seed; no next chunk runs (before: up to three), the case ends."""
    ctx, state = case_ctx(tmp_path, fakes.double_well())
    ctx.rt = replace(ctx.rt, qm=NoOpt(tmp_path, ctx.rt.qm.pes))
    a, b = ctx.ends
    frames = [a + t * (b - a) for t in np.linspace(0.0, 1.0, 5)]
    ctx.work.path = Profile(frames, tuple(e * K for e in (0, -3, 1, 2, 4)), "string")
    state = record_profile(state, BarrierVerdict(verdict="intermediate", source="string"))
    state = act(ctx, state, Action.VALIDATE_INTERMEDIATE, "path_intermediate")
    assert (state.intermediate, state.seeds, state.path_runs) == ("relax_failed", (), 1)
    assert decide(ctx.case, state, ctx.rules) == EXHAUSTED
    tried = replace(state, saddle_attempts=1)  # had its seed been tried, the next chunk runs
    assert decide(ctx.case, tried, ctx.rules) == Decision(Action.FIND_PATH, "dft_path")

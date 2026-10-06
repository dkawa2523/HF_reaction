"""A reaction case's saddle search, TS checks and QRC (design §7.3) on the fake surfaces and case
builders of test_reaction_case_actions: one VALIDATE_AND_CONNECT action from the TS freq to QRC
(G7-P5), the certification of a TS (X1), the one continuation of a search (U6-P3: a stalled
search, a higher-order or non-stationary saddle), and every failure as one token (U6-P4)."""

from dataclasses import replace

import fakes
import numpy as np
import pytest
from scipy.spatial.distance import pdist
from unit.drivers.test_reaction_case_actions import (
    DFT,
    EXHAUSTED,
    OFF_AXIS,
    FailingSaddle,
    K,
    NoOpt,
    act,
    case_ctx,
)

from hfauto.chemistry.modes import BOUNDS_A
from hfauto.chemistry.profile import Profile
from hfauto.core.evidence import Failure, FailureKind
from hfauto.core.records import BarrierVerdict, CaseOutcome
from hfauto.drivers.reaction_case import connection
from hfauto.drivers.reaction_case.state import (
    Action,
    CaseRules,
    Decision,
    Seed,
    decide,
    record_profile,
)

VALIDATE = Decision(Action.VALIDATE_AND_CONNECT, "saddle_converged")
CONTINUE = Decision(Action.REFINE_SADDLE, "seed:continuation")


class HigherOrderQM(fakes.FakeQM):  # every freq job also reports a second imaginary mode
    def frequencies(self, mol, method, **kw):
        ev = super().frequencies(mol, method, **kw)
        nu = sorted(ev.frequencies_cm1)
        return ev.model_copy(update={"frequencies_cm1": (nu[0], -200.0, *nu[2:]),
                                     "imaginary_modes": (*ev.imaginary_modes, OFF_AXIS)})


class SteppedQM(HigherOrderQM):  # ν2 at the first freq only: the saddle's
    def frequencies(self, mol, method, **kw):
        self.freqs = getattr(self, "freqs", 0) + 1
        if self.freqs == 1:
            return super().frequencies(mol, method, **kw)
        return fakes.FakeQM.frequencies(self, mol, method, **kw)


class OffAxisQM(fakes.FakeQM):  # every freq job's one imaginary mode moves H off the A-H-B axis
    def frequencies(self, *args, **kw):
        return super().frequencies(*args, **kw).model_copy(update={"imaginary_modes": (OFF_AXIS,)})


class ScriptedSaddle(fakes.FakeSaddle):
    """Each search takes the next ``script`` entry (else "stall"): "ok" converges as the fake
    does, "stall" stops at maxiter like NWChem, its last frame the start moved 0.01 Å and its
    last energy ``energy``."""

    def __init__(self, root, pes, script=(), energy=None):
        super().__init__(root, pes)
        self.script, self.energy, self.starts, self.guesses = list(script), energy, [], []

    def refine(self, seed, method, *, scf_guess=None, **kw):
        self.starts.append(np.array(seed.xyz.coords))
        self.guesses.append(scf_guess)
        if (self.script.pop(0) if self.script else "stall") == "ok":
            return super().refine(seed, method, scf_guess=scf_guess, **kw)
        self.calls.append("refine")
        last = fakes.write_geometry(self.root, f"last{len(self.calls)}.xyz", seed.xyz.symbols,
                                    seed.xyz.coords + 0.01)
        return Failure(kind=FailureKind.GEOMETRY_MAXITER, reason="maxiter", final=last,
                       energy_hartree=self.energy)


def seeded(ctx, state, depth=0, attempts=0):
    """The state with one seed at the double well's TS (no mode: ρ, the xTB Hessian)."""
    ts = ctx.geometry("ts", ctx.rt.qm.pes.points["ts"])
    return replace(state, seeds=(Seed(ts, "screen_hei", depth=depth),), attempts=attempts)


def off_ridge(ctx, bend=0.1):
    """The converged saddle replaced by the TS with H bent ``bend`` Å off the A-H-B axis, with
    its own gradient: a first-order saddle that fails the certification."""
    pes = ctx.rt.qm.pes
    x = pes.points["ts"] + [[0, 0, 0], [0, bend, 0], [0, 0, 0]]
    opt = ctx.rt.qm.optimize(ctx.mol(x), DFT)
    final = fakes.write_geometry(ctx.rt.qm.root, f"bent{bend}.xyz", pes.symbols, x)
    ctx.work.saddle = opt.model_copy(update={
        "task": "saddle", "final": final, "energy_hartree": pes.energy(x),
        "gradient": tuple(pes.gradient(x) * fakes.BOHR)})
    return x


# VALIDATE_AND_CONNECT ------------------------------------------------------------------------

def test_saddle_hessians_validation_and_qrc_on_a_double_well(tmp_path) -> None:
    ctx, state = case_ctx(tmp_path, fakes.double_well())
    state = act(ctx, state, Action.FIND_PATH)
    seeds, n_freq = state.seeds, ctx.rt.qm.calls.count("frequencies")
    ctx.log = (logged := []).append
    ctx.rt = replace(ctx.rt, screen_qm=(xtb := HigherOrderQM(tmp_path, ctx.rt.qm.pes)))
    state = act(ctx, state, Action.REFINE_SADDLE)  # xTB whenever there is one, 2 negative modes too
    assert state.saddle_pending and xtb.calls == ["frequencies"]
    assert ctx.rt.qm.calls.count("frequencies") == n_freq and ctx.rt.saddle.calls == ["refine"]
    ctx.rt = replace(ctx.rt, screen_qm=None)  # else DFT
    state = act(ctx, replace(state, seeds=seeds), Action.REFINE_SADDLE)
    assert ctx.rt.qm.calls.count("frequencies") == n_freq + 1 and state.attempts == 2
    notes = [r["note"] for r in logged if "saddle_hessian" in r["note"]]
    assert notes == ["saddle_hessian:xtb:rho", "saddle_hessian:dft:rho"]
    assert decide(ctx.case, state, ctx.rules) == VALIDATE
    jobs = len(ctx.rt.qm.calls)
    state = act(ctx, state, Action.VALIDATE_AND_CONNECT)  # a separate freq job, then QRC
    claim, connected = state.claim, ctx.work.connection
    assert claim.imag_cm1 < -50 and not claim.notes and claim.freq_calc != claim.saddle_calc
    assert (state.saddle_pending, state.connection) == (False, "elementary")
    assert set(connected.minima) == set(ctx.case.minima)
    # one amplitude; both sides from the TS Hessian as the positive model, run at once, each
    # in a known basin (no freq)
    assert ctx.rt.qm.calls[jobs:] == ["frequencies"] + ["optimize+init_hessian"] * 2
    assert ctx.rt.qm.models[-2:] == ["positive"] * 2
    assert len(set(connected.side_calcs)) == 2 == ctx.rt.map[-1]
    assert "minus_is_image" not in str(logged)


def test_a_symmetric_ts_optimizes_one_qrc_side_and_carries_its_image(tmp_path) -> None:
    """CA-1b: at the F-H-F TS the minus start is the plus start with F1 and F2 exchanged, so
    only the plus side is optimized; its image is the minus side as labelled and makes a
    degenerate step with distinct sides, both claimed by the one optimization."""
    pes = fakes.symmetric_double_well()
    ctx, state = case_ctx(tmp_path, pes)
    freq = ctx.rt.qm.frequencies(pes.molecule("ts"), DFT)
    jobs, ctx.log = len(ctx.rt.qm.calls), (logged := []).append
    state = connection.connect(ctx, state, freq)
    assert ctx.rt.qm.calls[jobs:] == ["optimize+init_hessian"] and state.connection == "degenerate"
    assert len(set(ctx.work.connection.side_calcs)) == 1 == ctx.rt.map[-1]
    assert {"note": "qrc:minus_is_image"} in logged


def test_a_qrc_side_in_a_new_basin_is_registered_from_its_own_optimization(tmp_path) -> None:
    """K5: the converged QRC side toward the intermediate is the new minimum's opt; registering
    it adds only its freq job, no second optimization."""
    pes = fakes.triple_well()
    ctx, state = case_ctx(tmp_path, pes)
    freq, jobs = ctx.rt.qm.frequencies(pes.molecule("ts1"), DFT), len(ctx.rt.qm.calls)
    state, claim = connection.connect(ctx, state, freq), ctx.work.connection
    assert ctx.rt.qm.calls[jobs:] == ["optimize+init_hessian"] * 2 + ["frequencies"]
    (well,) = ctx.work.minima.values()
    assert state.connection == "reassigned" and well.minimum_id in claim.minima
    assert well.opt_calc in claim.side_calcs


def test_a_qrc_side_on_a_saddle_goes_down_once_more(tmp_path) -> None:
    """W3 (S5 split2): a QRC side that stops on a saddle (here the triple well's side from ts1
    toward the intermediate, stopped on ts2) continues down once, one side along its imaginary
    mode away from the TS, so the case answers deterministically instead of unassigned_side."""
    pes = fakes.triple_well()
    ctx, state = case_ctx(tmp_path, pes)
    freq = ctx.rt.qm.frequencies(pes.molecule("ts1"), DFT)
    ts1, ts2 = (pes.points[n][1, 0] for n in ("ts1", "ts2"))

    class OnTheTop(fakes.FakeQM):  # a side heading for the intermediate stops on ts2
        def optimize(self, mol, method, **kw):
            onward = kw.get("init_hessian") is freq and mol.xyz.coords[1, 0] > ts1
            return super().optimize(pes.molecule("ts2") if onward else mol, method, **kw)

    ctx.rt = replace(ctx.rt, qm=OnTheTop(tmp_path, pes))
    state = connection.connect(ctx, state, freq)
    assert state.connection == "elementary" and set(ctx.work.connection.minima) == set(
        ctx.case.minima)
    assert ctx.rt.qm.calls.count("frequencies") == 1  # ts2's: the sides are known basins
    side = ctx.work.calcs[ctx.work.connection.side_calcs[1]]
    assert ctx.coords(side.final)[1, 0] > ts2  # past ts2, away from ts1


@pytest.mark.parametrize("screen,seeds,action", [
    (True, 0, Action.SCREEN), (False, 0, Action.FIND_PATH), (False, 1, Action.REFINE_SADDLE)])
def test_a_saddle_without_an_imaginary_mode_is_a_failed_attempt(tmp_path, screen, seeds, action):
    """G3-P1: a saddle search that converged onto a minimum (the reactant) is rejected like any
    other saddle: no QRC, no continuation, no well relaxed from it, and SCREEN, a string or the
    next seed goes on."""
    ctx, state = case_ctx(tmp_path, fakes.double_well())
    ctx.rules, ctx.log = CaseRules(screen=screen), (logged := []).append
    ctx.work.saddle = ctx.rt.qm.optimize(ctx.mol(ctx.ends[0]), DFT)
    seed = Seed(ctx.work.saddle.final, "path_hei")
    state = replace(state, seeds=(seed,) * seeds, attempts=1, saddle_pending=True)
    assert decide(ctx.case, state, ctx.rules) == VALIDATE
    jobs = len(ctx.rt.qm.calls)
    state = act(ctx, state, Action.VALIDATE_AND_CONNECT)
    assert logged == [{"note": "ts_rejected:no_imaginary_mode"}]
    assert ctx.rt.qm.calls[jobs:] == ["frequencies"] and state.claim is None
    assert (state.saddle_pending, state.attempts, state.connection) == (False, 1, None)
    assert decide(ctx.case, state, ctx.rules).action is action


def test_a_saddle_whose_mode_leaves_the_bond_change_gets_no_qrc(tmp_path) -> None:
    """G1-P2: chi on the labelled ends' bond change (F-H-F: one basin); H off the axis gets no
    QRC and no continuation, its chi kept in the token."""
    ctx, state = case_ctx(tmp_path, fakes.symmetric_double_well())
    state = act(ctx, act(ctx, state, Action.FIND_PATH), Action.REFINE_SADDLE)
    assert ctx.change() == (frozenset({(1, 2)}), frozenset({(0, 1)}))
    ctx.rt, ctx.log = replace(ctx.rt, qm=OffAxisQM(tmp_path, ctx.rt.qm.pes)), (logged := []).append
    state = act(ctx, state, Action.VALIDATE_AND_CONNECT)
    [note] = [r["note"] for r in logged]
    assert note.startswith("ts_rejected:not_reaction_mode:") and float(note[-5:]) < 0.01
    assert state.attempts == 1 and not state.seeds
    assert ctx.rt.qm.calls == ["frequencies"] and state.claim is None
    assert decide(ctx.case, state, ctx.rules) == Decision(Action.FIND_PATH, "dft_path")


# the one continuation -----------------------------------------------------------------------

def test_a_second_order_saddle_reaches_a_stationary_ts_through_the_continuation(tmp_path):
    """U6-P3: a stationary second-order saddle is continued once: its verified freq is the
    seed's Hessian and SCF guess, the imaginary mode of largest χ its mode, and the other mode
    below -saddle_cm1 pushes the start once (modes.off_saddle). The continuation is the next,
    counted seed; its saddle is the stationary TS, connected."""
    pes = fakes.double_well()
    ctx, state = case_ctx(tmp_path, pes, saddle=ScriptedSaddle(tmp_path, pes, ["ok"] * 2))
    state = act(ctx, act(ctx, state, Action.FIND_PATH), Action.REFINE_SADDLE)
    x = ctx.coords(ctx.work.saddle.final)
    ctx.rt, ctx.log = replace(ctx.rt, qm=SteppedQM(tmp_path, pes)), (logged := []).append
    state = act(ctx, state, Action.VALIDATE_AND_CONNECT)
    (seed,) = state.seeds
    assert (state.claim, seed.source, seed.depth) == (None, "continuation", 1)
    assert seed.hessian.task == "freq" and seed.mode == seed.hessian.imaginary_modes[0]
    push = (ctx.coords(seed.geometry) - x).ravel()  # H along y: the second mode
    assert np.flatnonzero(np.abs(push) > 1e-9).tolist() == [4]
    assert BOUNDS_A[0] < abs(push[4]) < BOUNDS_A[1]
    assert decide(ctx.case, state, ctx.rules) == CONTINUE
    jobs = ctx.rt.qm.calls.count("frequencies"), len(ctx.rt.screen_qm.calls)
    state = act(ctx, state, Action.REFINE_SADDLE)  # within 0.5 Å of its freq: no new Hessian
    assert state.saddle_pending and state.attempts == 2 and ctx.work.depth == 1
    assert "saddle_hessian:seed_freq:mode" in [r.get("note") for r in logged]
    assert (ctx.rt.qm.calls.count("frequencies"), len(ctx.rt.screen_qm.calls)) == jobs
    assert ctx.rt.saddle.guesses == [None, seed.hessian]
    state = act(ctx, state, Action.VALIDATE_AND_CONNECT)
    assert state.connection == "elementary" and not state.claim.notes


def test_a_non_stationary_ts_is_continued_from_where_it_stopped(tmp_path):
    """X1: a first-order saddle off its ridge (H bent 0.1 Å) passes the frequency gates but not
    the certification (ΔE_TR > BASIN_DE_HARTREE): continued once from that structure (its
    gradient leads the search), it reaches the stationary TS, connected."""
    pes = fakes.double_well()
    ctx, state = case_ctx(tmp_path, pes)
    state = act(ctx, act(ctx, state, Action.FIND_PATH), Action.REFINE_SADDLE)
    x, ctx.log = off_ridge(ctx), (logged := []).append
    state = act(ctx, state, Action.VALIDATE_AND_CONNECT)
    (seed,) = state.seeds
    assert {"note": "ts_rejected:not_stationary"} in logged and state.claim is None
    assert (seed.source, seed.depth) == ("continuation", 1) and np.allclose(
        ctx.coords(seed.geometry), x)
    assert decide(ctx.case, state, ctx.rules) == CONTINUE
    state = act(ctx, act(ctx, state, Action.REFINE_SADDLE), Action.VALIDATE_AND_CONNECT)
    assert state.connection == "elementary" and "not_stationary" not in state.claim.notes
    assert pdist(ctx.coords(ctx.work.saddle.final)) == pytest.approx(pdist(pes.points["ts"]),
                                                                    abs=1e-4)


def test_a_continued_saddle_is_continued_no_more(tmp_path):
    """MAX_DEPTH is gone: a continuation (depth 1) that is still not stationary is accepted,
    noted not_stationary (thermo blocks it); one that is higher-order is rejected without a
    seed; one that stalls is not continued."""
    pes = fakes.double_well()
    ctx, state = case_ctx(tmp_path, pes)
    state = act(ctx, seeded(ctx, state, depth=1), Action.REFINE_SADDLE)
    off_ridge(ctx)
    state = act(ctx, state, Action.VALIDATE_AND_CONNECT)
    assert "not_stationary" in state.claim.notes and state.connection == "elementary"
    ctx, state = case_ctx(tmp_path / "b", pes)
    state = act(ctx, seeded(ctx, state, depth=1), Action.REFINE_SADDLE)
    ctx.rt, ctx.log = replace(ctx.rt, qm=HigherOrderQM(tmp_path / "b", pes)), (logged := []).append
    state = act(ctx, state, Action.VALIDATE_AND_CONNECT)
    assert logged == [{"note": "ts_rejected:higher_order"}] and not state.seeds
    ctx, state = case_ctx(tmp_path / "c", pes, saddle=ScriptedSaddle(tmp_path / "c", pes))
    state = act(ctx, seeded(ctx, state, depth=1), Action.REFINE_SADDLE)
    assert len(ctx.rt.saddle.calls) == 1 and not state.saddle_pending


def test_a_stalled_search_continues_once_within_its_attempt_from_a_dft_freq(tmp_path) -> None:
    """U6-P3: maxiter continues the search once, inside its attempt, from its last frame with
    the DFT freq there (no xTB Hessian): the last attempt still continues; a continuation that
    stalls again ends the attempt, and a used budget ends the case without a string."""
    pes = fakes.double_well()
    ctx, state = case_ctx(tmp_path, pes, saddle=ScriptedSaddle(tmp_path, pes))
    xtb, dft = (q.calls.count("frequencies") for q in (ctx.rt.screen_qm, ctx.rt.qm))
    ctx.log = (logged := []).append
    state = act(ctx, seeded(ctx, state, attempts=1), Action.REFINE_SADDLE)
    saddle = ctx.rt.saddle
    assert saddle.calls == ["refine"] * 2 and np.allclose(saddle.starts[1], saddle.starts[0] + 0.01)
    assert ctx.rt.screen_qm.calls.count("frequencies") == xtb + 1  # the fresh seed's only
    assert ctx.rt.qm.calls.count("frequencies") == dft + 1  # at the last frame
    assert saddle.guesses[1] is not None and saddle.guesses[1].task == "freq"
    assert [r["note"] for r in logged] == ["saddle_hessian:xtb:rho",
                                           "saddle:geometry_maxiter:maxiter",
                                           "saddle_hessian:seed_freq:mode",
                                           "saddle:geometry_maxiter:maxiter"]
    assert (state.seeds, state.attempts, state.saddle_pending) == ((), 2, False)
    assert decide(ctx.case, state, ctx.rules) == EXHAUSTED


@pytest.mark.parametrize("above_kcal,profile,continues", [
    (1.5, "string", False),  # climbed past the profile's maximum
    (0.5, "string", True),  # within a resolution of it
    (None, "string", True),  # a stored stalled search without its energy: no bound
    (1.5, None, True),  # no DFT profile yet (a low-level TS or ts_calc): no bound
])
def test_a_stalled_search_above_the_latest_profile_is_not_continued(tmp_path, above_kcal,
                                                                    profile, continues):
    """G1-P3: a continuous DFT path bounds the saddle from above, so a search that stalled more
    than a resolution above the latest profile's maximum has climbed past the barrier: no
    continuation, noted with the profile that bounds it. VAL7 s5: frames 29-32 kcal/mol above
    their HEI seed."""
    peak = 6.0 * K  # the profile's maximum, an interior node
    energy = None if above_kcal is None else peak + above_kcal * K
    pes = fakes.double_well()
    ctx, state = case_ctx(tmp_path, pes, saddle=ScriptedSaddle(tmp_path, pes, energy=energy))
    a, b = ctx.ends
    frames = [a + t * (b - a) for t in np.linspace(0.0, 1.0, 5)]
    ctx.work.path = None if profile is None else Profile(
        frames, tuple(e * K for e in (0, 3, 6, 2, 1)), profile)
    ctx.log = (logged := []).append
    state = act(ctx, seeded(ctx, state), Action.REFINE_SADDLE)
    assert len(ctx.rt.saddle.calls) == 1 + continues and state.attempts == 1
    bounded = {"note": "saddle:above_path_bound:string"} in logged
    assert bounded is (above_kcal == 1.5 and profile is not None)


# one failure token --------------------------------------------------------------------------

class SideFails(fakes.FakeQM):  # every QRC side optimization fails
    def optimize(self, mol, method, *, init_hessian=None, **kw):
        if init_hessian is not None:
            return Failure(kind=FailureKind.SCF_NOT_CONVERGED, reason="scf_unavailable:scf")
        return super().optimize(mol, method, **kw)


class ScfFails(fakes.FakeQM):  # every freq job fails its SCF
    def frequencies(self, mol, method, **kw):
        return Failure(kind=FailureKind.SCF_NOT_CONVERGED, reason="scf_unavailable:scf")


@pytest.mark.parametrize("token", ["connection", "ts_rejected", "unconverged", "path_engine",
                                   "scf_unavailable"])
def test_every_failure_is_one_token_that_clears_claim_and_connection(tmp_path, token):
    """U6-P4 (R2 never reached the old connection_failed row): a connection failure, a TS
    rejection, a search that does not converge, a failed string and an unavailable SCF each
    leave neither claim nor connection and spend one attempt; the case then goes on to its next
    seed or path while the budget lasts, and ends unresolved when it is spent."""
    pes = fakes.double_well()
    ctx, state = case_ctx(tmp_path, pes)
    ctx.log = (logged := []).append
    if token == "path_engine":
        ctx.rt = replace(ctx.rt, path=fakes.FakePath(tmp_path, pes, ["failed"]))
        state = act(ctx, state, Action.FIND_PATH)
        assert (state.attempts, state.path_runs) == (1, 1)
        assert decide(ctx.case, state, ctx.rules) == Decision(Action.FIND_PATH, "dft_path")
        return
    state = act(ctx, state, Action.FIND_PATH)  # a single peak: its seed
    ctx.rt = replace(ctx.rt, **{
        "connection": {"qm": SideFails(tmp_path, pes)},
        "ts_rejected": {"qm": OffAxisQM(tmp_path, pes)},
        "unconverged": {"saddle": FailingSaddle(tmp_path, pes)},
        "scf_unavailable": {"qm": ScfFails(tmp_path, pes)}}[token])
    state = act(ctx, replace(state, claim="earlier", connection="elementary"), Action.REFINE_SADDLE)
    if state.saddle_pending:
        state = act(ctx, state, Action.VALIDATE_AND_CONNECT)
    assert (state.claim, state.connection, state.attempts) == (None, None, 1)
    assert decide(ctx.case, state, ctx.rules) == Decision(Action.FIND_PATH, "dft_path")
    assert decide(ctx.case, replace(state, attempts=2), ctx.rules) == EXHAUSTED
    assert ctx.work.connection is None and logged


def test_a_string_that_left_no_seed_ends_the_case(tmp_path) -> None:
    """G3-P2: a string whose well cannot be relaxed and whose profile rises monotonically from
    it to the product leaves no seed; no next chunk runs, the case ends."""
    ctx, state = case_ctx(tmp_path, fakes.double_well())
    ctx.rt = replace(ctx.rt, qm=NoOpt(tmp_path, ctx.rt.qm.pes))
    a, b = ctx.ends
    frames = [a + t * (b - a) for t in np.linspace(0.0, 1.0, 5)]
    ctx.work.path = Profile(frames, tuple(e * K for e in (0, -3, 1, 2, 4)), "string")
    state = record_profile(state, BarrierVerdict(verdict="intermediate", source="string"))
    state = act(ctx, state, Action.VALIDATE_INTERMEDIATE)
    assert (state.intermediate, state.seeds, state.path_runs) == ("relax_failed", (), 1)
    assert decide(ctx.case, state, ctx.rules) == EXHAUSTED
    tried = replace(state, attempts=1)  # had its seed been tried, the next chunk runs
    assert decide(ctx.case, tried, ctx.rules) == Decision(Action.FIND_PATH, "dft_path")
    assert CaseOutcome.UNRESOLVED is EXHAUSTED.outcome

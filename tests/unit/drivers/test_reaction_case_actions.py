"""Reaction-case actions on fake surfaces (§7.3); K cases of the old path/decision tests."""

from dataclasses import replace
from pathlib import Path

import fakes
import numpy as np
import pytest

from hfauto.chemistry.xyz import composition_key
from hfauto.core.evidence import Failure, FailureKind
from hfauto.core.method import Deadline, MethodSpec
from hfauto.core.records import ReactionRecord, SpeciesRecord
from hfauto.drivers.minimum import Registry, relax_to_minimum
from hfauto.drivers.reaction_case.actions import HANDLERS
from hfauto.drivers.reaction_case.driver import CaseRuntime, open_case
from hfauto.drivers.reaction_case.state import Action, CasePolicy, CaseState, Decision, Seed, decide

DFT = MethodSpec(id="pbe0", kind="dft", functional="pbe0", basis="def2-svp")
XTB = MethodSpec(id="gfn2", kind="xtb", gfn=2)


class HigherOrderQM(fakes.FakeQM):  # every freq job also reports a second imaginary mode
    def frequencies(self, mol, method, *, deadline=None):
        ev = super().frequencies(mol, method, deadline=deadline)
        nu = sorted(ev.frequencies_cm1)
        return ev.model_copy(update={"frequencies_cm1": (nu[0], -200.0, *nu[2:]),
                                     "imaginary_modes": (*ev.imaginary_modes, tuple(np.eye(9)[4]))})


class FailingSaddle(fakes.FakeSaddle):
    def refine(self, seed, method, **kw):
        return Failure(kind=FailureKind.GEOMETRY_MAXITER, reason="maxiter")


def case_ctx(root: Path, pes, *, script=(), saddle=None, screen_pes=None):
    qm, load = fakes.FakeQM(root, pes), fakes.xyz_loader(root)
    registry, minima, species = Registry([], load), {}, {}
    for name in ("reactant", "product"):
        geo = fakes.write_geometry(root, f"in/{name}.xyz", pes.symbols, pes.points[name])
        species[name] = SpeciesRecord(species_id=name, state_label="x", charge=0,
                                      multiplicity=1, geometry=geo, source="input",
                                      composition_id=composition_key(pes.symbols, 0, 1))
        out = relax_to_minimum(pes.molecule(name), DFT, qm, load_xyz=load)
        record, _ = registry.add(out, species[name], tier="dft")
        minima[record.minimum_id] = (record, out.opt.final)
    rt = CaseRuntime(
        qm=qm, saddle=saddle or fakes.FakeSaddle(root, pes), path=fakes.FakePath(root, pes, script),
        screen_qm=fakes.FakeQM(root, screen_pes or pes), screen_path=fakes.FakePath(root, pes),
        method=DFT, screen_method=XTB, registry=registry, load_xyz=load, case_dir=root / "cases",
        file_ref=lambda p: fakes._ref(root, p), deadline=Deadline.after,
        resolve=lambda ref: root / ref.path, minima=minima, species=species)
    case = ReactionRecord(reaction_id="rx", reactants=(), products=(), minima=tuple(minima),
                          endpoints=("reactant", "product"), source="declared")
    ctx = open_case(case, rt, CasePolicy(screen=False), Deadline.after(600), root, lambda _: None)
    return ctx, CaseState(minima=tuple(record for record, _ in minima.values()))


def act(ctx, state, action, reason="test"):
    return HANDLERS[action](ctx, state, Decision(action, reason))


def test_monotonic_string_is_confirmed_by_a_second_discretization(tmp_path) -> None:
    ctx, state = case_ctx(tmp_path, fakes.flat_uphill())
    state = act(ctx, state, Action.FIND_PATH, "no_dft_path")
    assert decide(ctx.case, state, ctx.policy) == Decision(Action.FIND_PATH, "confirm_monotonic")
    state = act(ctx, state, Action.FIND_PATH, "confirm_monotonic")
    assert len(ctx.work.path[0]) == ctx.policy.confirm_barrierless_beads == 15
    assert decide(ctx.case, state, ctx.policy).reason == "dft_path_monotonic"


@pytest.mark.parametrize("pes,runs", [(fakes.flat_uphill, ("failed",)),
                                      (fakes.double_well, ("single_max",))])
def test_unconverged_string_runs_in_chunks_and_its_maximum_is_only_a_seed(tmp_path, pes, runs):
    ctx, state = case_ctx(tmp_path, pes(), script=["unconverged"] * 3)
    state = act(ctx, state, Action.FIND_PATH)
    assert state.path_runs == runs and ctx.rt.path.calls == ["find_path:unconverged"] * 3
    assert [s.source for s in state.seeds] == (["path_hei"] if runs == ("single_max",) else [])


def test_saddle_hessians_validation_and_qrc_on_a_double_well(tmp_path) -> None:
    ctx, state = case_ctx(tmp_path, fakes.double_well())
    state = act(ctx, state, Action.FIND_PATH)
    seeds, n_freq = state.seeds, ctx.rt.qm.calls.count("frequencies")
    state = act(ctx, state, Action.REFINE_SADDLE)  # 1st: xTB Hessian along the tangent
    assert state.last_saddle == "converged" and ctx.rt.screen_qm.calls == ["frequencies"]
    assert ctx.rt.qm.calls.count("frequencies") == n_freq and ctx.rt.saddle.calls == ["refine:0"]
    state = act(ctx, replace(state, seeds=seeds), Action.REFINE_SADDLE)  # 2nd: DFT Hessian
    assert ctx.rt.qm.calls.count("frequencies") == n_freq + 1 and state.saddle_attempts == 2
    state = act(ctx, state, Action.VALIDATE_TS)  # a separate freq job on the saddle
    assert state.ts_check == "ok" and state.claim.imag_cm1 < -50 and not state.claim.notes
    assert state.claim.freq_calc != state.claim.saddle_calc
    state = act(ctx, state, Action.CONNECT)
    claim = ctx.work.connection
    assert state.connection == "elementary" and set(claim.minima) == set(ctx.case.minima)
    state = act(ctx, replace(state, connection=None), Action.CONNECT)  # 2nd amplitude: × 2
    assert ctx.work.connection.amplitude_A == pytest.approx(2 * claim.amplitude_A)


def test_higher_order_retry_and_failed_saddle_routing(tmp_path) -> None:
    ctx, state = case_ctx(tmp_path, fakes.double_well())
    state = act(ctx, act(ctx, state, Action.FIND_PATH), Action.REFINE_SADDLE)
    ctx.rt = replace(ctx.rt, qm=HigherOrderQM(tmp_path, ctx.rt.qm.pes))
    state = act(ctx, state, Action.VALIDATE_TS)
    assert state.ts_check == "higher_order" and state.seeds[0].source == "higher_order_retry"
    assert decide(ctx.case, state, ctx.policy).action is Action.REFINE_SADDLE
    ctx, state = case_ctx(tmp_path / "b", fakes.double_well(), saddle=FailingSaddle(tmp_path, None))
    ts = fakes.write_geometry(tmp_path / "b", "ts.xyz", ("N", "H", "O"), ctx.rt.qm.pes.points["ts"])
    state = act(ctx, replace(state, seeds=(Seed(ts, "discovery_ts", ctx.chord(ctx.coords(ts))),)),
                Action.REFINE_SADDLE)  # a rejected seed routes to the DFT string
    assert decide(ctx.case, state, ctx.policy) == Decision(Action.FIND_PATH, "no_dft_path")


def test_screen_shortcut_and_collapsed_low_level_endpoints(tmp_path) -> None:
    ctx, state = case_ctx(tmp_path, fakes.double_well())
    ts = fakes.write_geometry(tmp_path, "ts.xyz", ("N", "H", "O"), ctx.rt.qm.pes.points["ts"])
    ctx.case = ctx.case.model_copy(update={"low_level_ts": ts})
    state = act(ctx, state, Action.SCREEN)
    assert state.screen.verdict == "proceed" and state.screen.n_dft_points == 3
    assert state.seeds[0].source == "discovery_ts" and ctx.rt.screen_path.calls == []
    # A single-well low-level surface collapses both endpoints: DFT IDPP, no seed.
    ctx, state = case_ctx(tmp_path / "b", fakes.double_well(), screen_pes=fakes.harmonic())
    state = act(ctx, state, Action.SCREEN)
    assert state.screen.verdict == "proceed" and not state.seeds and ctx.work.initial is not None
    assert ctx.rt.screen_path.calls == [] and state.screen.n_dft_points == 11

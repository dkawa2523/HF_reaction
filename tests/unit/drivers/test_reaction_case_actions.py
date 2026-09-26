"""Reaction-case actions on fake surfaces (§7.3); K cases of the old path/decision tests."""

from dataclasses import replace
from itertools import pairwise
from pathlib import Path

import fakes
import numpy as np
import pytest

from hfauto.chemistry.identity import mapped_rmsd
from hfauto.chemistry.interpolation import align_mapped
from hfauto.chemistry.xyz import XYZ, composition_key
from hfauto.chemistry.xyz_trajectory import read_xyz_trajectory, write_xyz_trajectory
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


class DriftingPath(fakes.FakePath):  # like NWChem freezeN: the last bead moves, image 3 rotates
    def __init__(self, root, pes, script):
        super().__init__(root, pes, script)
        self.initial = []

    def find_path(self, start, end, method, *, initial_path=None, **kw):
        self.initial.append(initial_path)
        run = super().find_path(start, end, method, initial_path=initial_path, **kw)
        x = [i.coords for i in read_xyz_trajectory(self.root / run.images.path)]
        x[-1][1] += [0.1, 0.0, 0.0]
        c, s = np.cos(np.pi / 6), np.sin(np.pi / 6)
        x[3] = x[3] @ np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]).T
        path = self.root / f"drift{len(self.initial)}.xyz"
        write_xyz_trajectory([XYZ(list(start.xyz.symbols), f) for f in x], path)
        return run.model_copy(update={"images": fakes._ref(self.root, path)})


def case_ctx(root: Path, pes, *, script=(), saddle=None, screen_pes=None):
    qm, load = fakes.FakeQM(root, pes), fakes.xyz_loader(root)
    registry, minima, species = Registry([], load), {}, {}
    for name in ("reactant", "product"):
        geo = fakes.write_geometry(root, f"in/{name}.xyz", pes.symbols, pes.points[name])
        species[name] = SpeciesRecord(species_id=name, state_label="x", charge=0,
                                      multiplicity=1, geometry=geo, source="input",
                                      composition_id=composition_key(pes.symbols, 0, 1))
        out = relax_to_minimum(pes.molecule(name), DFT, qm, load_xyz=load)
        record = registry.add(out, species[name], tier="dft")
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


def test_one_monotonic_string_with_settled_energies_is_barrierless(tmp_path) -> None:
    ctx, state = case_ctx(tmp_path, fakes.flat_uphill())  # gmax flat, bead energies settled
    state = act(ctx, state, Action.FIND_PATH, "no_dft_path")
    assert state.path_runs == ("monotonic",) and ctx.rt.path.calls == ["find_path:pes"]
    assert len(ctx.work.path[0]) == ctx.policy.string_beads and not state.seeds
    assert decide(ctx.case, state, ctx.policy).reason == "dft_path_monotonic"


@pytest.mark.parametrize("pes,runs", [(fakes.flat_uphill, ("monotonic",)),
                                      (fakes.double_well, ("single_max",))])
def test_unsettled_string_runs_in_chunks_and_its_maximum_is_only_a_seed(tmp_path, pes, runs):
    ctx, state = case_ctx(tmp_path, pes(), script=["unconverged"] * 3)
    state = act(ctx, state, Action.FIND_PATH)
    assert state.path_runs == runs and ctx.rt.path.calls == ["find_path:unconverged"] * 3
    assert [s.source for s in state.seeds] == (["path_hei"] if runs == ("single_max",) else [])


def test_string_chunks_restore_the_frozen_ends_and_align_the_images(tmp_path) -> None:
    ctx, state = case_ctx(tmp_path, fakes.double_well())
    ctx.rt = replace(ctx.rt, path=DriftingPath(tmp_path, ctx.rt.qm.pes, ["unconverged"]))
    logged = []
    ctx.log = logged.append
    state = act(ctx, state, Action.FIND_PATH)
    _, second = ctx.rt.path.initial  # the 2nd chunk starts from the true minima
    seam = ctx.frames(second)
    assert np.allclose(seam[0], ctx.ends[0], atol=1e-6) and mapped_rmsd(seam[-1], ctx.ends[1]) < 1e-6
    notes = [r["note"] for r in logged if r["note"].startswith("string0:c")]
    assert notes[0].startswith("string0:c1:end_drift:0.0") and float(notes[0][-5:]) > 0
    frames = ctx.work.path[0]  # the final path too, without the rigid rotation of image 3
    assert mapped_rmsd(frames[-1], ctx.ends[1]) < 1e-6 and state.path_runs == ("single_max",)
    assert all(np.allclose(align_mapped(a, b), b, atol=1e-6) for a, b in pairwise(frames))


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
    assert ctx.rt.qm.calls[-2:] == ["optimize+init_hessian"] * 2  # both sides: the TS Hessian
    state = act(ctx, replace(state, connection=None), Action.CONNECT)  # 2nd amplitude: × 2
    assert ctx.work.connection.amplitude_A == pytest.approx(2 * claim.amplitude_A)
    ctx.policy = replace(ctx.policy, qrc_bounds_A=(0.03, 0.05))  # the first one is capped
    state, amplitudes = replace(state, connection_attempts=0), []
    for _ in range(2):  # C19: × 2, then clipped, so the 2nd amplitude equals the 1st
        state = act(ctx, replace(state, connection=None), Action.CONNECT)
        amplitudes.append((state.connection, ctx.work.connection.amplitude_A))
    assert amplitudes == [("elementary", pytest.approx(0.05))] * 2


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

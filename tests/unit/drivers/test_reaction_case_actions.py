"""Reaction-case actions on fake surfaces (§7.3); K cases of the old path/decision tests."""

from dataclasses import replace
from itertools import pairwise
from pathlib import Path

import fakes
import numpy as np
import pytest

from hfauto.chemistry.geometry import declared_coordinate_gradient
from hfauto.chemistry.identity import mapped_rmsd
from hfauto.chemistry.interpolation import align_mapped
from hfauto.chemistry.modes import overlap
from hfauto.chemistry.xyz import XYZ, composition_key, read_xyz
from hfauto.chemistry.xyz_trajectory import read_xyz_trajectory, write_xyz_trajectory
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Failure, FailureKind
from hfauto.core.method import Deadline, MethodSpec
from hfauto.core.records import (
    BarrierVerdict,
    CoordinateTerm,
    ReactionRecord,
    SaddleClaim,
    SpeciesRecord,
)
from hfauto.drivers.minimum import Registry, relax_to_minimum
from hfauto.drivers.reaction_case.actions import HANDLERS, Profile
from hfauto.drivers.reaction_case.driver import CaseRuntime, open_case
from hfauto.drivers.reaction_case.state import (
    Action,
    CasePolicy,
    CaseState,
    Decision,
    Seed,
    decide,
    record_profile,
)

K = 1.0 / HARTREE_TO_KCAL_MOL
DFT = MethodSpec(id="pbe0", kind="dft", functional="pbe0", basis="def2-svp")
XTB = MethodSpec(id="gfn2", kind="xtb", gfn=2)
# Ill-conditioned ends for Kabsch: linear HCN -> HNC (C N H), planar cis -> trans HONO (H O N O)
LINEAR = ("CNH", [[0, 0, 0], [0, 0, 1.156], [0, 0, -1.066]],
          [[0, 0, 0], [0, 0, 1.17], [0, 0, 2.17]])
PLANAR = ("HONO", [[-0.324178, -0.094409, 0], [0.395107, -0.759978, 0], [1.561309, -0.078391, 0],
                   [1.429662, 1.092978, 0]],
          [[-0.197748, -0.943876, 0], [0.024911, -0.001454, 0], [1.407066, 0.008329, 0],
           [1.827671, 1.097200, 0]])
ACAC = Path(__file__).resolve().parents[3] / "configs" / "systems" / "xyz" / "acac"


class HigherOrderQM(fakes.FakeQM):  # every freq job also reports a second imaginary mode
    def frequencies(self, mol, method, *, deadline=None):
        ev = super().frequencies(mol, method, deadline=deadline)
        nu = sorted(ev.frequencies_cm1)
        return ev.model_copy(update={"frequencies_cm1": (nu[0], -200.0, *nu[2:]),
                                     "imaginary_modes": (*ev.imaginary_modes, tuple(np.eye(9)[4]))})


class FailingSaddle(fakes.FakeSaddle):  # a failure without a last frame
    def refine(self, seed, method, **kw):
        return Failure(kind=FailureKind.NONZERO_EXIT, reason="scripted")


class StalledSaddle(fakes.FakeSaddle):  # NWChem at maxiter: a Failure with the last frame
    def refine(self, seed, method, **kw):
        self.calls.append("refine")
        last = fakes.write_geometry(self.root, f"last{len(self.calls)}.xyz", seed.xyz.symbols,
                                    seed.xyz.coords + 0.01)
        return Failure(kind=FailureKind.GEOMETRY_MAXITER, reason="maxiter", final=last)


class DriftingPath(fakes.FakePath):  # like NWChem freezeN: the last bead moves, image 3 rotates
    def __init__(self, root, pes, script):
        super().__init__(root, pes, script)
        self.initial = []

    def find_path(self, start, end, method, *, initial_path, **kw):
        self.initial.append(initial_path)
        run = super().find_path(start, end, method, initial_path=initial_path, **kw)
        if isinstance(run, Failure):
            return run
        x = [i.coords for i in read_xyz_trajectory(self.root / run.images.path)]
        x[-1][1] += [0.1, 0.0, 0.0]
        c, s = np.cos(np.pi / 6), np.sin(np.pi / 6)
        x[3] = x[3] @ np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]).T
        path = self.root / f"drift{len(self.initial)}.xyz"
        write_xyz_trajectory([XYZ(list(start.xyz.symbols), f) for f in x], path)
        return run.model_copy(update={"images": fakes._ref(self.root, path)})


def case_ctx(root: Path, pes, *, script=(), saddle=None, screen_pes=None, neb=()):
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
        screen_qm=fakes.FakeQM(root, screen_pes or pes),
        screen_path=fakes.FakePath(root, pes, neb, tsopt=True),
        method=DFT, screen_method=XTB, registry=registry, load_xyz=load, case_dir=root / "cases",
        file_ref=lambda p: fakes._ref(root, p), deadline=Deadline.after,
        resolve=lambda ref: root / ref.path, minima=minima, species=species)
    case = ReactionRecord(reaction_id="rx", reactants=(), products=(), minima=tuple(minima),
                          endpoints=("reactant", "product"), source="declared")
    ctx = open_case(case, rt, CasePolicy(screen=False), Deadline.after(600), root, lambda _: None)
    return ctx, CaseState(minima=tuple(record for record, _ in minima.values()))


def act(ctx, state, action, reason="test"):
    return HANDLERS[action](ctx, state, Decision(action, reason))


def test_one_barrierless_string_closes_the_case(tmp_path) -> None:
    ctx, state = case_ctx(tmp_path, fakes.flat_uphill())
    state = act(ctx, state, Action.FIND_PATH, "no_dft_path")
    assert (state.screen.verdict, state.screen.source, state.path_runs) == (
        "barrierless", "string", 1)
    assert ctx.rt.path.calls == ["find_path:pes"] and not state.seeds
    assert len(ctx.work.path.frames) == ctx.policy.string_beads
    assert decide(ctx.case, state, ctx.policy).reason == "string:barrierless"


def test_one_chunk_per_call_and_the_next_one_continues_the_last_path(tmp_path) -> None:
    """U5-P4: a chunk is classified at once and its peak seeds the saddle search; once its seeds
    fail, the next chunk starts from its path, up to string_chunks."""
    ctx, state = case_ctx(tmp_path, fakes.double_well(), saddle=FailingSaddle(tmp_path, None))
    ctx.rt = replace(ctx.rt, path=DriftingPath(tmp_path, ctx.rt.qm.pes, ["single"] * 3))
    state = act(ctx, state, Action.FIND_PATH, "no_dft_path")
    first = ctx.work.path.frames
    assert ctx.rt.path.calls == ["find_path:single"] and state.seeds[0].source == "path_hei"
    for runs in (2, 3):
        state = act(ctx, state, Action.REFINE_SADDLE)  # its peak fails
        assert decide(ctx.case, state, ctx.policy) == Decision(Action.FIND_PATH, "next_chunk")
        state = act(ctx, state, Action.FIND_PATH, "next_chunk")
        assert state.path_runs == runs and len(ctx.rt.path.calls) == runs
    seam = ctx.frames(ctx.rt.path.initial[1])  # the 1st chunk's path, true minima at its ends
    assert all(np.allclose(x, y, atol=1e-4) for x, y in zip(seam, first, strict=True))
    assert decide(ctx.case, state, ctx.policy).reason == "attempts_exhausted"


def test_a_string_gets_the_dft_minima_back_and_its_images_aligned(tmp_path) -> None:
    ctx, state = case_ctx(tmp_path, fakes.double_well())
    ctx.rt = replace(ctx.rt, path=DriftingPath(tmp_path, ctx.rt.qm.pes, ["pes"]))
    state = act(ctx, state, Action.FIND_PATH)
    frames = ctx.work.path.frames  # without the moved end and the rigid rotation of image 3
    assert np.allclose(frames[0], ctx.ends[0]) and mapped_rmsd(frames[-1], ctx.ends[1]) < 1e-6
    assert all(np.allclose(align_mapped(a, b), b, atol=1e-6) for a, b in pairwise(frames))
    assert state.screen.verdict == "single"


@pytest.mark.parametrize("symbols,a,b", [LINEAR, PLANAR])
def test_first_string_chunk_starts_from_an_idpp_without_rigid_jumps(tmp_path, symbols, a, b):
    ctx, state = case_ctx(tmp_path, fakes.double_well())  # M9: NWChem rotates frozen bead N
    ctx.rt = replace(ctx.rt, path=DriftingPath(tmp_path, ctx.rt.qm.pes, ["failed"]))
    ctx.symbols, a = list(symbols), np.array(a, dtype=float)
    ctx.ends = (a, align_mapped(a, np.array(b, dtype=float)))
    act(ctx, state, Action.FIND_PATH)
    frames = ctx.frames(ctx.rt.path.initial[0])
    assert np.allclose(frames[0], a, atol=1e-6) and mapped_rmsd(frames[-1], ctx.ends[1]) < 1e-6
    rms = [np.sqrt(np.mean(np.sum((y - x) ** 2, axis=1))) for x, y in pairwise(frames)]
    assert rms == pytest.approx([mapped_rmsd(x, y) for x, y in pairwise(frames)], abs=1e-6)


def test_saddle_hessians_validation_and_qrc_on_a_double_well(tmp_path) -> None:
    ctx, state = case_ctx(tmp_path, fakes.double_well())
    state = act(ctx, state, Action.FIND_PATH)
    seeds, n_freq, logged = state.seeds, ctx.rt.qm.calls.count("frequencies"), []
    ctx.log = logged.append
    xtb = HigherOrderQM(tmp_path, ctx.rt.qm.pes)  # two negative modes do not matter
    ctx.rt = replace(ctx.rt, screen_qm=xtb)
    state = act(ctx, state, Action.REFINE_SADDLE)  # a negative xTB mode along the tangent
    assert state.last_saddle == "converged" and xtb.calls == ["frequencies"]
    assert ctx.rt.qm.calls.count("frequencies") == n_freq and ctx.rt.saddle.calls == ["refine"]
    ctx.rt = replace(ctx.rt, screen_qm=fakes.FakeQM(tmp_path, fakes.harmonic()))  # none
    state = act(ctx, replace(state, seeds=seeds), Action.REFINE_SADDLE)  # xTB first, then DFT
    assert ctx.rt.qm.calls.count("frequencies") == n_freq + 1 and state.saddle_attempts == 2
    notes = [r["note"].rsplit(":", 1)[0] for r in logged if "saddle_hessian" in r["note"]]
    assert notes == ["saddle_hessian:xtb:overlap", "saddle_hessian:dft:overlap"]
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


def test_a_higher_order_saddle_is_pushed_once_and_refined_from_its_ts_hessian(tmp_path):
    """U6-P3: the reaction mode is the negative mode along the direction; the most negative other
    one pushes the saddle once by the energy target, and the TS freq is the seed's Hessian."""
    ctx, state = case_ctx(tmp_path, fakes.double_well())
    state = act(ctx, act(ctx, state, Action.FIND_PATH), Action.REFINE_SADDLE)
    x = ctx.coords(ctx.work.saddle.final)
    ctx.rt = replace(ctx.rt, qm=HigherOrderQM(tmp_path, ctx.rt.qm.pes))
    state = act(ctx, state, Action.VALIDATE_TS)
    seed = state.seeds[0]
    assert (state.last_saddle, state.ts_check, seed.source) == ("failed", None,
                                                                "higher_order_retry")
    assert seed.hessian.task == "freq" and seed.tangent == seed.hessian.imaginary_modes[0]
    push = (ctx.coords(seed.geometry) - x).ravel()  # H along y: the second mode
    assert np.flatnonzero(np.abs(push) > 1e-9).tolist() == [4]
    assert ctx.policy.qrc_bounds_A[0] < abs(push[4]) < ctx.policy.qrc_bounds_A[1]
    decision = decide(ctx.case, state, ctx.policy)
    assert decision == Decision(Action.REFINE_SADDLE, "seed:higher_order_retry")
    logged, jobs = [], (ctx.rt.qm.calls.count("frequencies"), len(ctx.rt.screen_qm.calls))
    ctx.log = logged.append
    state = act(ctx, state, Action.REFINE_SADDLE)  # within 0.5 Å of the TS freq: no new Hessian
    assert state.last_saddle == "converged" and state.saddle_attempts == 2
    assert [r["note"] for r in logged] == ["saddle_hessian:ts_freq"]
    assert (ctx.rt.qm.calls.count("frequencies"), len(ctx.rt.screen_qm.calls)) == jobs


def test_a_stalled_saddle_restarts_once_from_its_last_frame(tmp_path) -> None:
    """U6-P6: maxiter queues the last frame as the next seed (same attempt budget), searched
    with a fresh Hessian; a failure without a last frame goes on to the string."""
    ctx, state = case_ctx(tmp_path, fakes.double_well(), saddle=StalledSaddle(tmp_path, None))
    ts = fakes.write_geometry(tmp_path, "ts.xyz", ("N", "H", "O"), ctx.rt.qm.pes.points["ts"])
    state = act(ctx, replace(state, seeds=(Seed(ts, "discovery_ts", None),)),
                Action.REFINE_SADDLE)
    restart = state.seeds[0]
    assert (restart.source, state.saddle_attempts, state.last_saddle) == (
        "saddle_restart", 1, "failed")
    assert np.allclose(ctx.coords(restart.geometry), ctx.rt.qm.pes.points["ts"] + 0.01)
    decision = decide(ctx.case, state, ctx.policy)
    assert decision == Decision(Action.REFINE_SADDLE, "seed:saddle_restart")
    xtb = ctx.rt.screen_qm.calls.count("frequencies")
    state = act(ctx, state, Action.REFINE_SADDLE)  # an xTB Hessian at the last frame
    assert ctx.rt.screen_qm.calls.count("frequencies") == xtb + 1
    assert not state.seeds and state.saddle_attempts == 2  # restarted once only
    assert decide(ctx.case, state, ctx.policy) == Decision(Action.FIND_PATH, "no_dft_path")
    root = tmp_path / "b"
    ctx, state = case_ctx(root, fakes.double_well(), saddle=FailingSaddle(root, None))
    ts = fakes.write_geometry(root, "ts.xyz", ("N", "H", "O"), ctx.rt.qm.pes.points["ts"])
    state = act(ctx, replace(state, seeds=(Seed(ts, "discovery_ts", None),)),
                Action.REFINE_SADDLE)
    assert not state.seeds
    assert decide(ctx.case, state, ctx.policy) == Decision(Action.FIND_PATH, "no_dft_path")


def test_a_collapsed_saddle_or_a_failed_soft_ts_is_validated_once(tmp_path) -> None:
    """U6-P6: VALIDATE_INTERMEDIATE consumes its trigger, so row 8 fires only for a new saddle
    (also after an earlier well was an endpoint)."""
    ctx, state = case_ctx(tmp_path, fakes.double_well())
    state = replace(act(ctx, state, Action.FIND_PATH), seeds=(), saddle_attempts=1)
    ctx.work.saddle = ctx.rt.qm.energy(ctx.mol(ctx.ends[0]), DFT)  # fell into the reactant
    soft = SaddleClaim(saddle_calc="s", freq_calc="f", imag_cm1=-30.0, energy_hartree=0.0)
    collapsed = Decision(Action.VALIDATE_INTERMEDIATE, "saddle_collapsed")
    for trigger in ({"ts_check": "collapsed"},
                    {"ts_check": "ok", "claim": soft, "connection": "failed",
                     "connection_attempts": 2}):
        s = replace(state, last_saddle="converged", intermediate="same_as_endpoint", **trigger)
        assert decide(ctx.case, s, ctx.policy) == collapsed
        s = act(ctx, s, Action.VALIDATE_INTERMEDIATE, "saddle_collapsed")
        assert (s.intermediate, s.ts_check, s.last_saddle, s.claim, s.connection) == (
            "same_as_endpoint", None, None, None, None)
        assert decide(ctx.case, s, ctx.policy) == Decision(Action.FIND_PATH, "next_chunk")


def test_reaction_direction_is_a_low_level_mode_the_reaction_centre_or_a_torsion(tmp_path):
    """U6-P2: acac PT, whose endpoint chord is dominated by the two methyl rotors."""
    ctx, _ = case_ctx(tmp_path, fakes.double_well())
    geo = ctx.rt.minima[ctx.case.minima[0]][1]  # direction() reads only the source and tangent
    mode = tuple(np.eye(9)[4])
    assert tuple(ctx.direction(ctx.ends[0], Seed(geo, "screen_ts", mode))) == mode
    r, p = (read_xyz(ACAC / f"{end}.xyz") for end in ("reactant", "product"))
    ctx.symbols, a = list(r.symbols), np.asarray(r.coords)
    ctx.ends = (a, align_mapped(a, np.asarray(p.coords)))
    chord, methyls = ctx.ends[1] - a, [7, 8, 9, 12, 13, 14]
    pt, rotors = np.zeros((15, 3)), np.zeros((15, 3))
    pt[10], rotors[methyls] = a[6] - a[2], chord[methyls]  # H10 from O2 to O6; methyl H
    direction = ctx.direction(a, Seed(geo, "higher_order_retry", None))
    assert overlap(chord, pt) < 0.3 and overlap(chord, rotors) > 0.9
    assert overlap(direction, pt) > 0.9 and overlap(direction, rotors) == 0.0
    ctx.symbols, a = list(PLANAR[0]), np.array(PLANAR[1])  # cis -> trans HONO: no bond change
    ctx.ends = (a, align_mapped(a, np.array(PLANAR[2])))
    torsion = ctx.direction(a, Seed(geo, "path_hei", tuple(np.ones(12))))
    assert np.allclose(torsion.reshape(4, 3)[:, :2], 0.0) and np.abs(torsion).max() > 0
    angle = (CoordinateTerm(kind="angle", atoms=(0, 1, 2)),)  # a declared coordinate first
    ctx.case = ctx.case.model_copy(update={"coordinate": angle})
    assert np.allclose(ctx.direction(a, Seed(geo, "path_hei", None)),
                       declared_coordinate_gradient(angle, a))


def test_screen_shortcut_from_a_low_level_ts(tmp_path) -> None:
    ctx, state = case_ctx(tmp_path, fakes.double_well())
    ts = fakes.write_geometry(tmp_path, "ts.xyz", ("N", "H", "O"), ctx.rt.qm.pes.points["ts"])
    ctx.case = ctx.case.model_copy(update={"low_level_ts": ts})
    energies = ctx.rt.qm.calls.count("energy")
    state = act(ctx, state, Action.SCREEN)
    assert state.screen.verdict == "single" and ctx.rt.qm.calls.count("energy") == energies + 1
    assert state.seeds[0].source == "discovery_ts" and ctx.rt.screen_path.calls == []


@pytest.mark.parametrize("xtb,neb,seed", [(None, (), "screen_ts"),
                                          (fakes.harmonic, (), "screen_hei"),
                                          (None, ("failed",), "screen_hei")])
def test_screen_takes_dft_energies_inside_the_neb_between_the_dft_minima(tmp_path, xtb, neb,
                                                                         seed):
    """U5-P2 / U5-P3: the xTB NEB from the IDPP between the DFT minima, DFT SPs on its interior
    only and no DFT freq; a single peak seeds at the NEB's TS when its xTB freq confirms it,
    else at the DFT peak; a failed NEB leaves the IDPP."""
    ctx, state = case_ctx(tmp_path, fakes.double_well(), screen_pes=xtb and xtb(), neb=neb)
    freq, energy = ctx.rt.qm.calls.count("frequencies"), ctx.rt.qm.calls.count("energy")
    logged = []
    ctx.log = logged.append
    state = act(ctx, state, Action.SCREEN)
    assert ctx.rt.screen_path.calls == [f"find_path:{(*neb, 'pes')[0]}"]
    assert ctx.rt.qm.calls.count("energy") == energy + ctx.policy.screen_images - 2
    assert ctx.rt.qm.calls.count("frequencies") == freq
    assert state.screen.verdict == "single" and [s.source for s in state.seeds] == [seed]
    frames = ctx.work.path.frames
    assert np.allclose(frames[0], ctx.ends[0]) and mapped_rmsd(frames[-1], ctx.ends[1]) < 1e-6
    assert ("screen_neb:nonzero_exit" in [r["note"] for r in logged]) is bool(neb)


def test_screen_closes_a_barrierless_case_between_the_dft_minima(tmp_path) -> None:
    ctx, state = case_ctx(tmp_path, fakes.flat_uphill())
    state = act(ctx, state, Action.SCREEN)
    assert state.screen.verdict == "barrierless" and state.screen.max_node_spacing_A > 0
    assert decide(ctx.case, state, ctx.policy).reason == "screen:barrierless"


def test_screen_finds_the_intermediate_of_a_two_step_path(tmp_path) -> None:
    ctx, state = case_ctx(tmp_path, fakes.triple_well())
    state = act(ctx, state, Action.SCREEN)
    assert state.screen.verdict == "intermediate" and not state.seeds
    decision = decide(ctx.case, state, ctx.policy)
    assert decision == Decision(Action.VALIDATE_INTERMEDIATE, "path_intermediate")
    state = HANDLERS[decision.action](ctx, state, decision)
    record, _ = ctx.work.intermediate
    assert state.intermediate == "distinct" and record.minimum_id not in ctx.case.minima
    assert decide(ctx.case, state, ctx.policy).reason == "intermediate_distinct"


@pytest.mark.parametrize("kcal,verdict,seeds", [
    ((0, 2, 0.5, 3, 6, 9, 10), "barrierless", []),  # nothing rises above the product
    ((0, 2, 0.5, 3, 6, 12, 10), "intermediate", ["path_hei"]),  # the highest peak
])
def test_a_well_that_is_an_endpoint_leaves_its_peak_once(tmp_path, kcal, verdict, seeds):
    ctx, state = case_ctx(tmp_path, fakes.double_well(), saddle=FailingSaddle(tmp_path, None))
    a, b = ctx.ends  # the well (node 2) relaxes back into the reactant
    frames = [a + t * (b - a) for t in np.linspace(0.0, 1.0, len(kcal))]
    ctx.work.path = Profile(frames, tuple(e * K for e in kcal))
    state = record_profile(state, BarrierVerdict(verdict="intermediate", source="string"))
    state = act(ctx, state, Action.VALIDATE_INTERMEDIATE, "path_intermediate")
    assert state.intermediate == "same_as_endpoint" and state.screen.verdict == verdict
    assert [s.source for s in state.seeds] == seeds
    if seeds:  # node 5, refined by the parabola a quarter step toward the product
        x = ctx.coords(state.seeds[0].geometry)
        assert np.allclose(x, frames[5] + 0.25 * (frames[6] - frames[5]), atol=1e-4)
        state = act(ctx, state, Action.REFINE_SADDLE)  # U6-P6: its failure adds no seed again
        assert decide(ctx.case, state, ctx.policy) == Decision(Action.FIND_PATH, "next_chunk")

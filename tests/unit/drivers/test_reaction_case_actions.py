"""Reaction-case path actions on fake surfaces (design §7.3); K cases of the old path/decision
tests. The saddle search, TS checks and QRC: test_reaction_case_saddle (on these builders)."""

from dataclasses import replace
from itertools import pairwise
from pathlib import Path

import fakes
import numpy as np
import pytest

from hfauto.chemistry import xyz
from hfauto.chemistry.geometry import declared_coordinate_gradient
from hfauto.chemistry.identity import mapped_rmsd
from hfauto.chemistry.interpolation import align_mapped
from hfauto.core import records as R
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Failure, FailureKind
from hfauto.core.method import Deadline, MethodSpec
from hfauto.drivers.minimum import Registry, relax_to_minimum
from hfauto.drivers.reaction_case.actions import Profile
from hfauto.drivers.reaction_case.driver import HANDLERS, CaseRuntime, open_case
from hfauto.drivers.reaction_case.paths import SCREEN_IMAGES, STRING_BEADS
from hfauto.drivers.reaction_case.state import (
    Action,
    CaseRules,
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
OFF_AXIS = tuple(np.eye(9)[4])  # A-H-B lies along x: H moves along y
EXHAUSTED = Decision(Action.COMPLETE, "attempts_exhausted", R.CaseOutcome.UNRESOLVED)
ACAC = Path(__file__).resolve().parents[3] / "configs" / "systems" / "xyz" / "acac"


class FailingSaddle(fakes.FakeSaddle):  # a failure without a last frame
    def refine(self, seed, method, **kw):
        return Failure(kind=FailureKind.NONZERO_EXIT, reason="scripted")


class NoOpt(fakes.FakeQM):  # every optimization fails at maxiter
    def optimize(self, mol, method, *, init_hessian=None, deadline=None):
        return Failure(kind=FailureKind.GEOMETRY_MAXITER, reason="scripted")


class Bumped(fakes.FakeQM):  # an SP within 0.01 Å of ``at`` rises by ``bump`` (None: fails)
    def __init__(self, root, pes, at, bump):
        super().__init__(root, pes)
        self.at, self.bump = at, bump

    def energy(self, mol, method, *, scf_guess=None, deadline=None):
        ev = super().energy(mol, method, scf_guess=scf_guess, deadline=deadline)
        if np.abs(mol.xyz.coords - self.at).max() > 0.01:
            return ev
        if self.bump is None:
            return Failure(kind=FailureKind.NONZERO_EXIT, reason="scripted")
        return ev.model_copy(update={"energy_hartree": ev.energy_hartree + self.bump})


class DriftingPath(fakes.FakePath):  # like NWChem freezeN: the last bead moves, image 3 rotates
    def __init__(self, root, pes, script):
        super().__init__(root, pes, script)
        self.initial = []

    def find_path(self, start, end, method, *, initial_path, **kw):
        self.initial.append(initial_path)
        run = super().find_path(start, end, method, initial_path=initial_path, **kw)
        if isinstance(run, Failure):
            return run
        x = [i.coords for i in xyz.read_xyz_trajectory(self.root / run.images.path)]
        x[-1][1] += [0.1, 0.0, 0.0]
        c, s = np.cos(np.pi / 6), np.sin(np.pi / 6)
        x[3] = x[3] @ np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]).T
        path = self.root / f"drift{len(self.initial)}.xyz"
        xyz.write_xyz_trajectory([xyz.XYZ(list(start.xyz.symbols), f) for f in x], path)
        return run.model_copy(update={"images": fakes._ref(self.root, path)})


class Batches(list):  # CaseRuntime.map: records each batch's size, runs it backwards (like
    def __call__(self, fn, items):  # threads finishing out of order), returns it in input order
        self.append(len(items))
        return [fn(x) for x in reversed(items)][::-1]


def case_ctx(root: Path, pes, *, script=(), saddle=None, screen_pes=None, neb=()):
    qm, load = fakes.FakeQM(root, pes), fakes.xyz_loader(root)
    registry, minima, species, ids = Registry([], load), {}, {}, []
    for name in ("reactant", "product"):
        geo = fakes.write_geometry(root, f"in/{name}.xyz", pes.symbols, pes.points[name])
        species[name] = R.SpeciesRecord(species_id=name, state_label="x", charge=0,
                                        multiplicity=1, geometry=geo, source="input",
                                        composition_id=xyz.composition_key(pes.symbols, 0, 1))
        out = relax_to_minimum(pes.molecule(name), DFT, qm, load_xyz=load)
        record = registry.add(out, species[name], tier="dft")
        minima.setdefault(record.minimum_id, (record, out.opt.final))  # image ends: one basin
        ids.append(record.minimum_id)
    rt = CaseRuntime(
        qm=qm, saddle=saddle or fakes.FakeSaddle(root, pes), path=fakes.FakePath(root, pes, script),
        screen_qm=fakes.FakeQM(root, screen_pes or pes),
        screen_path=fakes.FakePath(root, pes, neb, tsopt=True),
        method=DFT, screen_method=XTB, registry=registry, load_xyz=load, case_dir=root / "cases",
        file_ref=lambda p: fakes._ref(root, p), resolve=lambda ref: root / ref.path,
        map=Batches(), minima=minima, species=species, calcs={})
    case = R.ReactionRecord(reaction_id="rx", reactants=(), products=(), minima=tuple(ids),
                            endpoints=("reactant", "product"), degenerate=ids[0] == ids[1],
                            source="declared")
    ctx = open_case(case, rt, CaseRules(screen=False), Deadline.after(600), root, lambda _: None)
    return ctx, CaseState(minima=tuple(minima[m][0] for m in ids))


def act(ctx, state, action, reason="test"):
    return HANDLERS[action](ctx, state, Decision(action, reason))


def test_one_barrierless_string_closes_the_case(tmp_path) -> None:
    ctx, state = case_ctx(tmp_path, fakes.flat_uphill())
    state = act(ctx, state, Action.FIND_PATH, "dft_path")
    assert (state.screen.verdict, state.screen.source, state.path_runs) == (
        "barrierless", "string", 1)
    assert ctx.rt.path.calls == ["find_path:pes"] and not state.seeds
    assert len(ctx.work.path.frames) == STRING_BEADS + 2  # densified at 2 midpoints
    assert decide(ctx.case, state, ctx.rules).reason == "string:barrierless"


def test_one_chunk_per_call_and_the_next_one_continues_the_last_path(tmp_path) -> None:
    """U5-P4 / G3-P2: a chunk is classified at once and its peak seeds the saddle search; once
    its seed fails, the next chunk starts from its path while attempts are left."""
    ctx, state = case_ctx(tmp_path, fakes.double_well(), saddle=FailingSaddle(tmp_path, None))
    ctx.rt = replace(ctx.rt, path=DriftingPath(tmp_path, ctx.rt.qm.pes, ["single"] * 2))
    state = act(ctx, state, Action.FIND_PATH, "dft_path")
    first = ctx.work.path.frames
    assert ctx.rt.path.calls == ["find_path:single"] and state.seeds[0].source == "path_hei"
    state = act(ctx, state, Action.REFINE_SADDLE)  # its peak fails
    assert decide(ctx.case, state, ctx.rules) == Decision(Action.FIND_PATH, "dft_path")
    state = act(ctx, state, Action.FIND_PATH, "dft_path")
    assert state.path_runs == 2 and len(ctx.rt.path.calls) == 2
    seam = ctx.frames(ctx.rt.path.initial[1])  # the 1st chunk's path, true minima at its ends
    assert all(np.allclose(x, y, atol=1e-4) for x, y in zip(seam, first, strict=True))
    state = act(ctx, state, Action.REFINE_SADDLE)
    assert decide(ctx.case, state, ctx.rules) == EXHAUSTED  # the default budget: no 3rd chunk


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


def test_failed_shortcut_seeds_go_on_to_the_screen_path_once(tmp_path) -> None:
    """R4 / G8-P7: every low-level TS of the pair whose xTB freq confirms its mode gets one DFT
    SP (in one batch, no NEB) and a three-point verdict; the single-step ones seed in order.
    Once they have failed, SCREEN runs the xTB NEB and DFT SPs (not the shortcut again)."""
    ctx, state = case_ctx(tmp_path, fakes.double_well(), saddle=FailingSaddle(tmp_path, None))
    pes, ctx.log, ctx.rules = ctx.rt.qm.pes, (logged := []).append, CaseRules()
    near = pes.points["ts"] + [[0, 0, 0], [0.02, 0, 0], [0, 0, 0]]
    ts, well, near = (fakes.write_geometry(tmp_path, f"{n}.xyz", pes.symbols, x) for n, x in
                      (("ts", pes.points["ts"]), ("well", pes.points["reactant"]), ("near", near)))
    ctx.case = ctx.case.model_copy(update={"low_level_ts": (ts, well, near)})  # well: no TS mode
    energies = ctx.rt.qm.calls.count("energy")
    state = act(ctx, state, Action.SCREEN)
    assert [s.geometry for s in state.seeds] == [ts, near] and not state.neb_done
    assert {s.source for s in state.seeds} == {"discovery_ts"} and state.screen.verdict == "single"
    assert ctx.rt.qm.calls.count("energy") == energies + 2 and ctx.rt.map == [2]
    assert {"note": "low_level_ts_rejected"} in logged
    for attempts in (1, 2):
        state = act(ctx, state, Action.REFINE_SADDLE)
        assert state.saddle_attempts == attempts
    assert decide(ctx.case, state, ctx.rules) == Decision(Action.SCREEN, "screen")
    state = act(ctx, state, Action.SCREEN)
    assert state.neb_done and ctx.rt.screen_path.calls == ["find_path:pes"]
    assert ctx.rt.qm.calls.count("energy") == energies + 2 + SCREEN_IMAGES - 2
    assert decide(ctx.case, state, ctx.rules) == EXHAUSTED


def test_reaction_direction_is_a_ts_mode_rho_or_a_coordinate(tmp_path):
    """G1-P1: a TS seed's own mode, else ρ = ∇(Σ_broken r − Σ_formed r) of the labelled ends (acac
    PT); with no bond change the chord (NH3 inversion), the dihedral or the declared coordinate."""
    ctx, _ = case_ctx(tmp_path, fakes.double_well())
    geo = ctx.rt.minima[ctx.case.minima[0]][1]  # direction() reads only the seed's mode
    kind, mode = ctx.direction(ctx.ends[0], Seed(geo, "screen_ts", OFF_AXIS))
    assert kind == "mode" and tuple(mode) == OFF_AXIS
    r, p = (xyz.read_xyz(ACAC / f"{end}.xyz") for end in ("reactant", "product"))
    ctx.symbols, a = list(r.symbols), np.asarray(r.coords)
    ctx.ends = (a, align_mapped(a, np.asarray(p.coords)))
    kind, rho = ctx.direction(a, Seed(geo, "path_hei"))  # H10 from O2 to O6, no methyl rotor
    e2, e6 = ((a[10] - a[i]) / np.linalg.norm(a[10] - a[i]) for i in (2, 6))
    stretches = np.zeros((15, 3))
    stretches[[2, 6, 10]] = -e2, e6, e2 - e6  # ∇r(O2-H10) − ∇r(O6-H10)
    assert kind == "rho" and np.allclose(rho, stretches.ravel(), atol=1e-6)
    ctx.symbols, nh3 = list(fakes.NH3_HF_SYMBOLS[:4]), fakes.NH3_HF[:4]  # no bond, no dihedral
    ctx.ends = (nh3, align_mapped(nh3, -nh3))
    kind, umbrella = ctx.direction(nh3)
    assert kind == "chord" and np.allclose(umbrella, (ctx.ends[1] - nh3).ravel())
    ctx.symbols, a = list(PLANAR[0]), np.array(PLANAR[1])  # cis -> trans HONO: no bond change
    ctx.ends = (a, align_mapped(a, np.array(PLANAR[2])))
    kind, torsion = ctx.direction(a)
    assert kind == "coordinate" and np.allclose(torsion.reshape(4, 3)[:, :2], 0.0)
    assert np.abs(torsion).max() > 0
    angle = (R.CoordinateTerm(kind="angle", atoms=(0, 1, 2)),)  # a declared coordinate first
    ctx.case = ctx.case.model_copy(update={"coordinate": angle})
    kind, gradient = ctx.direction(a)
    assert kind == "coordinate" and np.allclose(gradient, declared_coordinate_gradient(angle, a))


@pytest.mark.parametrize("xtb,neb,seed", [(None, (), "screen_ts"),
                                          (fakes.harmonic, (), "screen_hei"),
                                          (None, ("failed",), "screen_hei")])
def test_screen_takes_dft_energies_inside_the_neb_between_the_dft_minima(tmp_path, xtb, neb, seed):
    """U5-P2 / U5-P3: the xTB NEB from the IDPP between the DFT minima, DFT SPs on its interior
    only and no DFT freq; a single peak seeds at the NEB's TS when its xTB freq confirms it,
    else at the DFT peak; a failed NEB leaves the IDPP."""
    ctx, state = case_ctx(tmp_path, fakes.double_well(), screen_pes=xtb and xtb(), neb=neb)
    freq, energy = ctx.rt.qm.calls.count("frequencies"), ctx.rt.qm.calls.count("energy")
    ctx.log = (logged := []).append
    state = act(ctx, state, Action.SCREEN)
    assert ctx.rt.screen_path.calls == [f"find_path:{(*neb, 'pes')[0]}"]
    assert ctx.rt.qm.calls.count("energy") == energy + SCREEN_IMAGES - 2
    assert ctx.rt.qm.calls.count("frequencies") == freq and ctx.rt.map == [SCREEN_IMAGES - 2]
    assert state.screen.verdict == "single" and [s.source for s in state.seeds] == [seed]
    frames = ctx.work.path.frames  # SP energies in frame order
    assert ctx.work.path.energies[1:-1] == tuple(ctx.rt.qm.pes.energy(f) for f in frames[1:-1])
    assert np.allclose(frames[0], ctx.ends[0]) and mapped_rmsd(frames[-1], ctx.ends[1]) < 1e-6
    assert ("screen_neb:nonzero_exit" in [r["note"] for r in logged]) is bool(neb)


def test_screen_closes_a_barrierless_case_between_the_dft_minima(tmp_path) -> None:
    ctx, state = case_ctx(tmp_path, fakes.flat_uphill())
    energies, ctx.log = ctx.rt.qm.calls.count("energy"), (logged := []).append
    state = act(ctx, state, Action.SCREEN)
    assert state.screen.verdict == "barrierless" and {"note": "screen_midpoints:barrierless"} in logged
    assert ctx.rt.qm.calls.count("energy") == energies + SCREEN_IMAGES - 2 + 2  # 2 midpoints
    assert len(ctx.work.path.frames) == len(ctx.work.path.energies) == SCREEN_IMAGES + 2
    assert decide(ctx.case, state, ctx.rules).reason == "screen:barrierless"


@pytest.mark.parametrize("bump,verdict", [(0.0, "barrierless"), (3 * K, "single"),
                                          (None, "unavailable")])
def test_a_barrierless_profile_is_densified_beside_its_highest_node(tmp_path, bump, verdict):
    """R9a: DFT SPs at the midpoints of the two segments beside the highest interior node; a
    barrier the nodes stepped over makes it a single step, a failed SP leaves no verdict."""
    ctx, _ = case_ctx(tmp_path, fakes.flat_uphill())
    a, b = ctx.ends
    frames = [a + t * (b - a) for t in np.linspace(0.0, 1.0, 7)]
    inner = [ev.energy_hartree for ev in ctx.sps(frames[1:-1])]
    k = 1 + int(np.argmax(inner))
    mids = [0.5 * (frames[i] + frames[i + 1]) for i in (k - 1, k)]
    ctx.rt = replace(ctx.rt, qm=Bumped(tmp_path, ctx.rt.qm.pes, mids[1], bump))
    v = ctx.verdict(frames, inner, "string")
    assert (v.verdict, v.source) == (verdict, "string") and ctx.rt.qm.calls == ["energy"] * 2
    assert v.reasons == (("midpoint_single_point",) if bump is None else ()) and ctx.rt.map[-1] == 2
    path = ctx.work.path
    assert len(path.frames) == 7 + 2 * (bump is not None)
    if bump is not None:
        assert np.allclose(path.frames[k], mids[0]) and np.allclose(path.frames[k + 2], mids[1])


def test_screen_finds_the_intermediate_of_a_two_step_path(tmp_path) -> None:
    ctx, state = case_ctx(tmp_path, fakes.triple_well())
    state = act(ctx, state, Action.SCREEN)
    assert state.screen.verdict == "intermediate" and not state.seeds
    decision = decide(ctx.case, state, ctx.rules)
    assert decision == Decision(Action.VALIDATE_INTERMEDIATE, "path_intermediate")
    state = HANDLERS[decision.action](ctx, state, decision)
    record, _ = ctx.work.intermediate
    assert state.intermediate == "distinct" and record.minimum_id not in ctx.case.minima
    assert decide(ctx.case, state, ctx.rules).reason == "intermediate_distinct"


@pytest.mark.parametrize("kcal,verdict,seeds", [
    ((0, 2, 0.5, 3, 6, 9, 10), "barrierless", []),  # nothing rises above the product
    ((0, 2, 0.5, 3, 6, 12, 10), "intermediate", ["path_hei"]),  # the highest peak
])
def test_a_well_that_is_an_endpoint_leaves_its_peak_once(tmp_path, kcal, verdict, seeds):
    ctx, state = case_ctx(tmp_path, fakes.double_well(), saddle=FailingSaddle(tmp_path, None))
    a, b = ctx.ends  # the well (node 2) relaxes back into the reactant
    frames = [a + t * (b - a) for t in np.linspace(0.0, 1.0, len(kcal))]
    ctx.work.path = Profile(frames, tuple(e * K for e in kcal), "string")
    state = record_profile(state, R.BarrierVerdict(verdict="intermediate", source="string"))
    state = act(ctx, state, Action.VALIDATE_INTERMEDIATE, "path_intermediate")
    assert state.intermediate == "same_as_endpoint" and state.screen.verdict == verdict
    assert [s.source for s in state.seeds] == seeds
    if seeds:  # node 5, refined by the parabola a quarter step toward the product
        x = ctx.coords(state.seeds[0].geometry)
        assert np.allclose(x, frames[5] + 0.25 * (frames[6] - frames[5]), atol=1e-4)
        state = act(ctx, state, Action.REFINE_SADDLE)  # U6-P6: its failure adds no seed again
        assert decide(ctx.case, state, ctx.rules) == Decision(Action.FIND_PATH, "dft_path")


def test_a_failed_well_relaxation_is_no_result_and_seeds_the_peak(tmp_path) -> None:
    """R9b: the well (node 2) cannot be relaxed: no endpoint and no barrierless claim; the
    highest peak (node 5) seeds the saddle search."""
    ctx, state = case_ctx(tmp_path, fakes.double_well())
    ctx.rt, ctx.log = replace(ctx.rt, qm=NoOpt(tmp_path, ctx.rt.qm.pes)), (logged := []).append
    a, b = ctx.ends
    frames = [a + t * (b - a) for t in np.linspace(0.0, 1.0, 7)]
    ctx.work.path = Profile(frames, tuple(e * K for e in (0, 2, 0.5, 3, 6, 9, 10)), "string")
    state = record_profile(state, R.BarrierVerdict(verdict="intermediate", source="string"))
    state = act(ctx, state, Action.VALIDATE_INTERMEDIATE, "path_intermediate")
    assert (state.intermediate, state.screen.verdict) == ("relax_failed", "intermediate")
    assert [s.source for s in state.seeds] == ["path_hei"] and {"note": "int:relax_failed"} in logged
    assert decide(ctx.case, state, ctx.rules) == Decision(Action.REFINE_SADDLE, "seed:path_hei")

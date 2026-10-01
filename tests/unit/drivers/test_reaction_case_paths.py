"""SCREEN's relaxed scan of an association (design X4, G2-P4) on a fake constrained optimizer:
the point order and count, each point from the previous optimum and its SCF, the profile from
the separated monomers' energy sum, its verdicts and a failed point."""

from pathlib import Path

import fakes
import numpy as np
import pytest

from hfauto.core import records as R
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Failure, FailureKind
from hfauto.core.method import Deadline, MethodSpec
from hfauto.drivers.minimum import Registry
from hfauto.drivers.reaction_case.driver import HANDLERS, CaseRuntime, open_case
from hfauto.drivers.reaction_case.paths import SCAN_POINTS, SCAN_REACH_A
from hfauto.drivers.reaction_case.state import Action, CaseRules, CaseState, Decision, decide

K = 1.0 / HARTREE_TO_KCAL_MOL
DFT = MethodSpec(id="pbe0", kind="dft", functional="pbe0", basis="def2-svp")
SCREEN = Decision(Action.SCREEN, "screen")
# H0 + O1-H2 -> H2O: the formed bond (0, 1); the OH fragment (holding j = 1) moves
SYMBOLS = ("H", "O", "H")
WATER = np.array([[-0.24, 0.93, 0.0], [0.0, 0.0, 0.0], [0.96, 0.0, 0.0]])
R_P = float(np.linalg.norm(WATER[1] - WATER[0]))
APART = WATER + np.outer([0, 1, 1], 1.5 * (WATER[1] - WATER[0]) / R_P)  # the complex, O-H0 2.46 A
E_H, E_OH, D_E, MORSE_A = -0.5, -75.7, 0.05, 3.0  # Eh, Eh, Eh, 1/A
TARGETS = [R_P + SCAN_REACH_A * (1 - k / (SCAN_POINTS - 1)) for k in range(SCAN_POINTS - 1)]


def energy(x: np.ndarray, bump_kcal: float) -> float:
    """A Morse well along r(O1-H0) below the separated monomers, with a Gaussian bump at the
    third scan point."""
    a = np.reshape(x, (-1, 3))
    r = float(np.linalg.norm(a[0] - a[1]))
    morse = D_E * ((1.0 - np.exp(-MORSE_A * (r - R_P))) ** 2 - 1.0)
    return E_H + E_OH + morse + bump_kcal * K * float(np.exp(-0.5 * ((r - TARGETS[2]) / 0.15) ** 2))


class ScanQM(fakes.FakeQM):
    """A constrained optimization at its fixed bond: H2 settles 0.01 A along z at every point
    (the free coordinates relax); ``fail`` scripts the point that fails."""

    def __init__(self, root, pes, fail=None):
        super().__init__(root, pes)
        self.points, self.results, self.fail = [], [], fail

    def optimize(self, mol, method, *, init_hessian=None, fixed_bond=None, scf_guess=None,
                 deadline=None):
        if fixed_bond is None:
            return super().optimize(mol, method, init_hessian=init_hessian, deadline=deadline)
        self.points.append((fixed_bond, scf_guess, np.array(mol.xyz.coords)))
        key = self._key("scan", mol.fingerprint(), fixed_bond)
        if len(self.points) - 1 == self.fail:
            return Failure(kind=FailureKind.GEOMETRY_MAXITER, reason="scripted", job_key=key)
        settled = mol.xyz.coords + [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.01]]
        self.results.append(self._evidence("opt", mol, method, key, self._start(mol, key),
                                           settled))
        return self.results[-1]


def association(root: Path, bump_kcal: float = 0.0, fail: int | None = None,
                collapsed: bool = False):
    """The case H0·OH (complex) -> H2O (adduct) with the separated H and OH as its monomers;
    ``collapsed``: the complex relaxed into the adduct's basin (minima[0] is the adduct's)."""
    pes = fakes.PES(SYMBOLS, lambda x: energy(x, bump_kcal), {})
    qm, load, notes = ScanQM(root, pes, fail), fakes.xyz_loader(root), []
    minima, species = {}, {}
    for mid, symbols, x, e in (("m_complex", SYMBOLS, APART, energy(APART, bump_kcal)),
                               ("m_adduct", SYMBOLS, WATER, E_H + E_OH - D_E),
                               ("m_h", ("H",), [[0.0, 0, 0]], E_H),
                               ("m_oh", ("O", "H"), [[0.0, 0, 0], [0.97, 0, 0]], E_OH)):
        geo = fakes.write_geometry(root, f"in/{mid}.xyz", symbols, np.asarray(x, dtype=float))
        species[mid] = R.SpeciesRecord(species_id=mid, composition_id=mid, charge=0,
                                       multiplicity=1, geometry=geo, source="input",
                                       state_label=mid)
        minima[mid] = (R.MinimumRecord(minimum_id=mid, basin_id=mid, composition_id=mid,
                                       species_id=mid, tier="dft", level_key="L", opt_calc="o",
                                       freq_calc="f", energy_hartree=e, state_label=mid), geo)
    rt = CaseRuntime(
        qm=qm, saddle=fakes.FakeSaddle(root, pes), path=fakes.FakePath(root, pes),
        screen_qm=None, screen_path=None, method=DFT, screen_method=None,
        registry=Registry([], load), load_xyz=load, file_ref=lambda p: fakes._ref(root, p),
        case_dir=root / "cases", resolve=lambda ref: root / ref.path,
        map=lambda fn, items: [fn(x) for x in items], minima=minima, species=species, calcs={})
    first = "m_adduct" if collapsed else "m_complex"
    case = R.ReactionRecord(reaction_id="assoc", reactants=(), products=(), source="declared",
                            minima=(first, "m_adduct"), endpoints=("m_complex", "m_adduct"),
                            monomers=("m_h", "m_oh"))
    ctx = open_case(case, rt, CaseRules(screen=False), Deadline.after(600), root,
                    lambda entry: notes.append(entry.get("note")))
    return ctx, CaseState(minima=(minima[first][0], minima["m_adduct"][0])), notes


def test_the_scan_runs_from_the_separated_monomers_inwards_to_the_adduct(tmp_path) -> None:
    """X4: without a low-level engine SCREEN scans the formed O1-H0 bond from r_P + 1.5 A in
    SCAN_POINTS - 1 constrained points, the first moved out of the adduct (H0 stays, the OH
    fragment translates rigidly), each next from the previous optimum and its SCF. The profile
    starts at the separated monomers' energy sum, not the complex's, and ends at the adduct; a
    monotonic one is barrierless and completes the case."""
    ctx, state, notes = association(tmp_path)
    assert ctx.energies == (E_H + E_OH, E_H + E_OH - D_E)
    assert decide(ctx.case, state, ctx.rules) == SCREEN
    state = HANDLERS[Action.SCREEN](ctx, state, SCREEN)
    points, results = ctx.rt.qm.points, ctx.rt.qm.results
    assert len(points) == SCAN_POINTS - 1 == 7
    assert [p[0][:2] for p in points] == [(0, 1)] * 7
    assert [p[0][2] for p in points] == pytest.approx(TARGETS)
    assert TARGETS[0] == pytest.approx(R_P + 1.5) and TARGETS[-1] > R_P
    starts = [p[2] for p in points]
    assert [np.linalg.norm(x[1] - x[0]) for x in starts] == pytest.approx(TARGETS)
    first = starts[0]  # the adduct with its OH translated rigidly, atom i in place
    assert np.linalg.norm(first[2] - first[1]) == pytest.approx(0.96)
    assert np.allclose(first[0], ctx.ends[1][0])
    assert [x[2, 2] for x in starts] == pytest.approx([0.01 * k for k in range(7)])
    assert [p[1] for p in points] == [None, *results[:-1]]  # the previous point's SCF
    path = ctx.work.path
    assert len(path.frames) == len(path.energies) == SCAN_POINTS + 1
    assert path.energies == (E_H + E_OH, *(ev.energy_hartree for ev in results),
                             E_H + E_OH - D_E)
    assert np.allclose(path.frames[0], path.frames[1])  # the longest point stands for the monomers
    assert state.screen == R.BarrierVerdict(verdict="barrierless", source="scan")
    assert state.neb_done and not state.seeds and state.path_runs == 0
    assert f"scan:0-1:{R_P:.3f}+1.5A:8_points" in notes and "scan:barrierless:" in notes
    assert decide(ctx.case, state, ctx.rules) == Decision(
        Action.COMPLETE, "scan:barrierless", R.CaseOutcome.BARRIERLESS)


def test_a_maximum_on_the_scan_seeds_the_saddle_search(tmp_path) -> None:
    """A 4 kcal/mol bump at the third point (H + C2H4-like): one peak, whose parabola-refined
    structure seeds the saddle search; no string follows a failed seed."""
    ctx, state, _ = association(tmp_path, bump_kcal=4.0)
    state = HANDLERS[Action.SCREEN](ctx, state, SCREEN)
    assert state.screen.verdict == "single" and state.seeds[0].source == "path_hei"
    seed = ctx.coords(state.seeds[0].geometry)
    assert TARGETS[2] < np.linalg.norm(seed[1] - seed[0]) < TARGETS[1]
    assert decide(ctx.case, state, ctx.rules) == Decision(Action.REFINE_SADDLE, "seed:path_hei")
    failed = CaseState(minima=state.minima, screen=state.screen, neb_done=True,
                       saddle_attempts=1, last_saddle="failed")
    assert decide(ctx.case, failed, ctx.rules).reason == "attempts_exhausted"


def test_a_failed_scan_point_leaves_the_profile_unavailable(tmp_path) -> None:
    """No point is dropped: the fourth point's failure ends the scan without a profile."""
    ctx, state, notes = association(tmp_path, fail=3)
    state = HANDLERS[Action.SCREEN](ctx, state, SCREEN)
    assert state.screen == R.BarrierVerdict(verdict="unavailable", source="scan",
                                            reasons=("scan_point",))
    assert len(ctx.rt.qm.points) == 4 and ctx.work.path is None and not state.seeds
    assert "scan3:geometry_maxiter" in notes
    assert decide(ctx.case, state, ctx.rules) == Decision(
        Action.COMPLETE, "attempts_exhausted", R.CaseOutcome.UNRESOLVED)


def test_a_complex_relaxed_into_its_adduct_ends_at_its_own_structure(tmp_path) -> None:
    """X4: a complex without a DFT minimum of its own (BH3 + NH3) leaves one basin; the case
    ends at the complex's input structure, which gives the formed bond, and scans as before."""
    ctx, state, _ = association(tmp_path, collapsed=True)
    assert np.allclose(ctx.raw[0], APART) and ctx.energies[0] == E_H + E_OH
    assert decide(ctx.case, state, ctx.rules) == SCREEN
    state = HANDLERS[Action.SCREEN](ctx, state, SCREEN)
    points = ctx.rt.qm.points
    assert [p[0][:2] for p in points] == [(0, 1)] * 7
    assert [p[0][2] for p in points] == pytest.approx(TARGETS)
    assert state.screen == R.BarrierVerdict(verdict="barrierless", source="scan")

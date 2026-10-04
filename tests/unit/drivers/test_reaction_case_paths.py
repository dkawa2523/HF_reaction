"""SCREEN's relaxed scan of an association (design X4, G2-P4) on a fake constrained optimizer:
the point order and count, each point from the previous optimum and its SCF, the profile from
the separated monomers' energy sum, its verdicts and a failed point; the approximate spin
projection of a low-spin-coupled pair's points (G2-P3). An open shell's profile SPs on its DFT
minima's SCF branch and the branch-jump verdict on the r9 audit's real profiles (G2-P1)."""

from pathlib import Path

import fakes
import numpy as np
import pytest
from scipy.spatial.distance import pdist

from hfauto.chemistry.gates import Policy
from hfauto.chemistry.xyz import XYZ, Molecule
from hfauto.core import records as R
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Failure, FailureKind
from hfauto.core.method import MethodSpec
from hfauto.drivers.minimum import Registry
from hfauto.drivers.reaction_case import paths
from hfauto.drivers.reaction_case.actions import branch_jump
from hfauto.drivers.reaction_case.driver import HANDLERS, CaseRuntime, open_case
from hfauto.drivers.reaction_case.paths import SCAN_POINTS, SCAN_REACH_A, projected
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
# probe G2-1v, S5 CH3 + O2 doublet on its 11 IDPP frames (kcal/mol, <S2>): the frame-by-frame
# continuation, the atomic guess, and the reactant end's vectors for 1-3 with the adduct's after
CHAIN = ((0.0, 0.124, 0.448, 0.802, 0.779, -0.496, -4.683, -12.841, -22.612, -31.848, -36.503),
         (1.7115, 1.6785, 1.6261, 1.5439, 1.4145, 1.2107, 0.9167, 0.7682, 0.7565, 0.7548, 0.7545))
ATOMIC = ((0.0, 28.930, 23.772, 18.160, 11.967, 5.108, -4.683, -12.841, -22.612, -31.848, -36.503),
          (1.7115, 0.7591, 0.7605, 0.7631, 0.7696, 0.8114, 0.9166, 0.7682, 0.7565, 0.7548, 0.7545))
MIXED = ((0.0, 0.124, 0.448, 0.802, 11.967, 5.108, -4.683, -12.841, -22.612, -31.848, -36.503),
         (1.7115, 1.6785, 1.6261, 1.5438, 0.7697, 0.8111, 0.9167, 0.7682, 0.7565, 0.7548, 0.7545))
# the r9 S4 audit's open-shell DFT profiles from minimum to minimum (S4/audit/cur_open_shell.txt,
# val7_open_shell.txt; kcal/mol, <S2>): SCREEN's, and the shortcut's [minimum, TS SP, minimum]
S18_SCREEN = ((0.0, 0.26, 0.85, 2.28, 5.75, 7.56, 5.93, 2.38, 0.86, 0.26, 0.0),
              (0.75, 0.75, 0.7502, 0.7511, 0.7565, 0.7674, 0.7570, 0.7512, 0.7502, 0.75, 0.75))
S4_SCREEN = ((0.0, 4.73, 18.52, 34.31, 34.31, 34.52, 14.93, -0.29, -3.99, -6.40, -7.75),
             (0.7543, 0.7539, 0.7539, 0.7582, 0.7583, 0.7585, 0.7547, 0.7543, 0.7544, 0.7544,
              0.7543))
S6_SCREEN = ((0.0, 0.93, 0.85, 3.02, 10.78, 24.13, 39.13, 37.26, 6.29, -12.40, -13.40),
             (0.7526, 0.7526, 0.7526, 0.7526, 0.7527, 0.7541, 0.7728, 0.7881, 0.7552, 0.7540,
              0.7544))
S6_SHORTCUT = ((0.0, 56.33, 25.14), (0.7544, 0.7722, 0.75))  # f9bc3a1cd8
S5_SHORTCUT = ((0.0, 51.83, -17.14), (0.7545, 0.7594, 0.7526))  # b3e2f366b6
VAL7_S5_SCREEN = ((0.0, 26.23, 18.17, 8.79, -4.51, -20.12, -31.81, -35.55, -35.54, -35.83, 0.40),
                  (1.7114, 0.7604, 0.7651, 0.8255, 0.7750, 0.7554, 0.7548, 0.7547, 0.7547,
                   0.7546, 1.7637))  # b3febc8052: both minima broken-symmetry
TOL = Policy().spin_tol
# CH3O -> CH2OH (C O H H H; H2 moves from C to O), both C1; the product optimized in a frame
# turned 90 degrees about z and shifted
CH3O = np.array([[0.0, 0.0, 0.0], [1.37, 0.02, 0.01], [-0.37, 1.03, 0.05], [-0.33, -0.55, 0.88],
                 [-0.40, -0.47, -0.91]])
CH2OH = np.array([[0.0, 0.0, 0.0], [1.36, 0.05, -0.02], [1.75, 0.90, 0.10],
                  [-0.45, -0.40, 0.92], [-0.50, -0.35, -0.90]])
TURN = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])


def energy(x: np.ndarray, bump_kcal: float) -> float:
    """A Morse well along r(O1-H0) below the separated monomers, with a Gaussian bump at the
    third scan point."""
    a = np.reshape(x, (-1, 3))
    r = float(np.linalg.norm(a[0] - a[1]))
    morse = D_E * ((1.0 - np.exp(-MORSE_A * (r - R_P))) ** 2 - 1.0)
    return E_H + E_OH + morse + bump_kcal * K * float(np.exp(-0.5 * ((r - TARGETS[2]) / 0.15) ** 2))


def bs_s2(r: float) -> float:  # a broken-symmetry pair's <S2>: 0.75 bonded, ~1.72 apart
    return 0.75 + 1.0 - float(np.exp(-(((r - R_P) / 0.8) ** 2)))


class ScanQM(fakes.FakeQM):
    """A constrained optimization at its fixed bond: H2 settles 0.01 A along z at every point
    (the free coordinates relax); ``fail`` scripts the point that fails. An open shell's points
    are broken-symmetry (``bs_s2``); its SPs are the high-spin ones, repulsive, <S2> 3.76 but
    ``dirty`` (that point's: 3.95)."""

    def __init__(self, root, pes, fail=None, dirty=None):
        super().__init__(root, pes)
        self.points, self.results, self.fail, self.dirty, self.high = [], [], fail, dirty, []

    def energy(self, mol, method, *, scf_guess=None):
        x = np.array(mol.xyz.coords)
        self.high.append((mol.multiplicity, x))
        r = float(np.linalg.norm(x[0] - x[1]))
        s2 = 3.95 if len(self.high) - 1 == self.dirty else 3.76
        return super().energy(mol, method).model_copy(update={
            "energy_hartree": E_H + E_OH + 0.03 * np.exp(-2.0 * (r - R_P)), "s2": s2})

    def optimize(self, mol, method, *, init_hessian=None, fixed_bond=None, scf_guess=None):
        if fixed_bond is None:
            return super().optimize(mol, method, init_hessian=init_hessian)
        self.points.append((fixed_bond, scf_guess, np.array(mol.xyz.coords)))
        key = self._key("scan", mol.fingerprint(), fixed_bond)
        if len(self.points) - 1 == self.fail:
            return Failure(kind=FailureKind.GEOMETRY_MAXITER, reason="scripted", job_key=key)
        settled = mol.xyz.coords + [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.01]]
        s2 = {"s2": bs_s2(fixed_bond[2])} if mol.multiplicity > 1 else {}
        self.results.append(self._evidence("opt", mol, method, key, self._start(mol, key),
                                           settled, **s2))
        return self.results[-1]


def minimum(root: Path, mid: str, symbols, x, e: float, multiplicity: int = 1, opt: str = "o"):
    """(species, (minimum, its geometry)) of a DFT minimum ``mid`` at ``x``."""
    geo = fakes.write_geometry(root, f"in/{mid}.xyz", symbols, np.asarray(x, dtype=float))
    species = R.SpeciesRecord(species_id=mid, composition_id=mid, charge=0, state_label=mid,
                              multiplicity=multiplicity, geometry=geo, source="input")
    return species, (R.MinimumRecord(minimum_id=mid, basin_id=mid, composition_id=mid,
                                      species_id=mid, tier="dft", level_key="L", opt_calc=opt,
                                      freq_calc="f", energy_hartree=e, state_label=mid), geo)


def context(root: Path, qm, minima, species, case, calcs=None):
    """The case context on ``qm`` (no low-level engine, ``map`` in order) and its notes."""
    load, pes, notes = fakes.xyz_loader(root), qm.pes, []
    rt = CaseRuntime(
        qm=qm, saddle=fakes.FakeSaddle(root, pes), path=fakes.FakePath(root, pes),
        screen_qm=None, screen_path=None, method=DFT, screen_method=None,
        registry=Registry(minima.values(), load), load_xyz=load,
        file_ref=lambda p: fakes._ref(root, p),
        case_dir=root / "cases", resolve=lambda ref: root / ref.path,
        map=lambda fn, items: [fn(x) for x in items], species=species,
        calcs=calcs or {})
    return open_case(case, rt, CaseRules(screen=False), root,
                     lambda entry: notes.append(entry.get("note"))), notes


def association(root: Path, bump_kcal: float = 0.0, fail: int | None = None,
                collapsed: bool = False, spins=(1, 1, 1), dirty: int | None = None):
    """The case H0·OH (complex) -> H2O (adduct) with the separated H and OH as its monomers;
    ``collapsed``: the complex relaxed into the adduct's basin (minima[0] is the adduct's).
    ``spins``: the multiplicities of the pair, H and OH ((2, 2, 3) couples H and a triplet
    partner to a doublet, as CH3 + O2)."""
    pes = fakes.PES(SYMBOLS, lambda x: energy(x, bump_kcal), {})
    found = {mid: minimum(root, mid, symbols, x, e, m) for mid, symbols, x, e, m in (
        ("m_complex", SYMBOLS, APART, energy(APART, bump_kcal), spins[0]),
        ("m_adduct", SYMBOLS, WATER, E_H + E_OH - D_E, spins[0]),
        ("m_h", ("H",), [[0.0, 0, 0]], E_H, spins[1]),
        ("m_oh", ("O", "H"), [[0.0, 0, 0], [0.97, 0, 0]], E_OH, spins[2]))}
    minima = {mid: m for mid, (_, m) in found.items()}
    first = "m_adduct" if collapsed else "m_complex"
    case = R.ReactionRecord(reaction_id="assoc", reactants=(), products=(), source="declared",
                            minima=(first, "m_adduct"), endpoints=("m_complex", "m_adduct"),
                            monomers=("m_h", "m_oh"))
    ctx, notes = context(root, ScanQM(root, pes, fail, dirty), minima,
                       {mid: s for mid, (s, _) in found.items()}, case)
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


def test_the_approximate_spin_projection() -> None:
    """Yamaguchi's AP on probe G2-1's S5 scan at C-O 2.6 A: alpha = 0.4285 lowers the
    broken-symmetry doublet by 1.77 kcal/mol (BS -0.98, AP -2.75 against CH3 + O2); a solution
    of pure spin is left as it is."""
    e_bs, e_hs = -189.820102684952, -189.813539102947
    e_ap = projected(e_bs, 1.6542, e_hs, 3.7643, 0.5)
    assert (e_ap - e_bs) / K == pytest.approx(-1.765, abs=1e-3)
    assert projected(e_bs, 0.75, e_hs, 3.75, 0.5) == e_bs


def test_a_low_spin_coupled_scan_is_spin_projected_point_by_point(tmp_path) -> None:
    """G2-P3: H + a triplet partner as a doublet: a quartet SP at each point's structure, AP
    energies on the scan points only (the monomers' sum and the adduct keep theirs), the verdict
    on them; a point whose quartet is itself contaminated keeps its BS energy, noted."""
    ctx, state, notes = association(tmp_path, spins=(2, 2, 3), dirty=2)
    state = HANDLERS[Action.SCREEN](ctx, state, SCREEN)
    qm = ctx.rt.qm
    assert [m for m, _ in qm.high] == [4] * 7 and "scan:ap:4" in notes
    assert all(np.allclose(x, ctx.coords(bs.final)) for (_, x), bs in zip(qm.high, qm.results))
    sps = [ev for ev in ctx.work.calcs.values() if ev.task == "sp"]
    assert len(sps) == 7 and "scan2:ap_skipped:s2=3.95" in notes
    expected = [projected(bs.energy_hartree, bs.s2, hs.energy_hartree, hs.s2, 0.5)
                for bs, hs in zip(qm.results, sps, strict=True)]
    expected[2] = qm.results[2].energy_hartree
    path = ctx.work.path
    assert path.energies == pytest.approx((E_H + E_OH, *expected, E_H + E_OH - D_E))
    assert all(e < bs.energy_hartree for e, bs in zip(path.energies[1:3], qm.results))
    assert state.screen == R.BarrierVerdict(verdict="barrierless", source="scan")


@pytest.mark.parametrize("spins", [(4, 2, 3), (2, 2, 1), (1, 1, 1)])
def test_no_projection_without_a_low_spin_coupling(tmp_path, spins) -> None:
    """The quartet declared for H + a triplet is the high-spin coupling, a doublet of H and a
    closed shell has one open-shell monomer, a closed-shell pair none: no SP, BS energies."""
    ctx, state, notes = association(tmp_path, spins=spins)
    HANDLERS[Action.SCREEN](ctx, state, SCREEN)
    assert not ctx.rt.qm.high and not any(n and n.startswith("scan:ap") for n in notes)
    assert ctx.work.path.energies[1:-1] == tuple(ev.energy_hartree for ev in ctx.rt.qm.results)


@pytest.mark.parametrize("profile", [VAL7_S5_SCREEN, ATOMIC, MIXED])
def test_the_branch_jump_rule_flags_the_maxima_off_their_branch(profile) -> None:
    """G2-P1: VAL7 S5's SCREEN maximum (0.7604 beside 1.7114 and 0.7651), probe G2-1v's
    atomic-guess maximum (0.7591 beside 1.7115 and 0.7605) and its mixed-guess one at frame 4
    (0.7697 beside 1.5438 and 0.8111): outside the neighbours' <S2> and 0.77-0.95 from one."""
    energies, s2 = profile
    assert branch_jump(energies, s2, TOL)
    k = 1 + int(np.argmax(energies[1:-1]))
    assert not branch_jump(energies, (*s2[:k], None, *s2[k + 1:]), TOL)  # unobserved: no evidence


@pytest.mark.parametrize("profile", [S18_SCREEN, S4_SCREEN, S6_SCREEN, S6_SHORTCUT, S5_SHORTCUT,
                                     CHAIN])
def test_the_branch_jump_rule_keeps_the_smooth_radical_maxima(profile) -> None:
    """The UKS contamination peaking at a genuine radical TS lies outside its neighbours' <S2>
    but within spin_tol of both (S18 0.7674 beside 0.7565 and 0.7570, at most 0.022: S6's
    shortcut 0.7722 beside 0.7544 and 0.7500); S6's SCREEN maximum and the continuation's
    Coulson-Fischer hill lie between theirs. No CUR profile is flagged."""
    assert not branch_jump(*profile, TOL)


class Branches(fakes.FakeQM):
    """SCF branches on a flat surface: an SP keeps its guess's <S2> (from scratch 0.76) and lies
    ``offset[guess job key]`` Eh up, unless ``script`` gives its (Eh, <S2>); each SP's guess and
    structure are recorded."""

    def __init__(self, root, pes, offset, script=()):
        super().__init__(root, pes)
        self.offset, self.script, self.sps = offset, list(script), []

    def energy(self, mol, method, *, scf_guess=None):
        self.sps.append((scf_guess, np.array(mol.xyz.coords)))
        ev = super().energy(mol, method, scf_guess=scf_guess)
        e, s2 = self.script.pop(0) if self.script else (
            self.offset.get(scf_guess and scf_guess.job_key, 0.0),
            scf_guess.s2 if scf_guess else 0.76)
        return ev.model_copy(update={"energy_hartree": e, "s2": s2})


def radical(root: Path, s2_ends=(0.754, 0.754), multiplicity: int = 2, offset=None,
            e_ends=(0.0, 0.0), script=()):
    """CH3O -> CH2OH (a doublet) on a flat surface, its DFT minima at ``e_ends`` (Eh) with
    opts of <S2> ``s2_ends``; ``offset`` (opt index -> Eh) raises the SPs started from that
    minimum's opt. The context, the opts and the notes."""
    symbols = ("C", "O", "H", "H", "H")
    pes, calcs, found = fakes.PES(symbols, lambda x: 0.0, {}), {}, {}
    for k, (mid, x) in enumerate((("m_a", CH3O), ("m_b", CH2OH @ TURN.T + 1.0))):
        mol = Molecule(XYZ(list(symbols), x), 0, multiplicity)
        calcs[f"opt_{mid}"] = fakes.FakeQM(root, pes).energy(mol, DFT).model_copy(
            update={"task": "opt", "s2": s2_ends[k]})
        found[mid] = minimum(root, mid, symbols, x, e_ends[k], multiplicity, f"opt_{mid}")
    offsets = {calcs[f"opt_m_{'ab'[k]}"].job_key: e for k, e in (offset or {}).items()}
    case = R.ReactionRecord(reaction_id="shift", reactants=(), products=(), source="declared",
                            minima=("m_a", "m_b"), endpoints=("m_a", "m_b"))
    ctx, notes = context(root, Branches(root, pes, offsets, script),
                         {m: v for m, (_, v) in found.items()},
                         {m: s for m, (s, _) in found.items()}, case, calcs)
    return ctx, calcs, notes


def path(ctx, count=5):
    a, b = ctx.ends
    return [a + t * (b - a) for t in np.linspace(0.0, 1.0, count)]


def test_only_an_open_shell_starts_its_sps_from_the_nearer_minimum(tmp_path) -> None:
    """G2-P1: a closed shell's SPs start from scratch at the frames. An open shell's start from
    the opt of the nearer DFT minimum, the frame carried into that opt's frame (NWChem reads
    the vectors unrotated): unchanged near the reactant, turned back with the product near it;
    energy-invariant, as internal distances are kept. The densifying midpoints too."""
    closed, _, _ = radical(tmp_path / "closed", multiplicity=1)
    frames = path(closed, 6)
    closed.sps(frames[1:-1])
    assert [g for g, _ in closed.rt.qm.sps] == [None] * 4
    assert all(np.allclose(x, f) for (_, x), f in zip(closed.rt.qm.sps, frames[1:-1]))
    ctx, calcs, _ = radical(tmp_path / "open")
    a, b = calcs["opt_m_a"], calcs["opt_m_b"]
    found = ctx.sps(frames := path(ctx, 6))
    assert [g for g, _ in ctx.rt.qm.sps] == [a, a, a, b, b, b]
    assert [ev.s2 for ev in found] == [0.754] * 6
    (_, near_a), (_, near_b) = ctx.rt.qm.sps[1], ctx.rt.qm.sps[4]
    assert np.allclose(near_a, frames[1], atol=1e-8)
    assert pdist(near_b) == pytest.approx(pdist(frames[4]))
    rms = [float(np.sqrt(np.mean((near_b - end) ** 2))) for end in (ctx.raw[1], ctx.ends[1])]
    assert rms[0] < 0.2 and rms[1] > 0.5
    ctx.rt.qm.sps.clear()
    ctx.verdict(frames, [0.0] * 4, "screen", [0.754] * 4)  # flat: barrierless, densified
    assert len(ctx.rt.qm.sps) == 2 and all(g is not None for g, _ in ctx.rt.qm.sps)


def test_minima_of_different_spin_coupling_start_each_sp_from_both(tmp_path) -> None:
    """S5's broken-symmetry complex (1.71) and its adduct (0.754): every frame is solved from
    both minima's opts and keeps the lower solution (here the adduct's branch), every
    solution kept."""
    ctx, calcs, _ = radical(tmp_path, s2_ends=(1.71, 0.754), offset={0: 0.002})
    found = ctx.sps(path(ctx)[1:-1])
    a, b = calcs["opt_m_a"], calcs["opt_m_b"]
    assert [g for g, _ in ctx.rt.qm.sps] == [a, b] * 3
    assert [(ev.energy_hartree, ev.s2) for ev in found] == [(0.0, 0.754)] * 3
    assert len(ctx.work.calcs) == 6


def ends(profile):
    """A profile's minima: their energies (Eh) and <S2>."""
    (e, s2) = profile
    return {"e_ends": (e[0] * K, e[-1] * K), "s2_ends": (s2[0], s2[-1])}


@pytest.mark.parametrize("profile,verdict", [
    (VAL7_S5_SCREEN, "unavailable"), (S18_SCREEN, "single"), (S4_SCREEN, "single")])
def test_a_branch_jump_leaves_the_profile_without_a_class(tmp_path, profile, verdict) -> None:
    """SCREEN's profile: a maximum off its SCF branch is unavailable (scf_branch_jump, so no
    seed), every point kept, whatever the minima's spin coupling (VAL7 S5's both
    broken-symmetry); a radical TS's smooth <S2> peak is no jump."""
    ctx, _, _ = radical(tmp_path, **ends(profile))
    energies, s2 = profile
    v = ctx.verdict(path(ctx, 11), [e * K for e in energies[1:-1]], "screen", s2[1:-1])
    assert v.verdict == verdict and len(ctx.work.path.energies) == 11
    assert v.reasons == (("scf_branch_jump",) if verdict == "unavailable" else ())


def test_the_densified_profile_is_judged_with_its_midpoints_s2(tmp_path) -> None:
    """A barrierless profile's midpoints carry their <S2>: one off the branch beside the highest
    node (as the atomic guess's 0.7591 beside 1.7115) leaves the densified profile unavailable."""
    ctx, _, notes = radical(tmp_path, s2_ends=(1.7115, 1.65), e_ends=(0.0, -4.0 * K),
                            script=[(-0.2 * K, 0.7591), (-1.5 * K, 1.66)])
    v = ctx.verdict(path(ctx), [e * K for e in (-1.0, -2.0, -3.0)], "screen",
                    (1.6785, 1.6261, 1.5439))
    assert len(ctx.work.path.energies) == 7 and "screen_midpoints:unavailable" in notes
    assert v.reasons == ("scf_branch_jump",)


@pytest.mark.parametrize("profile,seeded", [
    (S5_SHORTCUT, True), (S6_SHORTCUT, True),
    (((0.0, 26.23, 0.40), (1.7114, 0.7604, 1.7637)), False)])  # VAL7 S5's maximum as a TS SP
def test_the_shortcut_judges_its_three_points_with_s2(tmp_path, monkeypatch, profile, seeded):
    """The shortcut's [minimum, TS SP, minimum] carries the TS SP's <S2>: S5's and S6's
    low-level TSs seed the saddle search; a TS SP off the minima's branch seeds nothing."""
    (e, s2) = profile
    ctx, _, notes = radical(tmp_path, **ends(profile), script=[(e[1] * K, s2[1])])
    monkeypatch.setattr(paths, "_xtb_ts_mode", lambda ctx, x: (1.0,) * 15)
    ts = ctx.geometry("ts", 0.5 * (ctx.ends[0] + ctx.ends[1]))
    ctx.case = ctx.case.model_copy(update={"low_level_ts": (ts,)})
    out = paths._shortcut(ctx)
    assert len(ctx.rt.qm.sps) == 1
    if seeded:
        assert out is not None and out[0].verdict == "single"
        assert [seed.source for seed in out[1]] == ["discovery_ts"]
    else:
        assert out is None and "shortcut:scf_branch_jump" in notes

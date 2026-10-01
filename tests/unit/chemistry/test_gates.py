from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from hfauto.chemistry import gates as g
from hfauto.chemistry.topology import bonds
from hfauto.chemistry.xyz import read_xyz
from hfauto.core import records as r
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Evidence, FileRef, Geometry, Level, Task
from hfauto.core.records import CaseOutcome

ACAC = Path(__file__).resolve().parents[3] / "configs" / "systems" / "xyz" / "acac"
REF = FileRef(path="jobs/x/out.txt", sha256="0" * 64)
LEVEL = Level(program="nwchem", version="7.2.3", method="pbe0", basis="def2-svpd",
              dispersion="d3bj", charge=0, multiplicity=1, grid="fine", scf_tol=1e-7)
E_TS = -100.0
K = 1.0 / HARTREE_TO_KCAL_MOL  # one kcal/mol in hartree


def modes(*low: float) -> tuple[float, ...]:  # 3N - 6 = 6 modes of a tetra-atomic
    return (*low, *(500.0 + 100.0 * i for i in range(6 - len(low))))


def ev(task: Task = "freq", freqs: tuple[float, ...] = modes(), *, fp: str = "g0",
       energy: float = E_TS, level: Level = LEVEL, s2: float | None = None) -> Evidence:
    geo = Geometry(file=REF, fingerprint=fp, symbols=("N", "H", "H", "H"))
    vibrations = {} if task != "freq" else {
        "frequencies_cm1": freqs, "n_external": 6, "hessian": REF,
        "imaginary_modes": tuple((0.0,) * 12 for f in freqs if f < 0)}
    return Evidence(engine="fake", task=task, level=level, start=geo, final=geo,
                    energy_hartree=energy, s2=s2, output=REF, job_key="k", **vibrations)


@pytest.mark.parametrize(("low", "notes"), [
    ((-27.09, 50.0), ("soft_imaginary_mode",)),
    ((6.83,), ()),
    ((-4.0,), ("noise_imaginary_mode",)),
])
def test_is_minimum_tiers(low, notes):
    gate = g.is_minimum(ev(freqs=modes(*low)), opt=ev("opt"))
    assert gate.ok and gate.notes == notes


def test_is_minimum_rejections():
    opt = ev("opt")
    assert g.is_minimum(ev(freqs=modes(-120.0)), opt=opt).reasons == ("imaginary_mode",)
    assert "geometry_mismatch" in g.is_minimum(ev(freqs=modes(), fp="g1"), opt=opt).reasons
    assert not g.is_minimum(ev(freqs=modes()), opt=ev("saddle"))
    fine_opt = ev("opt", level=LEVEL.model_copy(update={"grid": "xfine"}))
    assert "pes_mismatch:grid" in g.is_minimum(ev(freqs=modes()), opt=fine_opt).reasons


def test_freq_and_parent_share_one_electronic_state():
    """U3-P7: at one geometry the freq job reproduces its parent's energy within the SCF noise
    floor max(1e-5, 20 x scf_tol); another SCF solution 1e-4 Eh away is state_mismatch."""
    freqs = modes()
    assert g.is_minimum(ev(freqs=freqs), opt=ev("opt", energy=E_TS - 5e-6))
    shifted = ev("opt", energy=E_TS - 1e-4)
    assert g.is_minimum(ev(freqs=freqs), opt=shifted).reasons == ("state_mismatch",)
    moved = g.is_minimum(ev(freqs=freqs, fp="g1"), opt=shifted)  # energies of two geometries
    assert moved.reasons == ("geometry_mismatch",)
    saddle = ev("saddle", energy=E_TS + 1e-4)
    assert g.is_first_order_saddle(ev(freqs=modes(-700.0)), saddle=saddle).reasons == (
        "state_mismatch",)


@pytest.mark.parametrize("imag", [-1131.6, -680.11, -757.0])
def test_first_order_saddle_passes(imag):
    gate = g.is_first_order_saddle(ev(freqs=modes(imag)), saddle=ev("saddle"))
    assert gate.ok and gate.notes == ()


def test_saddle_mode_rules():  # one negative eigenvalue of any size above the noise
    saddle = ev("saddle")
    gate = g.is_first_order_saddle(ev(freqs=modes(-120.0, -60.0)), saddle=saddle)
    assert gate.reasons == ("higher_order",)
    gate = g.is_first_order_saddle(ev(freqs=modes(-700.0, -30.0)), saddle=saddle)
    assert gate.ok and gate.notes == ("soft_secondary_mode",)
    soft = g.is_first_order_saddle(ev(freqs=modes(-45.0)), saddle=saddle)  # torsion or not
    assert soft.ok and soft.notes == ()
    noise = g.is_first_order_saddle(ev(freqs=modes(-8.0)), saddle=saddle)
    assert noise.reasons == ("no_imaginary_mode",)


AHB = np.array([[-1.4, 0.0, 0.0], [0.04, 0.0, 0.0], [1.4, 0.0, 0.0]])  # A-H broken, H-B formed
AHB_BONDS = frozenset({(0, 1), (1, 2)})


def test_reaction_mode_character_on_the_changed_bonds():
    """G1-P2: chi = |Q_B^T q| on the orthonormalized stretches of the changed bonds."""
    antisymmetric = np.array([1, 0, 0, -2, 0, 0, 1, 0, 0])  # A-H shortens as H-B stretches
    h_alone = np.eye(9)[3]  # H moves along the axis, A and B stay: sqrt(2/3)
    rotation = np.cross([0.3, -0.2, 1.0], AHB + [0, 0.5, 0]).ravel()  # rigid, about any point
    bent = AHB + [[0, 0.5, 0], [0, 0, 0], [0, 0, 0]]
    rotor = (np.cross([0, 0, 1.0], bent - bent[1]) * [[1], [0], [0]]).ravel()  # A turns about H
    chi = [g.reaction_mode_chi(-m, x, AHB_BONDS) for m, x in
           ((antisymmetric, AHB), (h_alone, AHB), (rotation, AHB), (rotor, bent))]
    assert chi[:2] == pytest.approx([1.0, np.sqrt(2 / 3)])
    assert chi[2:] == pytest.approx([0.0, 0.0], abs=1e-12)
    ring = AHB_BONDS | {(0, 2)}  # three collinear stretches span two directions, not three
    assert g.reaction_mode_chi(np.tile([1.0, 0, 0], 3), AHB, ring) == pytest.approx(0.0)
    assert g.reaction_mode_character(h_alone, AHB, AHB_BONDS)
    assert g.reaction_mode_character(rotor, bent, AHB_BONDS).reasons == ("not_reaction_mode",)


def test_reaction_mode_character_without_a_bond_change():
    """A declared coordinate only: |cos(q, grad q)|; with neither the gate is not applied."""
    gradient = np.eye(9)[4]
    mode = np.eye(9)[4] + np.eye(9)[1] * np.tan(np.radians(80))  # 80 degrees off the gradient
    assert g.reaction_mode_chi(-mode, AHB, frozenset(), gradient) == pytest.approx(
        np.cos(np.radians(80)))
    assert g.reaction_mode_character(mode, AHB, frozenset(), gradient).reasons == (
        "not_reaction_mode",)
    assert g.reaction_mode_chi(mode, AHB, frozenset()) is None
    assert g.reaction_mode_character(mode, AHB, frozenset(), np.zeros(9))


def test_same_pes_and_spin():
    xfine = LEVEL.model_copy(update={"grid": "xfine"})
    assert g.same_pes(LEVEL, LEVEL)
    assert g.same_pes(LEVEL, xfine).reasons == ("pes_mismatch:grid",)
    assert g.same_pes(LEVEL, xfine, numerics=False)
    assert not g.same_pes(LEVEL, xfine.model_copy(update={"charge": 1}), numerics=False)
    doublet = LEVEL.model_copy(update={"multiplicity": 2})
    anion = doublet.model_copy(update={"charge": -1})  # complex vs monomer: state not compared
    assert g.same_pes(LEVEL, anion, state=False) and not g.same_pes(LEVEL, anion)
    assert g.same_pes(LEVEL, xfine, state=False).reasons == ("pes_mismatch:grid",)
    assert g.spin_ok(ev(level=doublet, s2=0.7523))
    assert g.spin_ok(ev(level=doublet, s2=0.90)).reasons == ("spin_contaminated",)
    assert g.spin_ok(ev(level=doublet, s2=0.90), g.Policy(spin_tol=0.2))
    assert g.spin_ok(ev())


DOWN = (E_TS - 1e-3, E_TS - 5e-3, E_TS - 1e-2)
AB = frozenset({"A", "B"})


def connect(sides=(DOWN, DOWN), assigned=("A", "B"), expected=AB, **kw):
    pair = (ev("opt", energy=sides[0][-1]), ev("opt", energy=sides[1][-1]))  # where each ends
    return g.connection(ev(freqs=modes(-700.0)), pair, assigned, expected,
                        degenerate=kw.pop("degenerate", False), **kw)


def test_connection_judges_where_each_side_ends():
    assert connect() == (g.Gate(True), "elementary")
    above = (E_TS + 1e-4, E_TS + 2e-4, E_TS - 1e-2)  # displaced above the TS, then down (H2O2)
    assert connect(sides=(above, above)) == (g.Gate(True), "elementary")
    gate, label = connect(sides=((E_TS - 1e-3, E_TS - 5e-6), DOWN))  # ends within drop of E_TS
    assert label == "failed" and gate.reasons == ("side0:no_descent",)


def test_a_ts_below_the_declared_product_survives_as_reassigned():
    """Endothermic A -> B with E_TS < E_B: the saddle of A -> I is a TS all the same (no endpoint
    energy enters the TS gate) and QRC names the minima it connects."""
    assert g.is_first_order_saddle(ev(freqs=modes(-700.0)), saddle=ev("saddle"))
    assert connect(assigned=("A", "I")) == (g.Gate(True), "reassigned")


def test_qrc_drop_is_the_scf_noise_floor():
    assert g.qrc_drop(LEVEL) == pytest.approx(1e-5)  # 20 x 1e-7 is below the minimum
    loose = LEVEL.model_copy(update={"scf_tol": 1e-5})
    assert g.qrc_drop(loose) == pytest.approx(2e-4)
    side = ev("opt", energy=E_TS - 1e-4, level=loose)  # 1e-4 below the TS: less than 2e-4
    gate, _ = g.connection(ev(freqs=modes(-700.0), level=loose), (side, side), ("A", "B"), AB,
                           degenerate=False)
    assert gate.reasons == ("side0:no_descent", "side1:no_descent")


def test_connection_assignment():
    assert connect(assigned=("B", "A"))[1] == "elementary"
    one = frozenset({"M"})
    assert connect(assigned=("M", "M"), expected=one, degenerate=True)[1] == "degenerate"
    undistinct = connect(assigned=("M", "M"), expected=one, degenerate=True, sides_distinct=False)
    assert undistinct[1] == "failed"
    assert connect(assigned=("A", "C")) == (g.Gate(True), "reassigned")
    assert connect(assigned=("A", None))[0].reasons == ("unassigned_side",)
    assert connect(assigned=("A", "A")) == (g.Gate(False, ("sides_same_basin",)), "failed")


def test_a_degenerate_proton_transfer_moves_the_proton():
    """U6-P7 (S16 acac enol): the PT's QRC sides carry the two endpoints' bond graphs. A methyl
    rotation TS connects the basin to itself through two rotamers of R, with the proton in
    place: no degenerate PT."""
    reactant, product = (read_xyz(ACAC / f"{name}.xyz") for name in ("reactant", "product"))
    r_bonds, p_bonds = (bonds(reactant.symbols, end.coords) for end in (reactant, product))
    assert r_bonds != p_bonds  # H10 on O2, then on O6
    kw = {"assigned": ("M", "M"), "expected": frozenset({"M"}), "degenerate": True}
    pt, rotation = (r_bonds, p_bonds, p_bonds, r_bonds), (r_bonds, p_bonds, r_bonds, r_bonds)
    assert connect(bond_sets=pt, **kw) == (g.Gate(True), "degenerate")
    assert connect(bond_sets=rotation, **kw) == (g.Gate(False, ("bond_change_missing",)), "failed")
    assert connect(bond_sets=(r_bonds,) * 4, **kw)[1] == "degenerate"  # inversion: no bond change
    assert connect(bond_sets=rotation, assigned=("A", "B"))[1] == "elementary"  # not degenerate


def reaction(**kw) -> r.ReactionRecord:
    term = (r.StoichTerm(composition_id="c", coefficient=1),)
    base = r.ReactionRecord(reaction_id="r", reactants=term, products=term, minima=("A", "B"),
                            endpoints=("sa", "sb"), source="declared",
                            outcome=CaseOutcome.ELEMENTARY_STEP)
    return base.model_copy(update=kw)


def thermo(dG_eff: float | None = 12.2, blockers: tuple[str, ...] = ()) -> r.ReactionThermo:
    return r.ReactionThermo(reaction_id="r", T_K=298.15, standard_state="1atm", dE_act_kcal=13.6,
                            dE_rxn_kcal=1.0, dG_act_kcal=12.2,
                            dG_rxn_kcal=1.0, dG_eff_kcal=dG_eff, blockers=blockers)


def test_discovery_verdict_window_edges():
    """U4-P6: the DFT reaction window (40 kcal/mol) and the 50 kcal/mol low-level barrier cap,
    both inclusive; an NT2 product needs its validated TS."""

    def verdict(act=10.0, rxn=5.0, ts=True, **policy):
        return g.discovery_verdict(ts_validated=ts, dE_act_kcal=act,
                                   dE_rxn_kcal=rxn, policy=g.Policy(**policy))

    assert verdict() is None and verdict(act=50.0, rxn=40.0) is None
    assert verdict(act=50.01) == verdict(rxn=40.01) == verdict(rxn=None) == "out_of_window"
    assert verdict(ts=False) == "ts_not_validated"
    assert verdict(rxn=45.0, reaction_window_kcal=50.0) is None


def test_rankable():  # any dZPE; a barrierless outcome ranks by its dG_eff
    assert g.rankable(reaction(), thermo())
    assert g.rankable(reaction(outcome=CaseOutcome.BARRIERLESS), thermo(1.0))
    assert not g.rankable(reaction(), thermo(blockers=("spin_contaminated",)))
    assert g.rankable(reaction(), None).reasons == ("thermo_record_missing",)
    lost = thermo(None, ("mixed_level_of_theory", "thermo_unavailable"))
    assert g.rankable(reaction(), lost).reasons == lost.blockers  # the thermo stage's facts only
    same = reaction(outcome=CaseOutcome.SAME_BASIN)  # the outcome and the thermo blockers
    assert g.rankable(same, None).reasons == ("outcome:same_basin", "thermo_record_missing")
    assert g.rankable(same, thermo(0.0, ("spin_contaminated",))).reasons == (
        "outcome:same_basin", "spin_contaminated")


def test_reaction_tier():
    saddle = r.SaddleClaim(saddle_calc="s", freq_calc="f", imag_cm1=-700.0, energy_hartree=E_TS)
    link = r.ConnectionClaim(side_calcs=("p", "m"), minima=("A", "B"))
    assert g.reaction_tier(reaction(outcome=None)) == "screening"
    assert g.reaction_tier(reaction(outcome=CaseOutcome.BLOCKED)) == "screening"
    assert g.reaction_tier(reaction(outcome=CaseOutcome.SAME_BASIN)) == "minima"
    assert g.reaction_tier(reaction(saddle=saddle)) == "saddle"
    assert g.reaction_tier(reaction(saddle=saddle, connection=link)) == "connected"


def test_barrier_verdict_classifies_a_profile_between_the_dft_minima():
    single = g.barrier_verdict((0.0, 2 * K, 3 * K, 1 * K, 0.5 * K), source="screen")
    assert (single.verdict, single.source) == ("single", "screen")
    well = g.barrier_verdict((0.0, 2 * K, 0.5 * K, 3 * K, 6 * K, 9 * K, 10 * K), source="string")
    assert well.verdict == "intermediate"
    flat = g.barrier_verdict((0.0, 0.2 * K, 0.4 * K, 0.5 * K), source="string")
    assert (flat.verdict, flat.reasons) == ("barrierless", ())
    none = g.barrier_verdict((0.0, 5 * K), source="screen")
    assert (none.verdict, none.reasons) == ("unavailable", ("too_few_points",))

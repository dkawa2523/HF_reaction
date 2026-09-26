from __future__ import annotations

import pytest

from hfauto.chemistry import gates as g
from hfauto.core import records as r
from hfauto.core.constants import CM1_TO_HARTREE, HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Evidence, FileRef, Geometry, Level, Task
from hfauto.core.records import CaseOutcome

REF = FileRef(path="jobs/x/out.txt", sha256="0" * 64)
LEVEL = Level(program="nwchem", version="7.2.3", method="pbe0", basis="def2-svpd",
              dispersion="d3bj", charge=0, multiplicity=1, grid="fine", scf_tol=1e-7)
E_TS = -100.0
K = 1.0 / HARTREE_TO_KCAL_MOL  # one kcal/mol in hartree


def modes(*low: float) -> tuple[float, ...]:  # 3N - 6 = 6 modes of a tetra-atomic
    return (*low, *(500.0 + 100.0 * i for i in range(6 - len(low))))


def ev(task: Task = "freq", freqs: tuple[float, ...] | None = None, *, fp: str = "g0",
       energy: float = E_TS, traj: tuple[float, ...] = (), level: Level = LEVEL,
       s2: float | None = None) -> Evidence:
    geo = Geometry(file=REF, fingerprint=fp, symbols=("N", "H", "H", "H"))
    return Evidence(engine="fake", task=task, level=level, start=geo, final=geo,
                    energy_hartree=energy, trajectory_energies_hartree=traj,
                    frequencies_cm1=freqs, n_external=6 if freqs is not None else None,
                    s2=s2, output=REF, job_key="k")


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
    assert not g.is_minimum(ev(freqs=modes()[:5]), opt=opt)
    assert not g.is_minimum(ev(freqs=modes()), opt=ev("saddle"))
    fine_opt = ev("opt", level=LEVEL.model_copy(update={"grid": "xfine"}))
    assert "pes_mismatch:grid" in g.is_minimum(ev(freqs=modes()), opt=fine_opt).reasons


@pytest.mark.parametrize("imag", [-1131.6, -680.11, -757.0])
def test_first_order_saddle_passes(imag):
    gate = g.is_first_order_saddle(ev(freqs=modes(imag)), saddle=ev("saddle"))
    assert gate.ok and gate.notes == ()


def test_saddle_mode_rules():
    saddle = ev("saddle")
    gate = g.is_first_order_saddle(ev(freqs=modes(-120.0, -60.0)), saddle=saddle)
    assert not gate and "higher_order" in gate.reasons
    gate = g.is_first_order_saddle(ev(freqs=modes(-700.0, -30.0)), saddle=saddle)
    assert gate.ok and gate.notes == ("soft_secondary_mode",)
    assert not g.is_first_order_saddle(ev(freqs=modes(-45.0)), saddle=saddle)
    assert g.is_first_order_saddle(ev(freqs=modes(-45.0)), saddle=saddle, torsional=True)


def test_saddle_prominence():
    freq, saddle = ev(freqs=modes(-700.0)), ev("saddle")
    low = g.is_first_order_saddle(freq, saddle=saddle, endpoint_energies=(E_TS - 1e-8,))
    assert low.reasons == ("low_prominence",)
    ends = (E_TS - 1e-3, E_TS - 2e-3)
    assert g.is_first_order_saddle(freq, saddle=saddle, endpoint_energies=ends)
    loose = LEVEL.model_copy(update={"scf_tol": 1e-5})  # floor becomes 20 x 1e-5 = 2e-4
    loose_freq = ev(freqs=modes(-700.0), level=loose)
    assert not g.is_first_order_saddle(loose_freq, saddle=ev("saddle", level=loose),
                                       endpoint_energies=(E_TS - 1e-4,))


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
    assert g.spin_ok(ev())


DOWN = (E_TS - 1e-3, E_TS - 5e-3, E_TS - 1e-2)
AB = frozenset({"A", "B"})


def connect(sides=(DOWN, DOWN), assigned=("A", "B"), expected=AB, **kw):
    pair = (ev("opt", traj=sides[0], energy=sides[0][-1]),
            ev("opt", traj=sides[1], energy=sides[1][-1]))
    return g.connection(ev(freqs=modes(-700.0)), pair, assigned, expected,
                        degenerate=kw.pop("degenerate", False), **kw)


def test_connection_trajectory_checks():
    assert connect() == (g.Gate(True), "elementary")
    gate, label = connect(sides=(DOWN, (E_TS - 1e-3, E_TS + 1e-4, E_TS - 1e-2)))
    assert label == "failed" and gate.reasons == ("side1:trajectory_above_ts",)
    gate, label = connect(sides=((E_TS - 1e-3, E_TS - 1.005e-3), DOWN))
    assert label == "failed" and gate.reasons == ("side0:no_descent",)
    gate, _ = connect(sides=((E_TS + 1e-4, E_TS - 1e-2), DOWN))
    assert "side0:no_initial_descent" in gate.reasons


def test_connection_assignment():
    assert connect(assigned=("B", "A"))[1] == "elementary"
    one = frozenset({"M"})
    assert connect(assigned=("M", "M"), expected=one, degenerate=True)[1] == "degenerate"
    undistinct = connect(assigned=("M", "M"), expected=one, degenerate=True, sides_distinct=False)
    assert undistinct[1] == "failed"
    assert connect(assigned=("A", "C")) == (g.Gate(True), "reassigned")
    assert connect(assigned=("A", None))[0].reasons == ("unassigned_side",)
    assert connect(assigned=("A", "A"))[1] == "failed"


def test_thermo_consistent():  # against the list GoodVibes got (chemistry.thermo_frequencies)
    real = (600.0, 800.0, 1300.0, 1700.0)
    sent = (*real, 2 * 0.018747 / CM1_TO_HARTREE - sum(real))  # TS without -680.11: ZPE 0.018747
    freq = ev(freqs=(-680.11, *sent))
    kw = {"frequencies_cm1": sent, "gv_energy_hartree": E_TS, "gv_n_real": 5, "scale": 0.985}
    assert g.thermo_consistent(freq, gv_zpe_hartree=0.0375115, **kw).reasons == ("zpe_mismatch",)
    assert g.thermo_consistent(freq, gv_zpe_hartree=0.985 * 0.018747, **kw)
    assert g.thermo_consistent(freq, gv_zpe_hartree=0.985 * 0.018747, **kw | {"gv_n_real": 6}
                               ).reasons == ("n_real:6!=5",)
    soft = modes(-30.0)
    flipped = (30.0, *soft[1:])  # a soft minimum's mode enters as |nu|
    gv_zpe = 0.5 * sum(flipped) * CM1_TO_HARTREE
    assert g.zpe_hartree(flipped) == pytest.approx(gv_zpe, abs=1e-15)
    assert g.zpe_hartree(soft) < gv_zpe  # negative modes never count
    kw = {"gv_zpe_hartree": gv_zpe, "gv_energy_hartree": E_TS, "gv_n_real": 6, "scale": 1.0}
    assert g.thermo_consistent(ev(freqs=soft), frequencies_cm1=flipped, **kw)
    assert not g.thermo_consistent(ev(freqs=soft), frequencies_cm1=soft, **kw)
    assert not g.thermo_consistent(ev(freqs=soft), frequencies_cm1=flipped,
                                   **kw | {"gv_energy_hartree": E_TS + 1e-5})


def reaction(**kw) -> r.ReactionRecord:
    term = (r.StoichTerm(composition_id="c", coefficient=1),)
    base = r.ReactionRecord(reaction_id="r", reactants=term, products=term, minima=("A", "B"),
                            endpoints=("sa", "sb"), source="declared", n_h_transferred=2,
                            outcome=CaseOutcome.ELEMENTARY_STEP)
    return base.model_copy(update=kw)


def thermo(dzpe: float, blockers: tuple[str, ...] = ()) -> r.ReactionThermo:
    return r.ReactionThermo(reaction_id="r", T_K=298.15, standard_state="1atm", dE_act_kcal=13.6,
                            dE_rxn_kcal=1.0, dzpe_act_kcal=dzpe, dG_act_kcal=12.2,
                            dG_rxn_kcal=1.0, blockers=blockers)


def test_rankable():
    assert g.rankable(reaction(), thermo(-9.0))
    assert g.rankable(reaction(), thermo(9.1)).reasons == ("dzpe_out_of_tolerance",)
    assert not g.rankable(reaction(n_h_transferred=0), thermo(5.5))
    assert not g.rankable(reaction(), thermo(1.0, ("spin_contaminated",)))
    assert not g.rankable(reaction(), thermo(1.0), participant_notes=("mixed_level_of_theory",))
    assert "thermo_unavailable" in g.rankable(reaction(), None).reasons
    assert not g.rankable(reaction(outcome=CaseOutcome.BARRIERLESS), thermo(1.0))


def test_reaction_tier():
    saddle = r.SaddleClaim(saddle_calc="s", freq_calc="f", imag_cm1=-700.0, energy_hartree=E_TS)
    link = r.ConnectionClaim(side_calcs=("p", "m"), minima=("A", "B"), amplitude_A=0.1)
    assert g.reaction_tier(reaction(outcome=None)) == "screening"
    assert g.reaction_tier(reaction(outcome=CaseOutcome.BLOCKED)) == "screening"
    assert g.reaction_tier(reaction(outcome=CaseOutcome.SAME_BASIN)) == "minima"
    assert g.reaction_tier(reaction(saddle=saddle)) == "saddle"
    assert g.reaction_tier(reaction(saddle=saddle, connection=link)) == "connected"


ENDS = (0.0, 0.5 * K)
MONOTONIC = (0.0, 0.2 * K, 0.4 * K, 0.5 * K)


def test_barrier_proceed():
    geo = ev().final
    v = g.barrier_verdict((0.0, 2 * K, 3 * K, 1 * K, 0.5 * K), dft_endpoints=ENDS,
                          low_profile=(0.0, K, 2 * K, 3 * K), seed=geo)
    assert v.verdict == "proceed" and v.seed == geo
    assert v.max_rel_dft_kcal == pytest.approx(2.5) and v.max_rel_low_kcal == pytest.approx(-1.0)
    low_max = g.barrier_verdict(MONOTONIC, dft_endpoints=ENDS, low_profile=(0.0, 2 * K, 0.0))
    assert low_max.verdict == "proceed"


def test_barrierless_and_unavailable():
    v = g.barrier_verdict(MONOTONIC, dft_endpoints=ENDS, low_profile=MONOTONIC,
                          tangent_mode_cm1=1000.0, max_node_spacing_A=0.12)
    assert (v.verdict, v.below_zpe, v.max_node_spacing_A, v.n_dft_points, v.seed) == (
        "barrierless", True, 0.12, 4, None)
    assert g.barrier_verdict((0.0, 5 * K), dft_endpoints=ENDS).verdict == "unavailable"

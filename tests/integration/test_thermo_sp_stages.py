"""sp -> thermo with FakeQM and fake_species_thermo: energy layer (parents, label, spin state),
band, the chain's separated zero (declared monomers NH + O(3P)), association, mixed LOT, m = 2,
the electronic term (the O atom's NIST levels, declared levels; E_SO in dE and G), the state G
(lowest spin-clean minimum, closed when one lacks its energy layer) and dG_eff (submerged
barrier; none for a barrierless step)."""

import dataclasses
import math
from pathlib import Path

import numpy as np
import pytest
from fakes import PES, FakeQM, double_well, fake_species_thermo, write_geometry, write_hessian

from hfauto.backends.protocols import Capability
from hfauto.chemistry import symmetry, thermo
from hfauto.chemistry.electronic_state import atom_levels
from hfauto.chemistry.topology import state_label
from hfauto.chemistry.vibrations import projected_frequencies
from hfauto.core import records as R
from hfauto.core.constants import HARTREE_TO_KCAL_MOL, R_KCAL_MOL_K
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec
from hfauto.core.system import CompositionInput, SpeciesInput, SystemConfig
from hfauto.drivers.minimum import calc_id
from hfauto.stages.single_point import SinglePointConfig, SinglePointStage
from hfauto.stages.thermochemistry import ThermoConfig, ThermoStage

T = R.ArtifactType
DFT, BIG = (MethodSpec(id=b, kind="dft", functional="xfake", basis=b) for b in ("svp", "tzvp"))
NEUTRAL = {"complex": (0, 1), "nh": (0, 1), "o": (0, 3)}  # (charge, multiplicity)
# the monomers' states: the reactant's fragments (N-H 1.0 A, O 2.8 A away) are NH and O
LABEL = {"nh": state_label(["N", "H"], np.eye(3)[:2]), "o": state_label(["O"], np.eye(3)[:1])}


@pytest.fixture(autouse=True)
def _closed_form_thermo(monkeypatch):
    monkeypatch.setattr(thermo, "species_thermo", fake_species_thermo)


def _state(ev, charge, multiplicity):
    level = ev.level.model_copy(update={"charge": charge, "multiplicity": multiplicity})
    return ev.model_copy(update={"level": level})


def _setup(fake_runtime, tmp_run, states=NEUTRAL, nh_levels=None):
    """The double well as composition c = NH + O (p2: the product at another grid; o2: a second
    minimum of the O state), rx1 = reactant -> product and rx2 = reactant -> p2; ``nh_levels``:
    NH's declared electronic_levels."""
    well = double_well()
    part = well.energy(well.points["reactant"]) / 2 + 0.05  # NH and O: 63 kcal/mol up
    pes = PES(well.symbols, lambda x: well.energy(x) if np.size(x) == 9 else part, well.points)
    qm = FakeQM(tmp_run, pes)
    ev = {p: _state(qm.frequencies(pes.molecule(p), DFT), *states["complex"])
          for p in ("reactant", "product", "ts")}
    medium = ev["product"].level.model_copy(update={"grid": "medium"})
    ev["p2"] = ev["product"].model_copy(update={"job_key": "p2", "level": medium})  # other LOT
    stretch = 0.25 * np.outer([-1, 1, 0], [-1, 1, 0])  # N-H along (H - N), 0.5 Eh/bohr²
    hessians = {"nh": np.kron([[1, -1], [-1, 1]], stretch), "o": np.zeros((3, 3))}
    for tag, symbols, external in (("nh", ["N", "H"], 5), ("o", ["O"], 3)):
        x = np.eye(3)[: len(symbols)]
        geom = write_geometry(tmp_run, f"{tag}.xyz", symbols, x)
        nu = tuple(map(float, projected_frequencies(hessians[tag], symbols, x)[0]))
        ev[tag] = _state(ev["reactant"], *states[tag]).model_copy(update={
            "start": geom, "final": geom, "job_key": tag, "frequencies_cm1": nu,
            "n_external": external, "energy_hartree": part,
            "hessian": write_hessian(tmp_run, f"{tag}_hessian.npy", hessians[tag])})
    ev["o2"] = ev["o"].model_copy(update={"job_key": "o2"})
    arts = [Artifact(artifact_id=calc_id(e), type=T.CALCULATION, payload=e) for e in ev.values()]
    arts += [Artifact(artifact_id=t, type=T.SPECIES, payload=R.SpeciesRecord(
        species_id=t, composition_id=t, charge=states[t][0], multiplicity=states[t][1],
        geometry=ev[t].final, source="input", state_label=LABEL[t])) for t in ("nh", "o")]
    arts += [Artifact(artifact_id=f"m_{k}", type=T.MINIMUM, payload=R.MinimumRecord(
        minimum_id=f"m_{k}", basin_id=k, composition_id={"nh": "nh", "o": "o", "o2": "o"}.get(
            k, "HNO"), species_id=k, tier="dft", level_key="x", opt_calc="o", freq_calc=calc_id(ev[k]),
        energy_hartree=0.0, state_label=LABEL.get({"o2": "o"}.get(k, k), k)))
        for k in ev if k != "ts"]
    term = (R.StoichTerm(composition_id="HNO", coefficient=1),)
    rx1 = R.ReactionRecord(
        reaction_id="rx1", reactants=term, products=term, minima=("m_reactant", "m_product"),
        endpoints=("a", "b"), source="declared", outcome=R.CaseOutcome.ELEMENTARY_STEP,
        saddle=R.SaddleClaim(saddle_calc="s", freq_calc=calc_id(ev["ts"]), imag_cm1=-900.0,
                             energy_hartree=0.0))
    rx2 = rx1.model_copy(update={"reaction_id": "rx2", "minima": ("m_reactant", "m_p2")})
    arts += [Artifact(artifact_id=r.reaction_id, type=T.REACTION, payload=r) for r in (rx1, rx2)]
    system = SystemConfig(system_id="s", species=[  # only the levels are read
        SpeciesInput(id="nh", xyz=Path("nh.xyz"), multiplicity=1, electronic_levels=nh_levels),
        SpeciesInput(id="o", xyz=Path("o.xyz"), multiplicity=3)],
        compositions=[CompositionInput(id="c", components={"nh": 1, "o": 1})])
    rt = fake_runtime(system, {(Capability.QM, "nwchem"): qm}, methods={"svp": DFT, "tzvp": BIG})
    return Manifest(run_id="r", stage_id="v", created_at="t", artifacts=arts), rt, ev


def _thermo(view, rt, method=None):  # <subject|rxn>_<T>... -> payload
    config = ThermoConfig(energy_method=method, standard_states=("1atm", "1M"))
    out = ThermoStage().run(view, config, rt)
    return {a.artifact_id.split("_thermo_")[1]: a.payload for a in out}


def _sp(inputs, rt):
    return SinglePointStage().run(inputs, SinglePointConfig(engine="nwchem", methods=["tzvp"]), rt)


def test_sp_then_thermo(fake_runtime, tmp_run):
    inputs, rt, ev = _setup(fake_runtime, tmp_run)
    sp = _sp(inputs, rt)
    p2 = next(a for a in sp if "m_p2" in a.parents)  # one geometry, but its own freq's SCF
    assert p2.parents == ("m_p2",)
    assert p2.payload.start.fingerprint == ev["product"].final.fingerprint
    view = inputs.model_copy(update={"artifacts": [*inputs.artifacts, *sp]})

    plain = _thermo(view, rt)
    rx = plain["rx1_298.15K_1atm"]
    assert rx.blockers == () and rx.band_kcal[0] < rx.dG_act_kcal < rx.band_kcal[1]
    assert rx.dG_assoc_kcal < 0 and plain["rx1_298.15K_1M"].dG_assoc_kcal == pytest.approx(
        rx.dG_assoc_kcal - 1.894, abs=1e-3)  # one molecule against two
    assert rx.dG_act_vs_separated_kcal == pytest.approx(rx.dG_assoc_kcal + rx.dG_act_kcal)
    assert rx.dG_eff_kcal == pytest.approx(rx.dG_act_kcal) and rx.notes == ()  # one conformer
    assert rx.reference == "complex"  # bound by 63 kcal/mol: the complex is the zero
    assert rx.band_kcal[0] <= rx.dG_eff_kcal <= rx.band_kcal[1]
    assert plain["rx2_298.15K_1atm"].blockers == ("mixed_level_of_theory",)
    assert rx.energy_level == plain["rx2_298.15K_1atm"].energy_level == "xfake/svp"  # no grid
    composite = _thermo(view, rt, "tzvp")
    assert composite["m_reactant_298.15K"].energy_calc  # the sp its parents name
    assert composite["rx1_298.15K_1atm"].energy_level == "xfake/tzvp"
    assert composite["rx1_298.15K_1atm"].blockers == ()
    uncertified = [a.model_copy(update={"payload": a.payload.model_copy(update={  # X1-2
        "saddle": a.payload.saddle.model_copy(update={"notes": ("not_stationary",)})})})
        if a.type == T.REACTION else a for a in view.artifacts]  # rx1 and rx2 share the TS
    blocked = _thermo(view.model_copy(update={"artifacts": uncertified}), rt, "tzvp")
    assert blocked["rx1_298.15K_1atm"].blockers == ("not_stationary",)

    ts = calc_id(ev["ts"])  # X3: the layer counts on its freq's spin state (S5's -83.15 row)

    def spin(freq_s2, sp_s2):  # the TS freq's and its layer's <S2> (a singlet: S(S+1) = 0)
        s2 = [a.model_copy(update={"payload": a.payload.model_copy(update={
            "s2": freq_s2 if a.artifact_id == ts else sp_s2})})
            if ts in (a.artifact_id, *a.parents) else a for a in view.artifacts]
        return _thermo(view.model_copy(update={"artifacts": s2}), rt, "tzvp")["rx1_298.15K_1atm"]

    assert spin(0.0, 0.05).blockers == ()
    for left_its_state in (spin(0.0, 0.75), spin(1.0, 1.02)):  # another state; a BS freq
        assert left_its_state.blockers == ("spin_contaminated",)
        assert left_its_state.dG_eff_kcal is not None

    # S9: the sp layer was asked for, the TS has none -> fail closed
    view = inputs.model_copy(update={"artifacts": [
        *inputs.artifacts, *(a for a in sp if ts not in a.parents)]})
    layered = _thermo(view, rt, "tzvp")
    assert layered[f"{ts}_298.15K"].G_hartree is None
    assert layered[f"{ts}_298.15K"].notes == ("energy_layer_missing",)
    assert "thermo_unavailable" in layered["rx1_298.15K_1atm"].blockers
    assert layered["rx1_298.15K_1atm"].dE_act_kcal is None  # no dE mixing the sp and freq LOTs
    assert layered["rx1_298.15K_1atm"].energy_level is None


def test_association_across_charge_and_spin_fails_closed(fake_runtime, tmp_run):  # S14
    anion = {"complex": (-1, 3), "nh": (-1, 1), "o": (0, 3)}  # F-.HF-like, with a triplet
    inputs, rt, _ = _setup(fake_runtime, tmp_run, anion)
    plain = _thermo(inputs, rt)
    rx, G = plain["rx1_298.15K_1atm"], {k: plain[f"m_{k}_298.15K"].G_hartree for k in (
        "reactant", "nh", "o", "o2")}
    assert rx.blockers == () and rx.dG_act_vs_separated_kcal is not None
    assert G["o"] == G["o2"]  # two minima of the O state count once, not as an ensemble
    assert rx.dG_assoc_kcal == pytest.approx(
        (G["reactant"] - G["nh"] - G["o"]) * HARTREE_TO_KCAL_MOL) and rx.dG_assoc_kcal < 0
    sp = [a.model_copy(update={"parents": tuple(p for p in a.parents if p not in ("m_o", "m_o2"))})
          for a in _sp(inputs, rt)]  # the O state loses its energy layer
    view = inputs.model_copy(update={"artifacts": [*inputs.artifacts, *sp]})
    out = _thermo(view, rt, "tzvp")
    assert out["m_nh_298.15K"].G_hartree is not None and out["m_o_298.15K"].G_hartree is None
    closed = out["rx1_298.15K_1atm"]  # the declared separated zero is a candidate: no fallback
    assert closed.dG_assoc_kcal is None and closed.dG_eff_kcal is None
    assert closed.blockers == ("thermo_unavailable",) and closed.dG_act_kcal is not None


def _with(view, *reactions):
    return view.model_copy(update={"artifacts": [*view.artifacts, *(
        Artifact(artifact_id=r.reaction_id, type=T.REACTION, payload=r) for r in reactions)]})


def _undeclared(view):  # no species records: no declared monomers, no separated points
    return view.model_copy(update={"artifacts": [a for a in view.artifacts
                                                 if a.type != T.SPECIES]})


def test_an_association_refers_to_its_separated_monomers(fake_runtime, tmp_run):
    """X4: NH + O -> HNO asked from its separated monomers (a triplet anion from a singlet anion
    and a triplet): dG_act and dG_rxn refer to them, with the shift of their count (dn = -1). A
    barrierless one has no dG_eff (capture-limited) and its precursor complex, an auxiliary
    point, gives dG_assoc. With a TS the complex is R, a candidate for the chain's zero, so a
    broken-symmetry complex (no clean G) closes dG_eff but blocks no barrierless row."""
    anion = {"complex": (-1, 3), "nh": (-1, 1), "o": (0, 3)}
    inputs, rt, ev = _setup(fake_runtime, tmp_run, anion)
    rx1 = inputs.get("rx1").payload
    separated = tuple(R.StoichTerm(composition_id=c, coefficient=1) for c in ("nh", "o"))
    bound = rx1.model_copy(update={"reaction_id": "bound", "reactants": separated,
                                   "monomers": ("m_nh", "m_o")})
    flat = bound.model_copy(update={"reaction_id": "flat", "saddle": None,
                                    "outcome": R.CaseOutcome.BARRIERLESS})
    view = _with(inputs, bound, flat)
    out = _thermo(view, rt)
    G = {k: out[f"m_{k}_298.15K"].G_hartree for k in ("reactant", "product", "nh", "o")}
    G_ts = out[f"{calc_id(ev['ts'])}_298.15K"].G_hartree
    dG_rxn = (G["product"] - G["nh"] - G["o"]) * HARTREE_TO_KCAL_MOL
    rx, molar = out["flat_298.15K_1atm"], out["flat_298.15K_1M"]
    assert rx.dG_rxn_kcal == pytest.approx(dG_rxn) and dG_rxn < 0
    assert molar.dG_rxn_kcal == pytest.approx(dG_rxn - 1.894, abs=1e-3)  # dn = -1
    assert (rx.dG_eff_kcal, rx.dG_act_kcal, rx.reference, rx.blockers) == (None, None, None, ())
    so = thermo.spin_orbit(atom_levels("O", 0, 3))  # the O(3P) atom's E_SO is in its E
    assert rx.dE_rxn_kcal == pytest.approx((ev["product"].energy_hartree - ev["nh"].energy_hartree
                                            - ev["o"].energy_hartree - so) * HARTREE_TO_KCAL_MOL)
    assert rx.dG_assoc_kcal == pytest.approx(
        (G["reactant"] - G["nh"] - G["o"]) * HARTREE_TO_KCAL_MOL)
    ts = out["bound_298.15K_1atm"]  # the TS lies below the monomers: submerged, dropped
    assert ts.dG_act_kcal == pytest.approx((G_ts - G["nh"] - G["o"]) * HARTREE_TO_KCAL_MOL)
    assert ts.dG_act_vs_separated_kcal == pytest.approx(ts.dG_act_kcal) and ts.dG_act_kcal < 0
    assert ts.notes == ("submerged_barrier",) and ts.reference == "complex"
    assert ts.dG_eff_kcal == pytest.approx(max(  # the chain [R_sep, R, P]
        (G["product"] - G["reactant"]) * HARTREE_TO_KCAL_MOL, 0.0))
    hot = _thermo(_relabel(view, {calc_id(ev["reactant"]): {"s2": 3.0}}), rt)  # a BS complex
    assert hot["bound_298.15K_1atm"].blockers == ("thermo_unavailable",)
    assert hot["flat_298.15K_1atm"].blockers == ()
    assert hot["flat_298.15K_1atm"].dG_assoc_kcal is None  # an auxiliary value only
    assert hot["rx1_298.15K_1atm"].blockers == ("thermo_unavailable", "spin_contaminated")


def _relabel(view, updates):
    arts = [a.model_copy(update={"payload": a.payload.model_copy(update=updates[a.artifact_id])})
            if a.artifact_id in updates else a for a in view.artifacts]
    return view.model_copy(update={"artifacts": arts})


def _lower_minima(ev, s2=None):
    """One more reactant and O minimum, 1 kcal/mol lower each (m_low_reactant, m_low_o), their
    freq's <S2> ``s2``."""
    arts = []
    for tag, composition in (("reactant", "HNO"), ("o", "o")):
        e = ev[tag].energy_hartree - 1 / HARTREE_TO_KCAL_MOL
        low = ev[tag].model_copy(update={"job_key": f"low_{tag}", "energy_hartree": e, "s2": s2})
        arts += [Artifact(artifact_id=calc_id(low), type=T.CALCULATION, payload=low),
                 Artifact(artifact_id=f"m_low_{tag}", type=T.MINIMUM, payload=R.MinimumRecord(
                     minimum_id=f"m_low_{tag}", basin_id=f"low_{tag}",
                     composition_id=composition, species_id=tag, tier="dft", level_key="x",
                     opt_calc="o", freq_calc=calc_id(low), energy_hartree=0.0,
                     state_label=LABEL.get(tag, tag)))]
    return arts


def test_state_g_is_the_lowest_spin_clean_minimum_of_the_state(fake_runtime, tmp_run):
    """dG_eff refers to the reactant state, dG_assoc compares it with the monomer states
    (Curtin-Hammett); a spin-contaminated minimum is no candidate."""
    inputs, rt, ev = _setup(fake_runtime, tmp_run)

    def with_lower_minima(s2):
        view = inputs.model_copy(update={"artifacts": [*inputs.artifacts,
                                                       *_lower_minima(ev, s2)]})
        return _thermo(view, rt)["rx1_298.15K_1atm"]

    before, clean = _thermo(inputs, rt)["rx1_298.15K_1atm"], with_lower_minima(None)
    assert clean.dG_act_kcal == pytest.approx(before.dG_act_kcal)  # seen from its own conformer
    assert clean.dG_eff_kcal == pytest.approx(before.dG_eff_kcal + 1.0)
    assert clean.dG_assoc_kcal == pytest.approx(before.dG_assoc_kcal)  # both states 1 lower
    hot = with_lower_minima(0.75)  # singlets with a doublet's <S2>
    assert (hot.dG_eff_kcal, hot.dG_assoc_kcal, hot.blockers) == (
        pytest.approx(before.dG_eff_kcal), pytest.approx(before.dG_assoc_kcal), ())


def test_a_state_g_fails_closed_without_the_layer_of_one_minimum(fake_runtime, tmp_run):
    """A minimum of the state without its energy layer may be its lowest: the state G is None
    and the reaction thermo_unavailable, though its own minima have their layer."""
    inputs, rt, ev = _setup(fake_runtime, tmp_run)
    view = inputs.model_copy(update={"artifacts": [*inputs.artifacts, *_lower_minima(ev)]})
    sp = [a.model_copy(update={"parents": tuple(p for p in a.parents if not p.startswith("m_low"))})
          for a in _sp(view, rt)]
    out = _thermo(view.model_copy(update={"artifacts": [*view.artifacts, *sp]}), rt, "tzvp")
    rx = out["rx1_298.15K_1atm"]
    assert out["m_reactant_298.15K"].G_hartree is not None
    assert "energy_layer_missing" in out["m_low_reactant_298.15K"].notes
    assert rx.blockers == ("thermo_unavailable",) and rx.dG_eff_kcal is None
    assert rx.dG_act_kcal is not None and rx.dG_assoc_kcal is None  # the O state is closed too


def test_a_submerged_barrier_drops_out_and_a_barrierless_step_has_no_dg_eff(fake_runtime,
                                                                           tmp_run):
    """Without declared monomers the chain is [R, TS, P]: a TS whose reverse dE0 is <= 0 leaves
    max(dG_rxn, 0); a barrierless step gives dG_rxn only. With them, the reverse step ends at
    the separated NH + O: their release from its product complex is its span."""
    inputs, rt, ev = _setup(fake_runtime, tmp_run)
    rx1 = inputs.get("rx1").payload
    down = rx1.model_copy(update={"reaction_id": "down", "minima": ("m_product", "m_reactant")})
    flat = [r.model_copy(update={"reaction_id": f"bl_{r.reaction_id}", "saddle": None,
                                 "outcome": R.CaseOutcome.BARRIERLESS}) for r in (rx1, down)]
    view = _with(inputs, down, *flat)
    plain = _thermo(_undeclared(view), rt)
    up, back = plain["bl_rx1_298.15K_1atm"], plain["bl_down_298.15K_1atm"]
    assert up.dG_rxn_kcal == pytest.approx(plain["rx1_298.15K_1atm"].dG_rxn_kcal)
    assert up.dG_rxn_kcal > 0 > back.dG_rxn_kcal
    assert up.dG_eff_kcal is None and back.dG_eff_kcal is None and up.dG_act_kcal is None
    assert up.blockers == () == back.blockers and plain["rx1_298.15K_1atm"].notes == ()
    # E_TS 1.6 kcal/mol above the product, whose ZPE is 2.3 kcal/mol higher than the TS's
    # (reaction mode dropped): the reverse dE0 is <= 0, the forward one stays > 0 (P4a).
    ts = calc_id(ev["ts"])
    sunk = _thermo(_undeclared(_relabel(view, {ts: {
        "energy_hartree": ev["product"].energy_hartree + 0.0025}})), rt)
    for name, value in (("rx1", up.dG_rxn_kcal), ("down", 0.0)):
        rx = sunk[f"{name}_298.15K_1atm"]
        assert rx.notes == ("submerged_barrier",) and rx.blockers == () and rx.dE_act_kcal > 0
        assert rx.dG_eff_kcal == pytest.approx(value)  # max(dG_rxn, 0)
    assert sunk["rx1_298.15K_1atm"].dG_act_kcal < sunk["rx1_298.15K_1atm"].dG_rxn_kcal
    assert sunk["down_298.15K_1atm"].dG_act_kcal < 0
    declared = _thermo(view, rt)
    assert declared["down_298.15K_1atm"].dG_eff_kcal == pytest.approx(
        -declared["rx1_298.15K_1atm"].dG_assoc_kcal)  # [R, TS, P, P_sep]: P_sep - P


def test_chiral_minimum_and_ts_gain_minus_rt_ln2(fake_runtime, tmp_run, monkeypatch):
    """The point group's m = 2 adds -RT ln 2 to G; nothing else of the group moves G here."""
    inputs, rt, ev = _setup(fake_runtime, tmp_run)
    plain = _thermo(inputs, rt)
    mirrored = [rt.load_xyz(ev[k].final).coords for k in ("ts", "product")]  # p2: the product's
    analyze = symmetry.analyze

    def chiral_group(symbols, coords, hessian, **kw):  # the TS and the product as if their
        sym = analyze(symbols, coords, hessian, **kw)  # group had no improper operation
        return dataclasses.replace(sym, m=2) if any(
            np.array_equal(coords, x) for x in mirrored) else sym

    monkeypatch.setattr(symmetry, "analyze", chiral_group)
    chiral = _thermo(inputs, rt)
    m = -R_KCAL_MOL_K * 298.15 * math.log(2)  # -0.4107 kcal/mol
    for name, shift in (("rx1_298.15K_1atm", (m, m)), ("rx2_298.15K_1atm", (m, m))):
        before, after = plain[name], chiral[name]
        assert after.dG_act_kcal - before.dG_act_kcal == pytest.approx(shift[0], abs=1e-9)
        assert after.dG_rxn_kcal - before.dG_rxn_kcal == pytest.approx(shift[1], abs=1e-9)


def test_the_electronic_term_enters_e_and_g_of_isolated_species_only(fake_runtime, tmp_run):
    """Declared levels on NH and the O(3P) atom's NIST levels: E_SO joins E (dE and G) and G_el
    replaces the spin multiplet in G of those minima; the complex, the product and the TS keep
    the multiplet alone. The point group of each subject is recorded."""
    levels = ((1, 0.0), (1, 200.0))
    plain_inputs, plain_rt, _ = _setup(fake_runtime, tmp_run)
    plain = _thermo(plain_inputs, plain_rt)
    inputs, rt, ev = _setup(fake_runtime, tmp_run / "declared", nh_levels=levels)
    out = _thermo(inputs, rt)
    T = 298.15
    nh = thermo.spin_orbit(levels) + thermo.electronic(levels, T) - thermo.electronic(((1, 0.0),), T)
    for key, shift in (("m_nh_298.15K", nh), ("m_reactant_298.15K", 0.0),
                       (f"{calc_id(ev['ts'])}_298.15K", 0.0)):
        assert out[key].G_hartree - plain[key].G_hartree == pytest.approx(shift, abs=1e-12)
    o = out["m_o_298.15K"]
    oxygen = atom_levels("O", 0, 3)
    assert (o.point_group, o.sigma, o.m) == ("Kh", 1, 1)
    assert (out["m_nh_298.15K"].point_group, out["m_nh_298.15K"].sigma) == ("Cinfv", 1)
    assert thermo.electronic(oxygen, T) * HARTREE_TO_KCAL_MOL == pytest.approx(
        -R_KCAL_MOL_K * T * math.log(6.73), abs=2e-3)
    rx, before = out["rx1_298.15K_1atm"], plain["rx1_298.15K_1atm"]
    assert rx.dE_act_kcal == pytest.approx(before.dE_act_kcal)  # complex and TS: quenched
    assert rx.dG_assoc_kcal - before.dG_assoc_kcal == pytest.approx(-nh * HARTREE_TO_KCAL_MOL)

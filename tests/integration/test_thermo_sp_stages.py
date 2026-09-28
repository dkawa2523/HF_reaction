"""sp -> thermo with FakeQM and fake_species_thermo: energy layer (parents, label, spin),
band, association, mixed LOT, m = 2, dG_eff (state reference, submerged barrier,
barrierless)."""

import math
from pathlib import Path

import numpy as np
import pytest
from fakes import PES, FakeQM, double_well, fake_species_thermo, write_geometry

from hfauto.backends.protocols import Capability
from hfauto.chemistry import thermo
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
NEUTRAL = {"complex": (0, 1), "nh": (0, 1), "o": (0, 1)}  # (charge, multiplicity)


@pytest.fixture(autouse=True)
def _closed_form_thermo(monkeypatch):
    monkeypatch.setattr(thermo, "species_thermo", fake_species_thermo)


def _state(ev, charge, multiplicity):
    level = ev.level.model_copy(update={"charge": charge, "multiplicity": multiplicity})
    return ev.model_copy(update={"level": level})


def _setup(fake_runtime, tmp_run, states=NEUTRAL):
    """The double well as composition c = NH + O (p2: the product at another grid; o2: a second
    minimum of the O state), rx1 = reactant -> product and rx2 = reactant -> p2."""
    well = double_well()
    part = well.energy(well.points["reactant"]) / 2 + 0.05  # NH and O: 63 kcal/mol up
    pes = PES(well.symbols, lambda x: well.energy(x) if np.size(x) == 9 else part, well.points)
    qm = FakeQM(tmp_run, pes)
    ev = {p: _state(qm.frequencies(pes.molecule(p), DFT), *states["complex"])
          for p in ("reactant", "product", "ts")}
    medium = ev["product"].level.model_copy(update={"grid": "medium"})
    ev["p2"] = ev["product"].model_copy(update={"job_key": "p2", "level": medium})  # other LOT
    for tag, symbols, nu, external in (("nh", ["N", "H"], (3000.0,), 5), ("o", ["O"], (), 3)):
        geom = write_geometry(tmp_run, f"{tag}.xyz", symbols, np.eye(3)[: len(symbols)])
        ev[tag] = _state(ev["reactant"], *states[tag]).model_copy(update={
            "start": geom, "final": geom, "job_key": tag, "frequencies_cm1": nu,
            "n_external": external, "energy_hartree": part})
    ev["o2"] = ev["o"].model_copy(update={"job_key": "o2"})
    arts = [Artifact(artifact_id=calc_id(e), type=T.CALCULATION, payload=e) for e in ev.values()]
    arts += [Artifact(artifact_id=t, type=T.SPECIES, payload=R.SpeciesRecord(
        species_id=t, composition_id=t, charge=states[t][0], multiplicity=states[t][1],
        geometry=ev[t].final, source="input", state_label=t)) for t in ("nh", "o")]
    arts += [Artifact(artifact_id=f"m_{k}", type=T.MINIMUM, payload=R.MinimumRecord(
        minimum_id=f"m_{k}", basin_id=k, composition_id={"nh": "nh", "o": "o", "o2": "o"}.get(
            k, "HNO"), species_id=k, tier="dft", level_key="x", opt_calc="o", freq_calc=calc_id(ev[k]),
        energy_hartree=0.0, state_label={"o2": "o"}.get(k, k))) for k in ev if k != "ts"]
    term = (R.StoichTerm(composition_id="HNO", coefficient=1),)
    rx1 = R.ReactionRecord(
        reaction_id="rx1", reactants=term, products=term, minima=("m_reactant", "m_product"),
        endpoints=("a", "b"), source="declared", outcome=R.CaseOutcome.ELEMENTARY_STEP,
        saddle=R.SaddleClaim(saddle_calc="s", freq_calc=calc_id(ev["ts"]), imag_cm1=-900.0,
                             energy_hartree=0.0))
    rx2 = rx1.model_copy(update={"reaction_id": "rx2", "minima": ("m_reactant", "m_p2")})
    arts += [Artifact(artifact_id=r.reaction_id, type=T.REACTION, payload=r) for r in (rx1, rx2)]
    system = SystemConfig(system_id="s", species=[SpeciesInput(id=t, xyz=Path(f"{t}.xyz"))
                                                  for t in ("nh", "o")],  # never read
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
    assert any({"m_product", "m_p2"} <= set(a.parents) for a in sp)  # one geometry, one sp
    view = inputs.model_copy(update={"artifacts": [*inputs.artifacts, *sp]})

    plain = _thermo(view, rt)
    rx = plain["rx1_298.15K_1atm"]
    assert rx.blockers == () and rx.band_kcal[0] < rx.dG_act_kcal < rx.band_kcal[1]
    assert rx.dG_assoc_kcal < 0 and plain["rx1_298.15K_1M"].dG_assoc_kcal == pytest.approx(
        rx.dG_assoc_kcal - 1.894, abs=1e-3)
    assert rx.dG_act_vs_separated_kcal == pytest.approx(rx.dG_assoc_kcal + rx.dG_act_kcal)
    assert rx.dG_eff_kcal == pytest.approx(rx.dG_act_kcal) and rx.notes == ()  # one conformer
    assert rx.band_kcal[0] <= rx.dG_eff_kcal <= rx.band_kcal[1]
    assert plain["rx2_298.15K_1atm"].blockers == ("mixed_level_of_theory",)
    assert rx.energy_level == plain["rx2_298.15K_1atm"].energy_level == "xfake/svp"  # no grid
    composite = _thermo(view, rt, "tzvp")
    assert composite["m_reactant_298.15K"].energy_calc  # the sp its parents name
    assert composite["rx1_298.15K_1atm"].energy_level == "xfake/tzvp"
    assert composite["rx1_298.15K_1atm"].blockers == ()

    ts = calc_id(ev["ts"])  # U8-P8: a UKS-like sp far from <S^2> = 0 blocks like a freq would
    spin = [a.model_copy(update={"payload": a.payload.model_copy(update={"s2": 0.75})})
            if ts in a.parents else a for a in view.artifacts]
    hot = _thermo(view.model_copy(update={"artifacts": spin}), rt, "tzvp")["rx1_298.15K_1atm"]
    assert hot.blockers == ("spin_contaminated",) and hot.dG_eff_kcal is not None

    # S9: the sp layer was asked for, the TS has none -> fail closed
    view = inputs.model_copy(update={"artifacts": [
        *inputs.artifacts, *(a for a in sp if ts not in a.parents)]})
    layered = _thermo(view, rt, "tzvp")
    assert layered[f"{ts}_298.15K"].G_hartree is None
    assert {"thermo_unavailable", "energy_layer_missing"} <= set(layered[f"{ts}_298.15K"].notes)
    assert "thermo_unavailable" in layered["rx1_298.15K_1atm"].blockers
    assert layered["rx1_298.15K_1atm"].dE_act_kcal is None  # no dE mixing the sp and freq LOTs
    assert layered["rx1_298.15K_1atm"].energy_level is None


def test_association_across_charge_and_spin_fails_closed(fake_runtime, tmp_run):  # S14
    anion = {"complex": (-1, 3), "nh": (-1, 1), "o": (0, 3)}  # F-.HF-like, with a triplet
    inputs, rt, _ = _setup(fake_runtime, tmp_run, anion)
    rx = _thermo(inputs, rt)["rx1_298.15K_1atm"]
    assert rx.blockers == () and rx.dG_assoc_kcal < 0 and rx.dG_act_vs_separated_kcal is not None
    sp = [a.model_copy(update={"parents": tuple(p for p in a.parents if p != "m_o2")})
          for a in _sp(inputs, rt)]  # m_o2 (the O geometry again) loses its energy layer
    view = inputs.model_copy(update={"artifacts": [*inputs.artifacts, *sp]})
    out = _thermo(view, rt, "tzvp")
    assert out["m_o2_298.15K"].G_hartree is None and out["m_o_298.15K"].G_hartree is not None
    assert out["rx1_298.15K_1atm"].dG_assoc_kcal is None  # no ensemble from m_o alone
    isomer = _relabel(view, {"m_o2": {"state_label": "o_isomer"}})  # another O state
    rt_ln2 = R_KCAL_MOL_K * 298.15 * math.log(2)  # O alone: its ensemble loses the o2 twin
    assert _thermo(isomer, rt, "tzvp")["rx1_298.15K_1atm"].dG_assoc_kcal == pytest.approx(
        rx.dG_assoc_kcal - rt_ln2)


def _relabel(view, updates):
    arts = [a.model_copy(update={"payload": a.payload.model_copy(update=updates[a.artifact_id])})
            if a.artifact_id in updates else a for a in view.artifacts]
    return view.model_copy(update={"artifacts": arts})


def test_dg_eff_takes_the_lowest_conformer_of_the_reactant_state(fake_runtime, tmp_run):
    inputs, rt, ev = _setup(fake_runtime, tmp_run)
    low = ev["reactant"].model_copy(update={
        "job_key": "low", "energy_hartree": ev["reactant"].energy_hartree - 1 / HARTREE_TO_KCAL_MOL})
    conformer = R.MinimumRecord(
        minimum_id="m_low", basin_id="low", composition_id="HNO", species_id="low", tier="dft",
        level_key="x", opt_calc="o", freq_calc=calc_id(low), energy_hartree=0.0,
        state_label="reactant")
    view = inputs.model_copy(update={"artifacts": [
        *inputs.artifacts, Artifact(artifact_id=calc_id(low), type=T.CALCULATION, payload=low),
        Artifact(artifact_id="m_low", type=T.MINIMUM, payload=conformer)]})
    before, after = _thermo(inputs, rt)["rx1_298.15K_1atm"], _thermo(view, rt)["rx1_298.15K_1atm"]
    assert after.dG_act_kcal == pytest.approx(before.dG_act_kcal)  # seen from its own conformer
    assert after.dG_eff_kcal == pytest.approx(before.dG_eff_kcal + 1.0)  # Curtin-Hammett


def test_a_submerged_barrier_and_a_barrierless_step_rank_by_max_dg_rxn_0(fake_runtime, tmp_run):
    inputs, rt, ev = _setup(fake_runtime, tmp_run)
    rx1 = inputs.get("rx1").payload
    down = rx1.model_copy(update={"reaction_id": "down", "minima": ("m_product", "m_reactant")})
    flat = [r.model_copy(update={"reaction_id": f"bl_{r.reaction_id}", "saddle": None,
                                 "outcome": R.CaseOutcome.BARRIERLESS}) for r in (rx1, down)]
    view = inputs.model_copy(update={"artifacts": [*inputs.artifacts, *(
        Artifact(artifact_id=r.reaction_id, type=T.REACTION, payload=r) for r in (down, *flat))]})
    plain = _thermo(view, rt)
    up, back = plain["bl_rx1_298.15K_1atm"], plain["bl_down_298.15K_1atm"]
    assert up.dG_eff_kcal == pytest.approx(plain["rx1_298.15K_1atm"].dG_rxn_kcal)
    assert up.dG_eff_kcal > 0 and up.dG_act_kcal is None
    assert back.dG_eff_kcal == 0.0 and back.dG_rxn_kcal < 0
    assert plain["rx1_298.15K_1atm"].notes == ()
    # E_TS 1.6 kcal/mol above the product, whose ZPE is 2.3 kcal/mol higher than the TS's
    # (reaction mode dropped): the reverse dE0 is <= 0, the forward one stays > 0 (P4a).
    ts = calc_id(ev["ts"])
    sunk = _thermo(_relabel(view, {ts: {"energy_hartree": ev["product"].energy_hartree
                                        + 0.0025}}), rt)
    for name, ref in (("rx1", up), ("down", back)):
        rx = sunk[f"{name}_298.15K_1atm"]
        assert rx.notes == ("submerged_barrier",) and rx.blockers == () and rx.dE_act_kcal > 0
        assert rx.dG_eff_kcal == pytest.approx(ref.dG_eff_kcal)  # max(dG_rxn, 0)
    assert sunk["rx1_298.15K_1atm"].dG_act_kcal < sunk["rx1_298.15K_1atm"].dG_rxn_kcal
    assert sunk["down_298.15K_1atm"].dG_act_kcal < 0


def test_chiral_minimum_and_ts_gain_minus_rt_ln2(fake_runtime, tmp_run):
    inputs, rt, ev = _setup(fake_runtime, tmp_run)
    plain = _thermo(inputs, rt)
    chfclbr = [[0.0, 0.0, 0.0], [0.63, 0.63, 0.63], [-0.8, -0.8, 0.8], [-1.0, 1.0, -1.0],
               [1.1, -1.1, -1.1]]
    ts = write_geometry(tmp_run, "chiral_ts.xyz", ["C", "H", "F", "Cl", "Br"], np.array(chfclbr))
    updates = {calc_id(ev["ts"]): {"final": ts},  # a C1 TS geometry; its thermo is unchanged
               "m_product": {"chiral": True}}
    chiral = _thermo(_relabel(inputs, updates), rt)
    m = -R_KCAL_MOL_K * 298.15 * math.log(2)  # -0.4107 kcal/mol
    for name, shift in (("rx1_298.15K_1atm", (m, m)), ("rx2_298.15K_1atm", (m, 0.0))):
        before, after = plain[name], chiral[name]  # rx2 ends at m_p2, which is achiral
        assert after.dG_act_kcal - before.dG_act_kcal == pytest.approx(shift[0], abs=1e-9)
        assert after.dG_rxn_kcal - before.dG_rxn_kcal == pytest.approx(shift[1], abs=1e-9)

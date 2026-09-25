"""sp -> thermo with FakeQM / FakeThermo: energy layer, band, association, mixed LOT."""

from dataclasses import replace

import numpy as np
import pytest
from fakes import FakeQM, FakeThermo, double_well, write_geometry

from hfauto.backends.protocols import Capability, ThermoResult
from hfauto.chemistry.gates import zpe_hartree
from hfauto.chemistry.thermo import settings_sha
from hfauto.core import records as R
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec
from hfauto.core.system import CompositionInput, Conditions, SpeciesInput, SystemConfig
from hfauto.drivers.minimum import calc_id
from hfauto.stages.single_point import SinglePointConfig, SinglePointStage
from hfauto.stages.thermochemistry import ThermoConfig, ThermoStage

T = R.ArtifactType
DFT, BIG = (MethodSpec(id=b, kind="dft", functional="xfake", basis=b) for b in ("svp", "tzvp"))


def _gv(freq, settings):  # consistent stand-in: G = E + ZPE - 1e-6 x cutoff x n_real
    e, zpe, n = freq.energy_hartree, zpe_hartree(nu := freq.frequencies_cm1), sum(f > 0 for f in nu)
    return [ThermoResult(settings_sha=settings_sha(s), T_K=t, E_hartree=e, H_hartree=e, S_rot=1,
                         G_hartree=e + zpe - 1e-6 * s.cutoff_cm1 * n, zpe_hartree=zpe, n_real=n,
                         notes=(), job_key="k") for s in settings for t in s.temperatures_K]


def test_sp_then_thermo(fake_runtime, tmp_run):
    qm = FakeQM(tmp_run, pes := double_well())
    ev = {p: qm.frequencies(pes.molecule(p), DFT) for p in ("reactant", "product", "ts")}
    medium = ev["product"].level.model_copy(update={"grid": "medium"})
    ev["p2"] = ev["product"].model_copy(update={"job_key": "p2", "level": medium})  # other LOT
    for tag, symbols, nu in (("nh", ["N", "H"], (3000.0,)), ("o", ["O"], ())):  # 63 kcal/mol up
        geom = write_geometry(tmp_run, f"{tag}.xyz", symbols, np.eye(3)[: len(symbols)])
        ev[tag] = ev["reactant"].model_copy(update={
            "start": geom, "final": geom, "job_key": tag, "frequencies_cm1": nu,
            "energy_hartree": ev["reactant"].energy_hartree / 2 + 0.05})
    arts = [Artifact(artifact_id=calc_id(e), type=T.CALCULATION, payload=e) for e in ev.values()]
    arts += [Artifact(artifact_id=t, type=T.SPECIES, payload=R.SpeciesRecord(
        species_id=t, composition_id=t, formula=t, charge=0, multiplicity=1,
        geometry=ev[t].final, source="input", state_label=t)) for t in ("nh", "o")]
    arts += [Artifact(artifact_id=f"m_{k}", type=T.MINIMUM, payload=R.MinimumRecord(
        minimum_id=f"m_{k}", basin_id=k, composition_id={"nh": "nh", "o": "o"}.get(k, "HNO"),
        species_id=k, tier="dft", level_key="x", opt_calc="o", freq_calc=calc_id(ev[k]),
        energy_hartree=0.0, state_label=k, n_fragments=1)) for k in ev if k != "ts"]
    term = (R.StoichTerm(composition_id="HNO", coefficient=1),)
    rx1 = R.ReactionRecord(
        reaction_id="rx1", reactants=term, products=term, minima=("m_reactant", "m_product"),
        endpoints=("a", "b"), source="declared", outcome=R.CaseOutcome.ELEMENTARY_STEP,
        saddle=R.SaddleClaim(saddle_calc="s", freq_calc=calc_id(ev["ts"]), imag_cm1=-900.0,
                             energy_hartree=0.0))
    rx2 = rx1.model_copy(update={"reaction_id": "rx2", "minima": ("m_reactant", "m_p2")})
    arts += [Artifact(artifact_id=r.reaction_id, type=T.REACTION, payload=r) for r in (rx1, rx2)]
    system = SystemConfig(system_id="s", species=[SpeciesInput(id="nh"), SpeciesInput(id="o")],
                          compositions=[CompositionInput(id="c", components={"nh": 1, "o": 1})])
    rt = replace(fake_runtime(system, {(Capability.QM, "nwchem"): qm, (
        Capability.THERMO, "goodvibes"): FakeThermo(_gv)}, methods={"svp": DFT, "tzvp": BIG}),
        conditions=Conditions(standard_states=("1atm", "1M")))
    inputs = Manifest(run_id="r", stage_id="v", created_at="t", artifacts=arts)

    sp = SinglePointStage().run(inputs, SinglePointConfig(engine="nwchem", methods=["tzvp"]), rt)
    assert sorted((a.payload.task, len(a.parents)) for a in sp) == [("sp", 1)] * 2 + [("sp", 2)]
    view = inputs.model_copy(update={"artifacts": [*arts, *sp]})  # p2 shares product's sp

    def thermo(method):
        out = ThermoStage().run(view, ThermoConfig(engine="goodvibes", energy_method=method), rt)
        return {a.artifact_id.split("_thermo_")[1]: a.payload for a in out}  # <subject|rxn>_<T>...

    plain = thermo(None)
    rx = plain["rx1_298.15K_1atm"]
    assert rx.blockers == () and rx.band_kcal[0] < rx.dG_act_kcal < rx.band_kcal[1]
    assert rx.dG_assoc_kcal < 0 and plain["rx1_298.15K_1M"].dG_assoc_kcal == pytest.approx(
        rx.dG_assoc_kcal - 1.894, abs=1e-3)
    assert rx.dG_act_vs_separated_kcal == pytest.approx(rx.dG_assoc_kcal + rx.dG_act_kcal)
    assert plain["rx2_298.15K_1atm"].blockers == ("mixed_level_of_theory",)
    pop = {k: plain[f"m_{k}_298.15K"].population for k in ("reactant", "product", "p2")}
    assert pop["reactant"] + pop["product"] == pytest.approx(1.0) and pop["p2"] == 1.0  # per LOT
    assert thermo("tzvp")["m_reactant_298.15K"].energy_calc  # composite G on the sp layer

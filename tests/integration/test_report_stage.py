import csv

from fakes import FakeQM, double_well

from hfauto.core import records as rec
from hfauto.core.evidence import Failure, FailureKind
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec
from hfauto.core.system import SpeciesInput, SystemConfig
from hfauto.stages.report import ReportConfig, ReportStage

T = rec.ArtifactType


def test_report_stage_writes_tables_from_a_fake_view(fake_runtime, tmp_run):
    rt = fake_runtime(SystemConfig(system_id="s", species=[SpeciesInput(id="w", smiles="O")]),
                      {}, stage_id="report")
    pes = double_well()
    qm, arts = FakeQM(tmp_run, pes), []
    for method in (MethodSpec(id="pbe0", kind="dft", functional="pbe0", grid="fine"),
                   MethodSpec(id="wb97", kind="dft", functional="wb97x", basis="def2-tzvpd")):
        arts += [Artifact(artifact_id=f"{method.id}_{p}", type=T.CALCULATION,
                          payload=qm.energy(pes.molecule(p), method))
                 for p in ("reactant", "ts", "product")]
    arts += [Artifact(artifact_id=side, type=T.MINIMUM, payload=rec.MinimumRecord(
        minimum_id=side, basin_id=side, composition_id="c", species_id=side, tier="dft",
        level_key="k", opt_calc=f"pbe0_{side}", freq_calc=f"pbe0_{side}", energy_hartree=0.0,
        state_label=side)) for side in ("reactant", "product")]
    saddle = rec.SaddleClaim(saddle_calc="pbe0_ts", freq_calc="pbe0_ts", imag_cm1=-800.0,
                             energy_hartree=0.0)
    reaction = rec.ReactionRecord(reaction_id="r1", reactants=(), products=(), saddle=saddle,
                                  minima=("reactant", "product"), endpoints=("a", "b"),
                                  source="declared", outcome=rec.CaseOutcome.ELEMENTARY_STEP)
    thermo = rec.ReactionThermo(reaction_id="r1", T_K=298.15, standard_state="1atm",
                                dE_act_kcal=6.0, dE_rxn_kcal=2.0, dzpe_act_kcal=-1.0,
                                dG_act_kcal=5.0, dG_rxn_kcal=2.0, band_kcal=(4.5, 5.5),
                                dG_eff_kcal=5.0)
    found = rec.DiscoveryRecord(discovery_id="d1", source_minimum="reactant", mechanism="nt2",
                                outcome="product")
    failure = Failure(kind=FailureKind.SCF_NOT_CONVERGED, reason="scf")
    arts += [Artifact(artifact_id="r1", type=T.REACTION, payload=reaction),
             Artifact(artifact_id="r1_thermo", type=T.REACTION_THERMO, payload=thermo),
             Artifact(artifact_id="d1", type=T.DISCOVERY, payload=found),
             Artifact(artifact_id="bad", type=T.DISCOVERY, status="failed", failure=failure)]
    view = Manifest(run_id="run", stage_id="report", created_at="now", artifacts=arts)

    [artifact] = ReportStage().run(view, ReportConfig(), rt)
    rows = artifact.payload.rows
    assert artifact.type == T.REPORT and [(r.reaction_id, r.rank, r.tier) for r in rows] == [
        ("r1", 1, "saddle")]
    tables = artifact.payload.tables
    assert sorted(tables) == ["coverage.csv", "method_panel.csv", "ranking.csv", "report.html"]
    page = (tmp_run / tables["report.html"].path).read_text(encoding="utf-8")
    assert "<svg class='energy'" in page and "href='#rxn-r1'" in page and "800.0i" in page
    read = {name: list(csv.DictReader((tmp_run / ref.path).read_text().splitlines()))
            for name, ref in tables.items() if name.endswith(".csv")}
    panel = read["method_panel.csv"]  # two levels on one PES: the same values
    assert len({row["level_key"] for row in panel}) == 2
    assert all(float(row["dE_rxn_kcal"]) > 0 for row in panel)  # the product lies higher
    assert {(r["metric"], r["key"], r["count"]) for r in read["coverage.csv"]} == {
        ("attempts", "nt2", "1"), ("products", "nt2", "1"),
        ("failure_kind", "scf_not_converged", "1")}
    ranked = read["ranking.csv"][0]
    assert (ranked["dG_eff_kcal"], ranked["band_low_kcal"], ranked["dG_act_kcal"]) == (
        "5.0", "4.5", "5.0")
    assert (ranked["T_K"], ranked["standard_state"], ranked["dG_rxn_kcal"]) == (
        "298.15", "1atm", "2.0")
    assert (ranked["torsional"], ranked["dG_act_vs_separated_kcal"], ranked["notes"]) == (
        "False", "", "")
    assert ranked["dE_act_panel_min_kcal"] == panel[0]["dE_act_min_kcal"]

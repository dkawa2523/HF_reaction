import csv
import re

from fakes import FakeQM, double_well

from hfauto.core import records as rec
from hfauto.core.evidence import Failure, FailureKind
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec
from hfauto.core.system import SpeciesInput, SystemConfig
from hfauto.stages.report import ReportConfig, ReportStage

T = rec.ArtifactType


def test_report_stage_writes_tables_from_a_fake_view(fake_runtime, tmp_run):
    rt = fake_runtime(SystemConfig(system_id="s", species=[SpeciesInput(id="w", smiles="O",
                                                                multiplicity=1)]),
                      {}, stage_id="report")
    pes = double_well()
    qm, arts = FakeQM(tmp_run, pes), []
    subjects = {"reactant": ("reactant",), "ts": ("pbe0_ts",), "product": ("product",)}
    for method in (MethodSpec(id="pbe0", kind="dft", functional="pbe0", grid="fine"),
                   MethodSpec(id="wb97", kind="dft", functional="wb97x", basis="def2-tzvpd")):
        arts += [Artifact(artifact_id=f"{method.id}_{p}", type=T.CALCULATION,
                          payload=qm.energy(pes.molecule(p), method),  # sp parents: subjects
                          parents=subjects[p] if method.id == "wb97" else ())
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
                                dE_act_kcal=6.0, dE_rxn_kcal=2.0,
                                dG_act_kcal=5.0, dG_rxn_kcal=2.0, band_kcal=(4.5, 5.5),
                                dG_eff_kcal=5.0, energy_level="wb97x/def2-tzvpd")
    hot = thermo.model_copy(update={"T_K": 400.0, "dG_eff_kcal": 7.0})  # a second temperature
    found = rec.DiscoveryRecord(discovery_id="d1", source_minimum="reactant", mechanism="nt2",
                                outcome="product")
    failure = Failure(kind=FailureKind.SCF_NOT_CONVERGED, reason="scf")
    arts += [Artifact(artifact_id="r1", type=T.REACTION, payload=reaction),
             Artifact(artifact_id="r1_thermo", type=T.REACTION_THERMO, payload=thermo),
             Artifact(artifact_id="r1_thermo_400", type=T.REACTION_THERMO, payload=hot),
             Artifact(artifact_id="d1", type=T.DISCOVERY, payload=found),
             Artifact(artifact_id="bad", type=T.DISCOVERY, status="failed", failure=failure)]
    view = Manifest(run_id="run", stage_id="report", created_at="now", artifacts=arts)

    [artifact] = ReportStage().run(view, ReportConfig(), rt)  # the first (T, state) of thermo
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
    assert list(read["ranking.csv"][0]) == [
        "rank", "reaction_id", "outcome", "tier", "rankable", "T_K", "standard_state",
        "energy_level", "dG_eff_kcal", "band_low_kcal", "band_high_kcal", "dG_act_kcal",
        "dG_rxn_kcal", "dG_act_vs_separated_kcal", "torsional", "blockers", "notes",
        "dE_act_panel_min_kcal", "dE_act_panel_max_kcal"]
    ranked = read["ranking.csv"][0]
    assert ranked["energy_level"] == "wb97x/def2-tzvpd"
    assert (ranked["dG_eff_kcal"], ranked["band_low_kcal"], ranked["dG_act_kcal"]) == (
        "5.0", "4.5", "5.0")
    assert (ranked["T_K"], ranked["standard_state"], ranked["dG_rxn_kcal"]) == (
        "298.15", "1atm", "2.0")
    assert (ranked["torsional"], ranked["dG_act_vs_separated_kcal"], ranked["notes"]) == (
        "False", "", "")
    assert ranked["dE_act_panel_min_kcal"] == panel[0]["dE_act_min_kcal"]
    [at_400] = ReportStage().run(view, ReportConfig(T_K=400.0), rt)
    assert at_400.payload.rows[0].dG_eff_kcal == 7.0
    page = (tmp_run / at_400.payload.tables["report.html"].path).read_text(encoding="utf-8")
    cells = re.findall(r"<td[^>]*>(.*?)</td>", page.split("</table>")[0])  # the summary row
    assert (cells[3], cells[-1]) == ("7.00", "wb97x/def2-tzvpd")  # at the ranked (T, state)
    assert "Energy diagram (400 K, 1atm)" in page

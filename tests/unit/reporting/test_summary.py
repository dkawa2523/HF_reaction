import csv
from dataclasses import astuple

from hfauto.core import records as rec
from hfauto.core.evidence import Evidence, Failure, FailureKind, FileRef, Geometry, Level
from hfauto.core.manifest import Artifact
from hfauto.reporting.summary import coverage, method_panel, rank_rows, write_tables

EL = rec.CaseOutcome.ELEMENTARY_STEP


def rxn(rid, outcome=EL, minima=("mr", "mp"), **kw):
    return rec.ReactionRecord(reaction_id=rid, reactants=(), products=(), minima=minima,
                              endpoints=("a", "b"), source="declared", outcome=outcome, **kw)


def th(rid, dg, band=None, T=298.15, blockers=()):
    return rec.ReactionThermo(reaction_id=rid, T_K=T, standard_state="1atm", dE_act_kcal=dg,
                              dE_rxn_kcal=1.0, dzpe_act_kcal=-0.5, dG_act_kcal=dg,
                              dG_rxn_kcal=1.0, dG_eff_kcal=dg, band_kcal=band, blockers=blockers,
                              energy_level="wb97x-d3/def2-tzvpd")


def calc(cid, fp, energy, grid, *subjects):  # an sp names its subjects in its parents
    level = Level(program="nwchem", version="7.2.3", method="pbe0", basis="def2-svpd",
                  charge=0, multiplicity=1, grid=grid)
    g = Geometry(file=FileRef(path=f"{fp}.xyz", sha256="0"), fingerprint=fp, symbols=("H",))
    freq = {} if subjects else {"frequencies_cm1": (), "n_external": 3}  # an atom: no mode
    payload = Evidence(engine="nwchem", task="freq" if not subjects else "sp", level=level,
                       start=g, final=g, energy_hartree=energy, output=g.file, job_key=cid,
                       **freq)
    return Artifact(artifact_id=cid, type=rec.ArtifactType.CALCULATION, payload=payload,
                    parents=subjects)


def mini(mid, calc):
    return rec.MinimumRecord(minimum_id=mid, basin_id=mid, composition_id="c", species_id=mid,
                             tier="dft", level_key="k", opt_calc=calc, freq_calc=calc,
                             energy_hartree=0.0, state_label="l")


def test_ranks_follow_dg_eff_with_overlapping_bands_sharing_a_rank():
    reactions = [rxn("a"), rxn("b"), rxn("c"), rxn("d", torsional=True),
                 rxn("bl", rec.CaseOutcome.BARRIERLESS), rxn("spin"), rxn("open", None)]
    thermo = [th("a", 10.0, (9.0, 11.0)), th("b", 10.5, (10.8, 12.0)), th("c", 20.0, (19.0, 21.0)),
              th("d", 12.5, (12.1, 13.0)), th("bl", 1.0), th("spin", 2.0, None, 298.15,
                                                             ("spin_contaminated",)),
              th("c", 1.0, T=500.0)]
    rows = rank_rows(reactions, thermo, 298.15, "1atm")
    assert [(r.reaction_id, r.rank) for r in rows] == [  # a barrierless step ranks too
        ("bl", 1), ("a", 2), ("b", 2), ("d", 4), ("c", 5), ("spin", None)]
    by_id = {r.reaction_id: r for r in rows}
    assert by_id["bl"].blockers == () and by_id["spin"].blockers == ("spin_contaminated",)
    assert by_id["a"].tier == "minima" and by_id["a"].band_kcal == (9.0, 11.0)
    assert by_id["d"].torsional and not by_id["a"].torsional  # a torsion stays ranked
    assert {(r.T_K, r.standard_state, r.dG_rxn_kcal) for r in rows} == {(298.15, "1atm", 1.0)}
    assert {r.energy_level for r in rows} == {"wb97x-d3/def2-tzvpd"}
    assert [r.rank for r in rank_rows(reactions, thermo, 298.15, "1M")] == [None] * 6
    new = {"T_K", "standard_state", "dG_rxn_kcal", "dG_eff_kcal", "dG_act_vs_separated_kcal",
           "torsional", "notes"}
    old = rows[0].model_dump(exclude=new)
    assert rec.RankRow.model_validate(old).T_K is None  # a report row of an older manifest


def test_coverage_counts_mechanisms_negative_reasons_and_failure_kinds():
    def disc(i, mech, outcome, reason=None):
        return rec.DiscoveryRecord(discovery_id=f"d{i}", source_minimum="m", mechanism=mech,
                                   outcome=outcome, reason=reason)

    found = [disc(1, "nt2", "product"), disc(2, "nt2", "negative", "collapsed_to:x_1"),
             disc(3, "nt2", "negative", "collapsed_to:y_2"), disc(4, "afir", "failed"),
             disc(5, "afir", "negative", "out_of_window")]
    failure = Failure(kind=FailureKind.TIMEOUT, reason="walltime")
    failed = [Artifact(artifact_id=f"f{i}", type=rec.ArtifactType.CALCULATION, status="failed",
                       failure=failure) for i in range(2)]
    assert {astuple(row) for row in coverage(found, failed)} == {
        ("attempts", "nt2", 3), ("products", "nt2", 1), ("negatives", "nt2", 2),
        ("attempts", "afir", 2), ("failed", "afir", 1), ("negatives", "afir", 1),
        ("negative_reason", "collapsed_to", 2), ("negative_reason", "out_of_window", 1),
        ("failure_kind", "timeout", 2)}


def test_panel_spread_is_two_ranking_columns_not_a_blocker(tmp_path):
    """dE_rxn changes sign between the levels (+6.28 / -6.28): no blocker, both still rank.
    Reference energies come from the freq calculations, the others through the sp parents
    (geometry fingerprints play no part)."""
    calcs = [calc("r", "R", -1.0, "fine"), calc("p", "P", -0.99, "fine"),
             calc("t", "T", -0.98, "fine"), calc("r2", "R", -1.0, "xfine", "mr"),
             calc("p2", "Q", -1.01, "xfine", "mp"), calc("t3", "T", 0, "xfine"),
             calc("t2", "T", -0.975, "xfine", "t"), calc("x", "X", 0, "xfine", "stranger")]
    minima = {"mr": mini("mr", "r"), "mp": mini("mp", "p")}
    saddle = rec.SaddleClaim(saddle_calc="s", freq_calc="t", imag_cm1=-900.0, energy_hartree=-0.98)
    reactions = [rxn("x", saddle=saddle), rxn("lost", minima=("mr", ""))]
    panel = method_panel(calcs, reactions, minima=minima)
    assert len({r.level_key for r in panel}) == 2 and {r.reaction_id for r in panel} == {"x"}
    fine, xfine = sorted(panel, key=lambda r: r.level)
    assert fine.level.endswith(" fine") and xfine.level.endswith(" xfine")
    assert (round(fine.dE_rxn_kcal, 2), round(fine.dE_act_kcal, 2)) == (6.28, 12.55)
    assert (round(xfine.dE_rxn_min_kcal, 2), round(xfine.dE_rxn_max_kcal, 2)) == (-6.28, 6.28)
    rows = rank_rows(reactions, [th("x", 5.0), th("lost", 6.0)], 298.15, "1atm")
    assert [(r.reaction_id, r.rank, r.blockers) for r in rows] == [("x", 1, ()), ("lost", 2, ())]

    def ranking(panel_rows):
        path = write_tables(tmp_path, rows, [], panel_rows)["ranking.csv"]
        return {r["reaction_id"]: r for r in csv.DictReader(path.read_text().splitlines())}

    x = ranking(panel)["x"]
    assert x["energy_level"] == "wb97x-d3/def2-tzvpd"
    assert (x["dE_act_panel_min_kcal"], x["dE_act_panel_max_kcal"]) == (
        str(fine.dE_act_min_kcal), str(fine.dE_act_max_kcal))
    assert round(float(x["dE_act_panel_max_kcal"]), 2) == 15.69
    assert ranking(panel)["lost"]["dE_act_panel_min_kcal"] == ""
    assert "dE_act_panel_min_kcal" not in ranking([fine])["x"]  # one level is no panel

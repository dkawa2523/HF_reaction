from dataclasses import astuple

from hfauto.core import records as rec
from hfauto.core.evidence import Evidence, Failure, FailureKind, FileRef, Geometry, Level
from hfauto.core.manifest import Artifact
from hfauto.reporting.summary import coverage, method_panel, participant_notes, rank_rows

EL = rec.CaseOutcome.ELEMENTARY_STEP


def rxn(rid, outcome=EL, minima=("mr", "mp"), **kw):
    return rec.ReactionRecord(reaction_id=rid, reactants=(), products=(), minima=minima,
                          endpoints=("a", "b"), source="declared", outcome=outcome, **kw)


def th(rid, dg, band=None, T=298.15):
    return rec.ReactionThermo(reaction_id=rid, T_K=T, standard_state="1atm", dE_act_kcal=dg,
                          dE_rxn_kcal=1.0, dzpe_act_kcal=-0.5, dG_act_kcal=dg, dG_rxn_kcal=1.0,
                          band_kcal=band)


def ev(fp, energy, grid):
    level = Level(program="nwchem", version="7.2.3", method="pbe0", basis="def2-svpd",
                  charge=0, multiplicity=1, grid=grid)
    g = Geometry(file=FileRef(path=f"{fp}.xyz", sha256="0"), fingerprint=fp, symbols=("H",))
    return Evidence(engine="nwchem", task="sp", level=level, start=g, final=g,
                    energy_hartree=energy, output=g.file, job_key=fp)


def mini(mid, calc, notes=()):
    return rec.MinimumRecord(minimum_id=mid, basin_id=mid, composition_id="c", species_id=mid,
                         tier="dft", level_key="k", opt_calc=calc, freq_calc=calc,
                         energy_hartree=0.0, state_label="l", n_fragments=1, notes=notes)


def test_overlapping_bands_share_a_rank_and_unrankable_reactions_are_listed_unranked():
    reactions = [rxn("a"), rxn("b"), rxn("c"), rxn("d"), rxn("bl", rec.CaseOutcome.BARRIERLESS),
                 rxn("spin"), rxn("open", None)]
    thermo = [th("a", 10.0, (9.0, 11.0)), th("b", 10.5, (10.8, 12.0)), th("c", 20.0, (19.0, 21.0)),
              th("d", 12.5, (12.1, 13.0)), th("bl", 1.0), th("spin", 2.0), th("c", 1.0, T=500.0)]
    notes = {"spin": ("spin_contaminated",)}
    rows = rank_rows(reactions, thermo, notes, 298.15, "1atm")
    assert [(r.reaction_id, r.rank) for r in rows] == [
        ("a", 1), ("b", 1), ("d", 3), ("c", 4), ("bl", None), ("spin", None)]
    by_id = {r.reaction_id: r for r in rows}
    assert by_id["bl"].blockers == ("outcome:barrierless_at_resolution",)
    assert "spin_contaminated" in by_id["spin"].blockers and not by_id["spin"].rankable
    assert by_id["a"].tier == "minima" and by_id["a"].band_kcal == (9.0, 11.0)
    assert [r.rank for r in rank_rows(reactions, thermo, notes, 298.15, "1M")] == [None] * 6
    by_rxn = rank_rows(reactions[:3], thermo[:3], {}, 298.15, "1atm", metric="dG_rxn")
    assert [r.rank for r in by_rxn] == [1, 1, 1]  # equal dG_rxn, bands apply to dG_act only


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


def test_panel_keys_by_full_level_and_sign_disagreement_blocks_ranking():
    calcs = {"r": ev("R", -1.0, "fine"), "p": ev("P", -0.99, "fine"), "t": ev("T", -0.98, "fine"),
             "r2": ev("R", -1.0, "xfine"), "p2": ev("P", -1.01, "xfine"), "x": ev("X", 0, "fine")}
    minima = {"mr": mini("mr", "r"), "mp": mini("mp", "p", ("soft_imaginary_mode",))}
    saddle = rec.SaddleClaim(saddle_calc="t", freq_calc="t", imag_cm1=-900.0, energy_hartree=-0.98)
    reactions = [rxn("x", saddle=saddle), rxn("lost", minima=("mr", ""))]
    rows, disagreements = method_panel(calcs, reactions, minima=minima)
    assert len({r.level_key for r in rows}) == 2 and {r.reaction_id for r in rows} == {"x"}
    fine, xfine = sorted(rows, key=lambda r: r.level)
    assert fine.level.endswith(" fine") and xfine.level.endswith(" xfine")
    assert round(fine.dE_rxn_kcal, 2) == 6.28 and round(fine.dE_act_kcal, 2) == 12.55
    assert xfine.dE_act_kcal is None and fine.dE_rxn_min_kcal == xfine.dE_rxn_kcal < 0
    assert disagreements == {"x"} and all(r.sign_disagreement for r in rows)
    notes = participant_notes(reactions, minima, disagreements)
    assert notes["x"] == ("soft_imaginary_mode", "method_sign_disagreement")
    lost, x = rank_rows(reactions, [th("x", 5.0), th("lost", 6.0)], notes, 298.15, "1atm")
    assert (lost.rank, x.rank, x.blockers) == (1, None, ("method_sign_disagreement",))

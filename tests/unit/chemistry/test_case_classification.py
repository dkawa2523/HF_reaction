"""finalize() and split() (design §7.3); split is ported from tests/test_reaction_segments.py."""

from __future__ import annotations

import pytest

from hfauto.chemistry.classification import finalize, split
from hfauto.core import records as r
from hfauto.core.evidence import FileRef, Geometry
from hfauto.drivers.reaction_case.state import Action, Decision

TERM = r.StoichTerm(composition_id="H3N_q0_m1", coefficient=1)
CASE = r.ReactionRecord(reaction_id="rxn", reactants=(TERM,), products=(TERM,), minima=("ma", "mb"),
                        endpoints=("sa", "sb"), source="declared", torsional=True,
                        reasons=("note",))
SADDLE = r.SaddleClaim(saddle_calc="s", freq_calc="f", imag_cm1=-900.0, energy_hartree=-56.0)
QRC = r.ConnectionClaim(side_calcs=("p", "m"), minima=("mb", "ma"), amplitude_A=0.1)


def done(outcome: r.CaseOutcome, reason: str = "r") -> Decision:
    return Decision(Action.COMPLETE, reason, outcome)


def test_finalize_attaches_outcome_reasons_and_claims():
    barrier = r.BarrierVerdict(verdict="single", source="screen", max_rel_kcal=12.0)
    rec = finalize(CASE, done(r.CaseOutcome.ELEMENTARY_STEP, "connection:elementary"),
                   barrier=barrier, claim=SADDLE, connection=QRC)
    assert rec.outcome is r.CaseOutcome.ELEMENTARY_STEP
    assert rec.reasons == ("note", "connection:elementary")
    assert (rec.barrier, rec.saddle, rec.connection) == (barrier, SADDLE, QRC)
    assert rec.minima == CASE.minima and rec.source == "declared"


def test_finalize_refuses_unsupported_conclusions():
    with pytest.raises(ValueError, match="not terminal"):
        finalize(CASE, Decision(Action.SCREEN, "screen"), barrier=None, claim=None,
                 connection=None)
    with pytest.raises(ValueError, match="claims"):
        finalize(CASE, done(r.CaseOutcome.ELEMENTARY_STEP), barrier=None, claim=SADDLE,
                 connection=None)
    rec = finalize(CASE, done(r.CaseOutcome.UNRESOLVED, "connection_failed"), barrier=None,
                   claim=SADDLE, connection=None)
    assert rec.outcome is r.CaseOutcome.UNRESOLVED and rec.saddle == SADDLE


def test_reassigned_step_adopts_the_connected_minima():
    qrc = QRC.model_copy(update={"minima": ("mc", "ma")})
    rec = finalize(CASE, done(r.CaseOutcome.REASSIGNED), barrier=None, claim=SADDLE,
                   connection=qrc)
    assert (rec.source, rec.minima, rec.reaction_id) == ("reassigned", ("ma", "mc"), "rxn")


def test_split_inherits_the_parent_stoichiometry():
    parent = finalize(CASE, done(r.CaseOutcome.MULTI_STEP, "intermediate_distinct"),
                      barrier=None, claim=None, connection=None)
    well = r.MinimumRecord(minimum_id="mi", basin_id="bi", composition_id="H3N_q0_m1",
                           species_id="si", tier="dft", level_key="L", opt_calc="o",
                           freq_calc="f", energy_hartree=-56.0, state_label="x")
    geo = Geometry(file=FileRef(path="si.xyz", sha256="0" * 64), fingerprint="si",
                   symbols=("N", "H", "H", "H"))
    well_species = r.SpeciesRecord(species_id="si", composition_id="H3N_q0_m1", charge=0,
                                   multiplicity=1, geometry=geo, source="intermediate",
                                   state_label="x")

    first, second = split(parent, well, well_species)
    assert (first.minima, second.minima) == (("ma", "mi"), ("mi", "mb"))
    assert (first.endpoints, second.endpoints) == (("sa", "si"), ("si", "sb"))
    assert (first.reaction_id, second.reaction_id) == ("rxn_split1", "rxn_split2")
    for child in (first, second):
        assert child.source == "split"
        assert child.reactants == parent.reactants and child.products == parent.products
        assert child.outcome is None and child.reasons == ()
        assert child.torsional
    other = well.model_copy(update={"composition_id": "HF_q0_m1"})
    with pytest.raises(ValueError, match="composition"):
        split(parent, other, well_species)

"""decide(): the 9 rows of design §7.3 and its one stop rule (every search without an answer
is one failure token)."""

from __future__ import annotations

from dataclasses import replace

import pytest

from hfauto.core import records as r
from hfauto.core.evidence import FileRef, Geometry
from hfauto.drivers.reaction_case import state as st
from hfauto.drivers.reaction_case.state import Action, CaseState, Decision, decide, record_profile

A, C, BV = Action, r.CaseOutcome, r.BarrierVerdict
TERM = r.StoichTerm(composition_id="CHN_q0_m1", coefficient=1)
CASE = r.ReactionRecord(reaction_id="rxn", reactants=(TERM,), products=(TERM,), minima=("ma", "mb"),
                        endpoints=("sa", "sb"), source="declared")
DEGENERATE = CASE.model_copy(update={"minima": ("ma", "ma"), "degenerate": True})
ASSOCIATION = CASE.model_copy(update={"monomers": ("m_ch3", "m_o2")})  # X4: from the monomers
GEO = Geometry(file=FileRef(path="seed.xyz", sha256="0" * 64), fingerprint="f", symbols=("H",))
SEED = st.Seed(geometry=GEO, source="screen_ts")
OTHER_TS = st.Seed(geometry=GEO, source="discovery_ts")
CONTINUED = st.Seed(geometry=GEO, source="continuation", depth=1)
CLAIM = r.SaddleClaim(saddle_calc="s", freq_calc="f", imag_cm1=-1131.0, energy_hartree=-93.0)
RULES = st.CaseRules()
BASE = CaseState()  # nothing done yet
S = replace(BASE, screen=BV(verdict="single", source="screen"), neb_done=True)  # one peak
SHORTCUT = replace(BASE, shortcut_done=True, attempts=1)  # the seed of a low-level TS failed
STRING = BV(verdict="single", source="string")
DFT_PATH = Decision(A.FIND_PATH, "dft_path")
VALIDATE = Decision(A.VALIDATE_AND_CONNECT, "saddle_converged")
EXHAUSTED = Decision(A.COMPLETE, "attempts_exhausted", C.UNRESOLVED)
ELEMENTARY = Decision(A.COMPLETE, "connection:elementary", C.ELEMENTARY_STEP)
DISTINCT = Decision(A.COMPLETE, "intermediate_distinct", C.MULTI_STEP)


ROW_CASES = [  # (row, id, case, state, rules, expected decision)
    (1, "elementary", CASE, replace(S, claim=CLAIM, connection="elementary"), RULES, ELEMENTARY),
    (1, "degenerate", DEGENERATE, replace(S, claim=CLAIM, connection="degenerate"), RULES,
     Decision(A.COMPLETE, "connection:degenerate", C.DEGENERATE)),
    (1, "reassigned", CASE, replace(S, claim=CLAIM, connection="reassigned"), RULES,
     Decision(A.COMPLETE, "connection:reassigned", C.REASSIGNED)),
    (2, "distinct", CASE, replace(S, intermediate="distinct"), RULES, DISTINCT),
    (3, "screen", CASE, replace(S, screen=BV(verdict="barrierless", source="screen"),
                                seeds=(SEED,)), RULES,
     Decision(A.COMPLETE, "screen:barrierless", C.BARRIERLESS)),
    (3, "string", CASE, replace(S, screen=BV(verdict="barrierless", source="string"),
                                path_runs=1), RULES,
     Decision(A.COMPLETE, "string:barrierless", C.BARRIERLESS)),
    (3, "association_scan", ASSOCIATION, replace(S, screen=BV(verdict="barrierless",
                                                              source="scan")), RULES,
     Decision(A.COMPLETE, "scan:barrierless", C.BARRIERLESS)),
    (4, "converged", CASE, replace(S, saddle_pending=True, attempts=1), RULES, VALIDATE),
    (4, "given_ts_calc", CASE.model_copy(update={"ts_calc": "calc_ts"}),
     replace(BASE, saddle_pending=True), RULES, VALIDATE),
    (5, "screen", CASE, BASE, RULES, Decision(A.SCREEN, "screen")),
    (5, "degenerate_goes_on", DEGENERATE, BASE, RULES, Decision(A.SCREEN, "screen")),
    (5, "after_the_shortcut_seeds", CASE, SHORTCUT, RULES, Decision(A.SCREEN, "screen")),
    (5, "association_without_a_low_level_engine", ASSOCIATION, BASE,
     replace(RULES, screen=False), Decision(A.SCREEN, "screen")),
    (6, "intermediate", CASE, replace(S, screen=BV(verdict="intermediate", source="screen"),
                                      seeds=(SEED,)), RULES,
     Decision(A.VALIDATE_INTERMEDIATE, "path_intermediate")),
    (7, "seed", CASE, replace(S, seeds=(SEED,)), RULES, Decision(A.REFINE_SADDLE, "seed:screen_ts")),
    (7, "continuation", CASE, replace(S, seeds=(CONTINUED,), attempts=1), RULES,
     Decision(A.REFINE_SADDLE, "seed:continuation")),
    (7, "the_next_shortcut_seed", CASE, replace(SHORTCUT, seeds=(OTHER_TS,)), RULES,
     Decision(A.REFINE_SADDLE, "seed:discovery_ts")),
    (8, "no_path", CASE, replace(S, attempts=1), RULES, DFT_PATH),
    (8, "screen_off", CASE, BASE, replace(RULES, screen=False), DFT_PATH),
    (8, "unavailable", CASE, replace(S, screen=BV(verdict="unavailable", source="screen")),
     RULES, DFT_PATH),
    (8, "after_the_screen_path", CASE, replace(SHORTCUT, neb_done=True), RULES, DFT_PATH),
    (8, "its_seed_failed", CASE, replace(S, screen=STRING, path_runs=1, attempts=1), RULES,
     DFT_PATH),
    # U6-P4: an unavailable string is a token itself, so the next chunk follows it
    (8, "after_an_unavailable_string", CASE, replace(S, screen=BV(
        verdict="unavailable", source="string"), path_runs=1, attempts=1), RULES, DFT_PATH),
    (8, "after_an_endpoint_well", CASE, replace(S, screen=BV(verdict="intermediate",
                                                             source="string"),
                                                intermediate="same_as_endpoint", path_runs=1,
                                                attempts=1), RULES, DFT_PATH),
    (9, "exhausted", CASE, replace(S, path_runs=1, seeds=(SEED,), attempts=2), RULES, EXHAUSTED),
    (9, "a_continuation_at_the_limit", CASE, replace(S, seeds=(CONTINUED,), attempts=2), RULES,
     EXHAUSTED),
    # G3-P2: the second string left no seed that was tried
    (9, "a_string_without_a_seed", CASE, replace(S, screen=STRING, path_runs=2, attempts=1),
     RULES, EXHAUSTED),
    (9, "no_string_without_attempts", CASE, replace(S, attempts=2), RULES, EXHAUSTED),
    (9, "no_chunk_without_attempts", CASE, replace(S, screen=STRING, path_runs=1, attempts=2),
     RULES, EXHAUSTED),
    # X4: no string joins an association's monomers, after a failed scan or its failed seed
    (9, "association_scan_point_failed", ASSOCIATION, replace(S, screen=BV(
        verdict="unavailable", source="scan", reasons=("scan_point",))), RULES, EXHAUSTED),
    (9, "association_seed_failed", ASSOCIATION, replace(S, screen=BV(verdict="single",
                                                                     source="scan"),
                                                        attempts=1), RULES, EXHAUSTED),
]


def test_table_has_9_rows_and_every_row_is_covered():
    assert len(st.ROWS) == 9
    assert {row for row, *_ in ROW_CASES} == set(range(1, 10))


@pytest.mark.parametrize(("row", "case", "state", "rules", "expected"),
                         [c[:1] + c[2:] for c in ROW_CASES],
                         ids=[f"row{c[0]}-{c[1]}" for c in ROW_CASES])
def test_row(row, case, state, rules, expected):
    assert decide(case, state, rules) == expected


def test_a_well_is_validated_before_any_seed_then_the_seeds_go_on():
    """eng P0-4: an intermediate profile adds no peak seed and validates its well first."""
    well = BV(verdict="intermediate", source="string")
    state = record_profile(replace(S, seeds=(SEED,)), well, (replace(SEED, source="path_hei"),))
    assert state.seeds == (SEED,) and state.path_runs == 1 and state.screen == well
    assert decide(CASE, state, RULES).action is A.VALIDATE_INTERMEDIATE
    after = replace(state, intermediate="same_as_endpoint")
    assert decide(CASE, after, RULES) == Decision(A.REFINE_SADDLE, "seed:screen_ts")


def test_a_profile_seeds_only_a_single_step_and_keeps_its_seeds_in_order():
    hei = replace(SEED, source="path_hei")
    assert record_profile(S, STRING, (hei,)).seeds == (hei,)
    assert record_profile(S, BV(verdict="single", source="screen"),
                          (SEED, OTHER_TS)).seeds == (SEED, OTHER_TS)
    for verdict in ("barrierless", "intermediate", "unavailable"):
        assert record_profile(S, BV(verdict=verdict, source="string"), (hei,)).seeds == ()


def test_a_new_profile_reopens_the_checks_of_an_earlier_well():
    """A string after an earlier well: its own well is validated."""
    old = replace(S, attempts=2, intermediate="same_as_endpoint")
    state = record_profile(old, BV(verdict="intermediate", source="string"))
    assert state.intermediate is None
    assert decide(CASE, state, RULES) == Decision(A.VALIDATE_INTERMEDIATE, "path_intermediate")
    assert record_profile(BASE, BV(verdict="single", source="screen")).path_runs == 0


def test_every_failure_is_one_token_and_the_search_goes_on_while_the_budget_lasts():
    """U6-P4: a rejected connection, a rejected TS or an unconverged search leaves the state
    without claim or connection, one token spent: the next seed, else the next path, while
    max_saddle_attempts lasts (R2 never reached the old connection_failed row)."""
    failed = replace(S, attempts=1)  # what any failed attempt leaves
    assert decide(CASE, failed, RULES) == DFT_PATH
    assert decide(CASE, replace(failed, seeds=(SEED,)), RULES) == Decision(
        A.REFINE_SADDLE, "seed:screen_ts")
    assert decide(CASE, replace(failed, path_runs=1), RULES) == DFT_PATH  # its next chunk
    assert decide(CASE, replace(failed, path_runs=2), RULES) == EXHAUSTED
    assert decide(CASE, replace(failed, attempts=2), RULES) == EXHAUSTED
    larger = replace(RULES, budget=st.ReactionPathsPolicy(max_saddle_attempts=3))
    assert decide(CASE, replace(failed, attempts=2, path_runs=1), larger) == DFT_PATH


def test_the_screen_path_runs_once_after_the_failed_shortcut_seeds():
    """R4 / G8-P7: the shortcut's seeds (the low-level TSs whose DFT SP lies at least one
    resolution above both minima) go first, in order, and give no verdict; once they have
    failed, the low-level path runs before any string, and never again."""
    seeded = replace(SHORTCUT, seeds=(SEED, OTHER_TS), attempts=0)
    assert decide(CASE, seeded, RULES) == Decision(A.REFINE_SADDLE, "seed:screen_ts")
    assert decide(CASE, SHORTCUT, RULES) == Decision(A.SCREEN, "screen")
    assert decide(CASE, replace(SHORTCUT, neb_done=True), RULES) == DFT_PATH
    assert decide(CASE, SHORTCUT, replace(RULES, screen=False)) == DFT_PATH


def test_strings_run_only_after_the_seed_of_each_earlier_string():
    """G3-P2: one FIND_PATH row. The first string, then a next chunk once the latest string's
    seed has failed, while attempts are left; a string whose verdict left no seed to try (an
    intermediate whose well relaxed to an end and rose no peak) ends the case."""
    first = replace(BASE, neb_done=True)
    assert decide(CASE, first, replace(RULES, screen=False)) == DFT_PATH
    unseeded = replace(first, screen=BV(verdict="intermediate", source="string"), path_runs=1,
                       intermediate="relax_failed")
    assert decide(CASE, unseeded, RULES) == EXHAUSTED
    tried = replace(unseeded, attempts=1)  # SCREEN's or the string's seed was tried
    assert decide(CASE, tried, RULES) == DFT_PATH
    assert decide(CASE, replace(tried, path_runs=2), RULES) == EXHAUSTED

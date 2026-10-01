"""decide(): the 14 rows of design §7.3 (ports tests/test_reaction_classification.py cases)."""

from __future__ import annotations

from dataclasses import replace

import pytest

from hfauto.core import records as r
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
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
HIGHER = st.Seed(geometry=GEO, source="higher_order_retry", depth=1)
CLAIM = r.SaddleClaim(saddle_calc="s", freq_calc="f", imag_cm1=-1131.0, energy_hartree=-93.0)
SOFT = CLAIM.model_copy(update={"imag_cm1": -30.0})  # |nu| < saddle_cm1 = 50
RULES = st.CaseRules()


def minimum(mid: str, basin: str, kcal: float = 0.0) -> r.MinimumRecord:
    return r.MinimumRecord(minimum_id=mid, basin_id=basin, composition_id="CHN_q0_m1",
                           species_id="s" + mid, tier="dft", level_key="L", opt_calc="o",
                           freq_calc="f", energy_hartree=-93.0 + kcal / HARTREE_TO_KCAL_MOL,
                           state_label="x")


MA, MB = minimum("ma", "A"), minimum("mb", "B", 10.0)
BASE = CaseState(minima=(MA, MB))  # nothing done yet
S = replace(BASE, screen=BV(verdict="single", source="screen"), neb_done=True)  # one peak
SHORTCUT = replace(BASE, screen=BV(verdict="single", source="screen"), saddle_attempts=1,
                   last_saddle="failed")  # the seed of a low-level TS has failed
STRING = BV(verdict="single", source="string")
BLOCKED = Decision(A.BLOCKED, "endpoint_without_dft_minimum", C.BLOCKED)
DFT_PATH = Decision(A.FIND_PATH, "dft_path")
VALIDATE = Decision(A.VALIDATE_AND_CONNECT, "saddle_converged")
EXHAUSTED = Decision(A.COMPLETE, "attempts_exhausted", C.UNRESOLVED)
ELEMENTARY = Decision(A.COMPLETE, "connection:elementary", C.ELEMENTARY_STEP)
DISTINCT = Decision(A.COMPLETE, "intermediate_distinct", C.MULTI_STEP)
WALLTIME = Decision(A.COMPLETE, "walltime", C.UNRESOLVED)
FAILED = Decision(A.COMPLETE, "connection_failed", C.UNRESOLVED)


ROW_CASES = [  # (row, id, case, state, rules, expected decision)
    (1, "reactant_missing", CASE, replace(BASE, minima=(None, MB)), RULES, BLOCKED),
    (1, "product_missing", CASE, replace(S, minima=(MA, None), path_runs=2), RULES, BLOCKED),
    (2, "same_basin", CASE, replace(BASE, minima=(MA, minimum("mb", "A")), claim=CLAIM), RULES,
     Decision(A.COMPLETE, "same_basin", C.SAME_BASIN)),
    # X4: a complex that relaxed into its adduct leaves one basin; the monomers are the reactant
    (10, "association_from_the_adduct_basin", ASSOCIATION, replace(BASE, minima=(MA, MA)), RULES,
     Decision(A.SCREEN, "screen")),
    (3, "uphill", CASE, replace(BASE, minima=(MA, minimum("mb", "B", 41.0))), RULES,
     Decision(A.COMPLETE, "out_of_window", C.OUT_OF_WINDOW)),
    (4, "elementary", CASE, replace(S, claim=CLAIM, connection="elementary"), RULES, ELEMENTARY),
    (4, "at_the_deadline", CASE, replace(S, expired=True, claim=CLAIM, connection="elementary"),
     RULES, ELEMENTARY),  # U9-P5: evidence in hand completes the case
    (4, "degenerate", DEGENERATE, replace(S, minima=(MA, MA), claim=CLAIM,
                                          connection="degenerate"), RULES,
     Decision(A.COMPLETE, "connection:degenerate", C.DEGENERATE)),
    (4, "reassigned", CASE, replace(S, claim=CLAIM, connection="reassigned"), RULES,
     Decision(A.COMPLETE, "connection:reassigned", C.REASSIGNED)),
    (5, "distinct", CASE, replace(S, intermediate="distinct"), RULES, DISTINCT),
    (5, "at_the_deadline", CASE, replace(S, expired=True, intermediate="distinct"), RULES,
     DISTINCT),
    (6, "screen", CASE, replace(S, screen=BV(verdict="barrierless", source="screen"),
                                seeds=(SEED,)), RULES,
     Decision(A.COMPLETE, "screen:barrierless", C.BARRIERLESS)),
    (6, "string_at_the_deadline", CASE, replace(S, screen=BV(verdict="barrierless",
                                                             source="string"),
                                                path_runs=1, expired=True), RULES,
     Decision(A.COMPLETE, "string:barrierless", C.BARRIERLESS)),
    (6, "association_scan", ASSOCIATION, replace(S, screen=BV(verdict="barrierless",
                                                              source="scan")), RULES,
     Decision(A.COMPLETE, "scan:barrierless", C.BARRIERLESS)),
    (7, "walltime", CASE, replace(S, expired=True, last_saddle="converged"), RULES, WALLTIME),
    (7, "a_rejected_connection", CASE, replace(S, expired=True, claim=CLAIM,
                                               connection="same_basin"), RULES, WALLTIME),
    (8, "failed", CASE, replace(S, claim=CLAIM, connection="failed"), RULES, FAILED),
    (8, "same_basin_without_attempts", CASE, replace(S, claim=CLAIM, connection="same_basin",
                                                     saddle_attempts=2), RULES, FAILED),
    (8, "same_state_without_attempts", CASE, replace(S, claim=CLAIM, connection="same_state",
                                                     saddle_attempts=2), RULES, FAILED),
    (9, "converged", CASE, replace(S, last_saddle="converged", saddle_attempts=1), RULES,
     VALIDATE),
    (9, "given_ts_calc", CASE.model_copy(update={"ts_calc": "calc_ts"}),
     replace(BASE, last_saddle="converged"), RULES, VALIDATE),
    (10, "screen", CASE, BASE, RULES, Decision(A.SCREEN, "screen")),
    (10, "degenerate_goes_on", DEGENERATE, replace(BASE, minima=(MA, MA)), RULES,
     Decision(A.SCREEN, "screen")),
    (10, "after_the_shortcut_seeds", CASE, SHORTCUT, RULES, Decision(A.SCREEN, "screen")),
    (10, "a_given_ts_failed", CASE, replace(BASE, last_saddle="failed"), RULES,
     Decision(A.SCREEN, "screen")),
    (10, "association_without_a_low_level_engine", ASSOCIATION, BASE,
     replace(RULES, screen=False), Decision(A.SCREEN, "screen")),
    (11, "intermediate", CASE, replace(S, screen=BV(verdict="intermediate", source="screen"),
                                       seeds=(SEED,)), RULES,
     Decision(A.VALIDATE_INTERMEDIATE, "path_intermediate")),
    (12, "seed", CASE, replace(S, seeds=(SEED,)), RULES,
     Decision(A.REFINE_SADDLE, "seed:screen_ts")),
    (12, "higher_order_retry", CASE, replace(S, seeds=(HIGHER,), last_saddle="failed",
                                             saddle_attempts=1), RULES,
     Decision(A.REFINE_SADDLE, "seed:higher_order_retry")),
    (12, "the_next_shortcut_seed", CASE, replace(SHORTCUT, seeds=(OTHER_TS,)), RULES,
     Decision(A.REFINE_SADDLE, "seed:discovery_ts")),
    (13, "no_path", CASE, replace(S, last_saddle="failed", saddle_attempts=1), RULES, DFT_PATH),
    (13, "screen_off", CASE, BASE, replace(RULES, screen=False), DFT_PATH),
    (13, "unavailable", CASE, replace(S, screen=BV(verdict="unavailable", source="screen")),
     RULES, DFT_PATH),
    (13, "after_the_screen_path", CASE, replace(SHORTCUT, neb_done=True), RULES, DFT_PATH),
    (13, "its_seed_failed", CASE, replace(S, screen=STRING, path_runs=1, saddle_attempts=1,
                                          last_saddle="failed"), RULES, DFT_PATH),
    (13, "after_an_endpoint_well", CASE, replace(S, screen=BV(verdict="intermediate",
                                                              source="string"),
                                                 intermediate="same_as_endpoint", path_runs=1,
                                                 saddle_attempts=1), RULES, DFT_PATH),
    (14, "exhausted", CASE, replace(S, path_runs=1, seeds=(SEED,), saddle_attempts=2,
                                    last_saddle="failed"), RULES, EXHAUSTED),
    (14, "a_push_at_the_limit", CASE, replace(S, seeds=(HIGHER,), saddle_attempts=2), RULES,
     EXHAUSTED),
    # G3-P2: the second string left no seed that was tried
    (14, "a_string_without_a_seed", CASE, replace(S, screen=STRING, path_runs=2,
                                                  saddle_attempts=1), RULES, EXHAUSTED),
    (14, "string_failed", CASE, replace(S, path_runs=1, saddle_attempts=1), RULES, EXHAUSTED),
    (14, "unavailable_string", CASE, replace(S, screen=BV(verdict="unavailable", source="string"),
                                             path_runs=1, saddle_attempts=1), RULES, EXHAUSTED),
    # R3: no string whose peak seed could never be refined
    (14, "no_string_without_attempts", CASE, replace(S, saddle_attempts=2, last_saddle="failed"),
     RULES, EXHAUSTED),
    (14, "no_chunk_without_attempts", CASE, replace(S, screen=STRING, path_runs=1,
                                                    saddle_attempts=2), RULES, EXHAUSTED),
    # X4: no string joins an association's monomers, after a failed scan or its failed seed
    (14, "association_scan_point_failed", ASSOCIATION, replace(S, screen=BV(
        verdict="unavailable", source="scan", reasons=("scan_point",))), RULES, EXHAUSTED),
    (14, "association_seed_failed", ASSOCIATION, replace(S, screen=BV(verdict="single",
                                                                      source="scan"),
                                                         saddle_attempts=1, last_saddle="failed"),
     RULES, EXHAUSTED),
]


def test_table_has_14_rows_and_every_row_is_covered():
    assert len(st.ROWS) == 14
    assert {row for row, *_ in ROW_CASES} == set(range(1, 15))


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


def test_a_new_profile_reopens_the_checks_of_an_earlier_saddle():
    """A string after a validated saddle and an earlier well: its own well is validated."""
    old = replace(S, saddle_attempts=2, claim=CLAIM, intermediate="same_as_endpoint")
    state = record_profile(old, BV(verdict="intermediate", source="string"))
    assert (state.last_saddle, state.intermediate) == (None, None)
    assert decide(CASE, state, RULES) == Decision(A.VALIDATE_INTERMEDIATE, "path_intermediate")
    assert record_profile(BASE, BV(verdict="single", source="screen")).path_runs == 0


def test_seeds_share_one_attempt_budget():
    """U6-P3: a higher-order push is an ordinary counted seed; a used budget ends the case
    without a string (R3). A stalled search's restart runs inside its attempt (refine_saddle)."""
    first = replace(S, seeds=(HIGHER,), last_saddle="failed", saddle_attempts=1, path_runs=1)
    assert decide(CASE, first, RULES) == Decision(A.REFINE_SADDLE, "seed:higher_order_retry")
    assert decide(CASE, replace(first, saddle_attempts=2), RULES) == EXHAUSTED


def test_the_screen_path_runs_once_after_the_failed_shortcut_seeds():
    """R4 / G8-P7: the shortcut's seeds (every low-level TS of the pair that gives a single
    step) go first, in order; once they have failed, the low-level path runs before any string,
    and never again."""
    seeded = replace(SHORTCUT, seeds=(SEED, OTHER_TS), saddle_attempts=0, last_saddle=None)
    assert decide(CASE, seeded, RULES) == Decision(A.REFINE_SADDLE, "seed:screen_ts")
    assert decide(CASE, SHORTCUT, RULES) == Decision(A.SCREEN, "screen")
    assert decide(CASE, replace(SHORTCUT, neb_done=True), RULES) == DFT_PATH
    screen_off = replace(RULES, screen=False)
    assert decide(CASE, SHORTCUT, screen_off) == DFT_PATH


def test_strings_run_only_after_the_seed_of_each_earlier_string():
    """G3-P2: one FIND_PATH row. The first string, then a next chunk once the latest string's
    seed has failed, while attempts are left; a string whose verdict left no seed to try (an
    intermediate whose well relaxed to an end and rose no peak) ends the case instead of running
    up to three chunks as before."""
    first = replace(BASE, neb_done=True)
    assert decide(CASE, first, replace(RULES, screen=False)) == DFT_PATH
    unseeded = replace(first, screen=BV(verdict="intermediate", source="string"), path_runs=1,
                       intermediate="relax_failed")
    assert decide(CASE, unseeded, RULES) == EXHAUSTED
    tried = replace(unseeded, saddle_attempts=1)  # SCREEN's or the string's seed was tried
    assert decide(CASE, tried, RULES) == DFT_PATH
    assert decide(CASE, replace(tried, path_runs=2), RULES) == EXHAUSTED


@pytest.mark.parametrize("claim", [CLAIM, SOFT])
def test_a_ts_whose_sides_join_one_basin_leaves_the_search_open(claim):
    """VAL7/s6_oh_ch4: a -73i reorientation saddle of the OH...CH4 complex (both QRC sides in
    its basin at both amplitudes) is not the abstraction TS; with a saddle attempt left the case
    goes on to the string (the W3 run found the -481i TS there) instead of stopping. G3-P1: a
    soft TS (|nu| < saddle_cm1) is a TS like any other; its failed QRC splits no well."""
    other = replace(S, claim=claim, saddle_attempts=1, connection="same_basin")
    assert decide(CASE, other, RULES) == DFT_PATH
    seeded = replace(other, seeds=(SEED,))
    assert decide(CASE, seeded, RULES) == Decision(A.REFINE_SADDLE, "seed:screen_ts")
    assert decide(CASE, replace(other, path_runs=1), RULES) == EXHAUSTED
    assert decide(CASE, replace(other, saddle_attempts=2), RULES) == FAILED
    failed = replace(other, connection="failed")  # a side did not optimize: no new evidence
    assert decide(CASE, failed, RULES) == FAILED


@pytest.mark.parametrize("claim", [CLAIM, SOFT])
def test_a_ts_whose_sides_join_one_state_is_not_retried_and_the_search_goes_on(claim):
    """X2-2: both QRC sides in one case key but two basins (one state of a bond-changing case):
    a wider displacement cannot help; with a saddle attempt left the search goes on, the next
    saddle search dropping the claim."""
    other = replace(S, claim=claim, saddle_attempts=1, connection="same_state")
    assert decide(CASE, other, RULES) == DFT_PATH
    seeded = replace(other, seeds=(SEED,))
    assert decide(CASE, seeded, RULES) == Decision(A.REFINE_SADDLE, "seed:screen_ts")

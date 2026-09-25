"""Structures for minima(dft) (§8.2); K cases ported from test_dft_minima_refinement."""

from hfauto.chemistry.selection import Candidate, rerank, select_for_refinement

KCAL = 1 / 627.5094740631


def ids(candidates):
    return [c.species_id for c in candidates]


def test_window_is_bounded_per_state_and_balanced_across_formulas():
    pool = [Candidate(f"nci_{i}", "H2_q0_m1", "s", -100.0 + i * KCAL) for i in range(3)]
    pool += [Candidate("h2_far", "H2_q0_m1", "s", -100.0 + 7 * KCAL),  # outside 6 kcal/mol
             Candidate("h2_bent", "H2_q0_m1", "t", -90.0),  # another state label: own budget
             Candidate("h3", "H3_q0_m1", "s", -1.0)]  # never compared with H2 energies
    assert ids(select_for_refinement(pool, per_state=2)) == ["nci_0", "nci_1", "h2_bent", "h3"]
    assert ids(select_for_refinement(pool)) == ["nci_0", "nci_1", "nci_2", "h2_bent", "h3"]


def test_discovery_endpoints_always_kept_unranked_dropped_and_rerank():
    pool = [Candidate("low", "A", "s", -1.0), Candidate("mid", "A", "s", -1.0 + 1 * KCAL),
            Candidate("source", "A", "s", -1.0 + 9 * KCAL, always=True),
            Candidate("product", "A", "p", None, always=True),
            Candidate("unscreened", "A", "u", None)]
    assert ids(select_for_refinement(pool, per_state=1)) == ["low", "source", "product"]
    sp = {"low": -2.0, "mid": -2.3}  # re-ranked by single points; no point, no refinement
    assert ids(rerank(pool, sp, keep_per_state=1)) == ["mid", "source", "product"]

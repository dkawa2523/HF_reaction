"""Structures for minima(dft) (§8.2): the DFT entry of the low-level edges (U4-P3) and the
refined structures; K cases ported from test_dft_minima_refinement."""

import pytest

from hfauto.chemistry.selection import (
    Candidate,
    Edge,
    admit,
    reacting_candidates,
    rerank,
    select_for_refinement,
)
from hfauto.core.evidence import FileRef, Geometry
from hfauto.core.records import DiscoveryRecord, MinimumRecord, SpeciesRecord
from hfauto.core.system import SystemConfig

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


def scorer(energies, asked):
    """Single points from ``energies``; ``asked`` records which candidates were scored."""
    def single_points(crowd):
        asked.extend(ids(crowd))
        return {sid: energies[sid] for sid in ids(crowd) if sid in energies}
    return single_points


def test_always_kept_candidates_unranked_dropped_and_rerank():
    pool = [Candidate("low", "A", "s", -1.0), Candidate("mid", "A", "s", -1.0 + 1 * KCAL),
            Candidate("source", "A", "s", -1.0 + 9 * KCAL, always=True),
            Candidate("product", "A", "p", None, always=True),
            Candidate("unscreened", "A", "u", None)]
    assert ids(select_for_refinement(pool, per_state=1)) == ["low", "source", "product"]
    chosen, asked = select_for_refinement(pool, per_state=8), []
    sp = {"low": -2.0, "mid": -2.3}  # re-ranked by single points
    assert ids(rerank(chosen, scorer(sp, asked), 1, window_kcal=6.0)) == [
        "mid", "source", "product"]
    assert asked == ["low", "mid"]  # always kept: never a single point


def test_rerank_cuts_crowded_groups_by_single_points_and_keeps_small_groups_whole():
    pool = [Candidate(f"c{i}", "A", "s", -1.0 + i * KCAL) for i in range(4)]
    pool += [Candidate("source", "A", "s", -1.0, always=True),
             Candidate("b0", "B", "s", -3.0), Candidate("b1", "B", "s", -3.0 + 5 * KCAL)]
    sp, asked = {"c0": -2.0, "c1": -2.0 + 5 * KCAL, "c2": -2.0 + 7 * KCAL}, []  # c3 failed
    assert ids(rerank(pool, scorer(sp, asked), 3, window_kcal=6.0)) == [
        "c0", "c1", "source", "b0", "b1"]  # c2 7 kcal/mol up
    assert asked == ["c0", "c1", "c2", "c3"]  # B holds only 2: no single point
    asked.clear()
    assert ids(rerank(pool, scorer({}, asked), 4, window_kcal=6.0)) == ids(pool)
    assert asked == []  # no crowded group


def spc(sid, composition, symbols, energy=None, state="s"):
    geometry = Geometry(file=FileRef(path=f"in/{sid}.xyz", sha256="0"), fingerprint=sid,
                        symbols=symbols)
    return SpeciesRecord(species_id=sid, composition_id=composition, charge=0, multiplicity=1,
                         geometry=geometry, source="input", state_label=state,
                         energy_hartree=energy)


def screen_minimum(sid, composition, energy, members=()):
    return MinimumRecord(minimum_id=f"min_{sid}", basin_id=f"basin_{sid}",
                         composition_id=composition, species_id=sid, tier="screen",
                         level_key="gfn2", opt_calc=f"opt_{sid}", freq_calc=f"freq_{sid}",
                         energy_hartree=energy, state_label="s", members=(sid, *members))


ASKED = []


def edge(did, source, product, *points, composition="A", first=True):
    return Edge(did, composition, source, product,
                lambda: ASKED.append(did) or tuple(points), first)


def test_an_edge_is_admitted_by_its_height_above_the_start_state():
    """U4-P3: the window is judged from the start state (explore's first generation), over the
    highest point on the way: a second-generation TS 15 kcal/mol above its own source but 35
    above the start is out of a 32 kcal/mol window; a degenerate step (dE_rxn 0) is judged by
    its barrier alone."""
    k = KCAL
    edges = [edge("e1", "s0", "p1", 0.0, 20 * k, 30 * k),
             edge("e2", "p1", "p2", 20 * k, 10 * k, 35 * k, first=False),
             edge("swap", "s0", "s0_swapped", 0.0, 0.0, 33 * k)]
    verdicts = admit(edges, window_kcal=32.0, per_composition=6)
    assert verdicts["e1"] == (None, pytest.approx(30.0))
    assert verdicts["e2"] == ("out_of_window", pytest.approx(35.0))
    assert verdicts["swap"] == ("out_of_window", pytest.approx(33.0))
    assert admit(edges, window_kcal=35.5, per_composition=6)["e2"][0] is None


def test_single_points_are_taken_only_for_edges_whose_source_is_within_the_window():
    """W4 S10: the height only grows along a way, so an edge leaving a state reached above the
    window is above it too; its single points are never taken (850 edges after the budget)."""
    k = KCAL
    ASKED.clear()
    edges = [edge("up", "s0", "high", 0.0, 30 * k, 45 * k),
             edge("past", "high", "far", 30 * k, 10 * k, 31 * k, first=False),
             edge("farther", "far", "end", 10 * k, 0.0, 12 * k, first=False),
             edge("near", "s0", "p1", 0.0, 2 * k, 5 * k)]
    verdicts = admit(edges, window_kcal=40.0, per_composition=6)
    assert verdicts == {"up": ("out_of_window", pytest.approx(45.0)),
                        "past": ("out_of_window", None), "farther": ("out_of_window", None),
                        "near": (None, pytest.approx(5.0))}
    assert sorted(ASKED) == ["near", "up"]


def test_a_relaxed_seed_is_a_start_and_the_lowest_start_height_counts():
    """S5: the seed (CH3 + O2, a start) relaxes without a TS into the well (CH3O2), itself a
    start of its own edges; a TS 36 kcal/mol above the well lies 6 above the seed."""
    k = KCAL
    edges = [edge("relax", "seed", "well", 30 * k, 0.0),
             edge("h_shift", "well", "ch2ooh", 0.0, 20 * k, 36 * k)]
    verdicts = admit(edges, window_kcal=10.0, per_composition=6)
    assert verdicts == {"relax": (None, 0.0), "h_shift": (None, pytest.approx(6.0))}
    alone = admit(edges[1:], window_kcal=10.0, per_composition=6)  # the well its only start
    assert alone == {"h_shift": ("out_of_window", pytest.approx(36.0))}


def test_a_state_is_one_node_at_its_lowest_conformer():
    """W4 fluoroethane: a second-generation edge leaves its state from another conformer than
    the one the first edge reached; it is reached all the same, and a start state is measured
    from its lowest conformer (each edge's end single points are at its own structures)."""
    k = KCAL
    edges = [edge("e3", "s0", "p3", 0.0, 1 * k, 4 * k),
             edge("e1", "s0", "p1", 2 * k, 10 * k, 30 * k),  # a conformer 2 above the lowest
             edge("e2", "p1", "p2", 12 * k, 5 * k, 33 * k, first=False)]
    verdicts = admit(edges, window_kcal=40.0, per_composition=6)
    assert verdicts == {"e1": (None, pytest.approx(30.0)), "e2": (None, pytest.approx(33.0)),
                        "e3": (None, pytest.approx(4.0))}


def test_below_a_high_start_edges_are_ordered_by_their_energetic_span():
    """W4 S5: the unrelaxed seed lies above everything it leads to, so every height is 0; the
    order (and with it the count cap) then follows the energetic span on the way: the downhill
    association first, then the climb out of the well, then a downhill step reached only over a
    higher climb (not the discovery-id order that cut the association in the first S5 run)."""
    k = KCAL
    edges = [edge("a_far", "well", "far", -60 * k, -10 * k, -10 * k),  # TS 50 above the well
             edge("b_down", "far", "deep", -10 * k, -70 * k, -9 * k, first=False),
             edge("c_shift", "well", "ooh", -60 * k, -40 * k, -25 * k),  # TS 35 above the well
             edge("z_relax", "seed", "well", 0.0, -60 * k)]
    verdicts = admit(edges, window_kcal=40.0, per_composition=2)
    assert {d: v for d, v in verdicts.items()} == {
        "z_relax": (None, 0.0), "c_shift": (None, 0.0),
        "a_far": ("over_cap", 0.0), "b_down": ("over_cap", 0.0)}
    # equal spans (every step past the climb out of the well): the step itself comes first
    after = [edge("a_on", "ooh", "x", -40 * k, -50 * k, -45 * k, first=False), *edges]
    assert admit(after, window_kcal=40.0, per_composition=2)["c_shift"] == (None, 0.0)


def test_a_failed_single_point_leaves_its_edge_and_what_it_leads_to_not_evaluated():
    k = KCAL
    edges = [edge("e1", "s0", "p1", 0.0, 5 * k, "scf_not_converged"),
             edge("e2", "p1", "p2", 5 * k, 1 * k, 9 * k, first=False),
             edge("e3", "s0", "p3", 0.0, 2 * k, 4 * k)]
    assert admit(edges, window_kcal=40.0, per_composition=6) == {
        "e1": ("not_evaluated:scf_not_converged", None),
        "e2": ("not_evaluated:unreached", None), "e3": (None, pytest.approx(4.0))}


def test_the_count_cap_takes_pairs_of_ends_per_composition_lowest_first():
    k = KCAL
    edges = [edge("high", "s", "p1", 0.0, 0.0, 20 * k), edge("low", "s", "p2", 0.0, 0.0, 10 * k),
             edge("low_again", "s", "p2", 0.0, 0.0, 15 * k),  # another TS of an admitted pair
             edge("other", "t", "q", 0.0, 0.0, 30 * k, composition="B")]
    verdicts = admit(edges, window_kcal=40.0, per_composition=1)
    assert {d: v[0] for d, v in verdicts.items()} == {
        "high": "over_cap", "low": None, "low_again": None, "other": None}


def found(discovery_id, source, product):
    return DiscoveryRecord(discovery_id=discovery_id, mechanism="nt2", outcome="product",
                           source_species=source, product_species=product)


def system(*endpoints):
    declared = [{"id": "hf", "xyz": "hf.xyz", "multiplicity": 1},
                {"id": "h2", "xyz": "h2.xyz", "multiplicity": 1},
                *({"id": e, "xyz": f"{e}.xyz", "multiplicity": 1, "role": "endpoint"}
                  for e in endpoints)]
    return SystemConfig(system_id="t", species=declared,
                        compositions=[{"id": "hf2", "components": {"hf": 2}},
                                      {"id": "hf_h2", "components": {"hf": 1, "h2": 1}}])


def test_reacting_candidates_are_the_screen_minima_and_admitted_ends_of_reacting_states():
    """(HF)2 reacts through its admitted edges, HF as its monomer; H2 and the HF·H2 complex do
    not, unless the complex is a declared endpoint (its monomers HF and H2 then react too). An
    end a screen basin holds stands for that basin; one none holds (a new product, a seed
    copy) for itself."""
    dimer, mixed = ("H", "F", "H", "F"), ("H", "F", "H", "H")
    species = {s.species_id: s for s in (
        spc("hf", "FH_q0_m1", ("H", "F")), spc("h2", "H2_q0_m1", ("H", "H")),
        spc("d1", "F2H2_q0_m1", dimer), spc("d2", "F2H2_q0_m1", dimer),
        spc("d3", "F2H2_q0_m1", dimer), spc("x", "FH3_q0_m1", mixed),
        spc("p", "F2H2_q0_m1", dimer, energy=-199.0, state="p"),
        spc("seed", "F2H2_q0_m1", dimer, state="q"))}
    screen = [screen_minimum("hf", "FH_q0_m1", -100.0),
              screen_minimum("d1", "F2H2_q0_m1", -200.0, members=("d2",)),
              screen_minimum("d3", "F2H2_q0_m1", -199.9), screen_minimum("x", "FH3_q0_m1", -50.0),
              screen_minimum("h2", "H2_q0_m1", -1.0)]
    admitted = [found("to_member", "d3", "d2"),  # d2 is held by d1's basin
                found("to_new", "d3", "p"),  # p skipped screen: it stands for itself
                found("relax", "seed", "d1")]  # a seed copy no screen basin holds

    pool, idle = reacting_candidates(species, admitted, screen, system())
    assert ids(pool) == ["hf", "d1", "d3", "p", "seed"]
    assert [c.always for c in pool] == [False, True, True, True, True]
    assert [c.energy_hartree for c in pool] == [-100.0, -200.0, -199.9, -199.0, None]
    assert idle == ["FH3_q0_m1", "H2_q0_m1"]

    pool, idle = reacting_candidates(species, admitted, screen, system("x"))
    assert ids(pool) == ["hf", "d1", "d3", "x", "h2", "p", "seed"] and idle == []
    assert [c.always for c in pool][3:5] == [True, False]  # the endpoint, not its monomer

    pool, idle = reacting_candidates(species, [], screen, system())
    assert pool == [] and idle == ["F2H2_q0_m1", "FH3_q0_m1", "FH_q0_m1", "H2_q0_m1"]

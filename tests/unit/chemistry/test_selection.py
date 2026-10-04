"""Structures for minima(dft) (§8.2); K cases ported from test_dft_minima_refinement."""

from hfauto.chemistry.selection import Candidate, reacting_candidates, rerank, select_for_refinement
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


def test_discovery_endpoints_always_kept_unranked_dropped_and_rerank():
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


def found(discovery_id, source_minimum, product=None, outcome="product"):
    return DiscoveryRecord(discovery_id=discovery_id, source_minimum=source_minimum,
                           mechanism="nt2", outcome=outcome, product_species=product)


def system(*endpoints):
    declared = [{"id": "hf", "xyz": "hf.xyz", "multiplicity": 1},
                {"id": "h2", "xyz": "h2.xyz", "multiplicity": 1},
                *({"id": e, "xyz": f"{e}.xyz", "multiplicity": 1, "role": "endpoint"}
                  for e in endpoints)]
    return SystemConfig(system_id="t", species=declared,
                        compositions=[{"id": "hf2", "components": {"hf": 2}},
                                      {"id": "hf_h2", "components": {"hf": 1, "h2": 1}}])


def test_reacting_candidates_are_the_screen_minima_discovery_ends_and_seeds_of_reacting_states():
    """(HF)2 reacts through its discoveries, HF as its monomer; H2 and the HF·H2 complex do not,
    unless the complex is a declared endpoint (its monomers HF and H2 then react too)."""
    dimer, mixed = ("H", "F", "H", "F"), ("H", "F", "H", "H")
    species = {s.species_id: s for s in (
        spc("hf", "FH_q0_m1", ("H", "F")), spc("h2", "H2_q0_m1", ("H", "H")),
        spc("d1", "F2H2_q0_m1", dimer), spc("d2", "F2H2_q0_m1", dimer),
        spc("d3", "F2H2_q0_m1", dimer), spc("x", "FH3_q0_m1", mixed),
        spc("p", "F2H2_q0_m1", dimer, energy=-199.0, state="p"))}
    screen = [screen_minimum("hf", "FH_q0_m1", -100.0),
              screen_minimum("d1", "F2H2_q0_m1", -200.0, members=("d2",)),
              screen_minimum("d3", "F2H2_q0_m1", -199.9), screen_minimum("x", "FH3_q0_m1", -50.0),
              screen_minimum("h2", "H2_q0_m1", -1.0)]
    discoveries = [found("to_member", "min_d3", "d2"),  # d2 is held by d1's basin
                   found("to_new", "min_d3", "p"),  # p skipped screen: it stands for itself
                   found("none", "min_x", outcome="negative")]  # a negative discovery: no end
    seed = species["p"].model_copy(update={"species_id": "spc_relax_q", "state_label": "q"})

    pool, idle = reacting_candidates(species, discoveries, screen, [seed], system())
    assert ids(pool) == ["hf", "d1", "d3", "p", "spc_relax_q"]
    assert [c.always for c in pool] == [False, True, True, True, True]
    assert [c.energy_hartree for c in pool] == [-100.0, -200.0, -199.9, -199.0, -199.0]
    assert idle == ["FH3_q0_m1", "H2_q0_m1"]

    pool, idle = reacting_candidates(species, discoveries, screen, [seed], system("x"))
    assert ids(pool) == ["hf", "d1", "d3", "x", "h2", "p", "spc_relax_q"] and idle == []

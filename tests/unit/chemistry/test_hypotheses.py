"""Reaction hypotheses (design §8.2): priority, degeneracy on basin structures (declared and
discovered), undeclared pairs only with a bond change (CH-07), one hypothesis per case key with
every distinct TS of it (G8-P7), R6, a discovery's own ends (X2: never a basin representative's
labelling), associations (X4) and the hypotheses closed by a static check (U5-P5)."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from fakes import NH3_HF, NH3_HF_EXCHANGED, NH3_HF_SYMBOLS

from hfauto.chemistry import identity, topology
from hfauto.chemistry.hypotheses import select
from hfauto.chemistry.xyz import XYZ
from hfauto.core import records as r
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import FileRef, Geometry, Level
from hfauto.core.system import ReactionInput

STORE: dict[str, XYZ] = {}
GEOMETRY: dict[str, Geometry] = {}  # species id -> its geometry, the basin structure in tests
REPRESENTATIVE: dict[str, str] = {}  # minimum id -> species id: a discovery's source (explore)
NH3 = np.array([[0.0, 0.0, 0.38], [0.94, 0.0, 0.0], [-0.47, 0.814, 0.0], [-0.47, -0.814, 0.0]])
CHFCLBR = np.array(
    [[0.0, 0.0, 0.0], [0.63, 0.63, 0.63], [-0.8, -0.8, 0.8], [-1.0, 1.0, -1.0], [1.1, -1.1, -1.1]]
)


def load(geometry: Geometry) -> XYZ:
    return STORE[geometry.fingerprint]


def species(sid: str, symbols: Sequence[str], coords) -> r.SpeciesRecord:
    STORE[sid] = XYZ(list(symbols), np.asarray(coords, dtype=float))
    geo = GEOMETRY[sid] = Geometry(file=FileRef(path=f"{sid}.xyz", sha256="0" * 64),
                                   fingerprint=sid, symbols=tuple(symbols))
    return r.SpeciesRecord(species_id=sid, composition_id="c", charge=0,
                           multiplicity=1, geometry=geo, source="input", state_label="x")


def minimum(mid: str, sid: str, energy: float = -1.0, members: tuple[str, ...] = ()
            ) -> r.MinimumRecord:
    REPRESENTATIVE[mid] = sid
    return r.MinimumRecord(minimum_id=mid, basin_id="b_" + mid, composition_id="c", species_id=sid,
                           tier="dft", level_key="L", opt_calc="o", freq_calc="f",
                           energy_hartree=energy, state_label="x", members=members)


def basins(*minima: r.MinimumRecord) -> list[tuple[r.MinimumRecord, Geometry]]:
    """Each minimum with its optimized structure: here its representative's geometry."""
    return [(m, GEOMETRY[m.species_id]) for m in minima]


def methanol(turn_deg: float = 0.0) -> np.ndarray:
    """C, O, H(O), then the methyl H staggered to H(O); turn_deg rotates the methyl rotor."""
    rows = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.43], [0.91, 0.0, 1.735]]
    for azimuth in np.radians(np.array([180.0, 300.0, 60.0]) + turn_deg):
        rows.append([1.028 * np.cos(azimuth), 1.028 * np.sin(azimuth), -0.364])
    return np.array(rows)


def product(did: str, source: str, sid: str, ts: Geometry | None = None) -> r.DiscoveryRecord:
    """An NT2 edge from the representative of the source minimum (as explore records it)."""
    return r.DiscoveryRecord(discovery_id=did, mechanism="nt2", outcome="product",
                             source_species=REPRESENTATIVE[source], product_species=sid, ts=ts,
                             generation=1)


def saddle(did: str, source: str, sid: str, ts_calc: str | None = "calc_saddle",
           side: str | None = None) -> r.DiscoveryRecord:
    """A mode-follow discovery from side 1 (default: the source minimum's representative) to
    side 2: ts_calc is its verified saddle at the DFT tier, None at screen."""
    return r.DiscoveryRecord(discovery_id=did, mechanism="mode_follow", outcome="product",
                             source_species=side or REPRESENTATIVE[source],
                             product_species=sid, ts=GEOMETRY[sid], ts_calc=ts_calc)


def h2o2(dihedral_deg: float) -> np.ndarray:
    t = np.radians(dihedral_deg)
    return np.array([[0.0, 0.0, 0.0], [1.45, 0.0, 0.0], [-0.2, 0.95, 0.0],
                     [1.65, 0.95 * np.cos(t), 0.95 * np.sin(t)]])


def test_declared_reaction_comes_first_and_borrows_the_discovery_ts():
    hcn = species("hcn", "HCN", [[0.0, 0, 0], [1.06, 0, 0], [2.22, 0, 0]])
    hnc = species("hnc", "HCN", [[3.22, 0, 0], [1.06, 0, 0], [2.22, 0, 0]])
    minima = [minimum("m_hcn", "hcn"), minimum("m_hnc", "hnc", -0.98)]
    negative = r.DiscoveryRecord(discovery_id="d2", mechanism="nt2", outcome="negative",
                                 reason="no_nt2_maximum", source_species="hcn")
    found = [product("d1", "m_hcn", "hnc", ts=hnc.geometry), negative]
    term = r.CoordinateTerm(kind="distance", atoms=(0, 2))
    iso = ReactionInput(id="iso", reactant="hcn", product="hnc", coordinate=[term])

    [rec] = select(basins(*minima), [hcn, hnc], found, [iso], load)
    assert (rec.reaction_id, rec.source) == ("iso", "declared")
    assert rec.minima == ("m_hcn", "m_hnc") and rec.endpoints == ("hcn", "hnc")
    assert rec.reactants[0].composition_id == rec.products[0].composition_id == "CHN_q0_m1"
    assert rec.coordinate == (term,) and not rec.torsional and not rec.degenerate
    assert rec.low_level_ts == (hnc.geometry,)  # d2 vetoes nothing

    [auto] = select(basins(*minima), [hcn, hnc], found, [], load)
    assert auto.source == "discovery" and auto.reaction_id.startswith("rxn_discovery_")
    assert auto.low_level_ts == (hnc.geometry,)
    [closed] = select(basins(*minima), [hcn, hnc], found, [], load, window_kcal=10.0)
    assert (closed.outcome, closed.reasons) == (r.CaseOutcome.OUT_OF_WINDOW, ("out_of_window",))
    assert closed.low_level_ts == (hnc.geometry,) and closed.log is None  # no job, no log
    unconnected = found[0].model_copy(update={"outcome": "unconnected", "generation": None})
    assert select(basins(*minima), [hcn, hnc], [unconnected], [], load) == []


def test_ammonia_inversion_is_a_degenerate_reaction_in_one_basin():
    up, down = species("nh3_up", "NHHH", NH3), species("nh3_down", "NHHH", NH3 * [1, 1, -1])
    copy = species("nh3_copy", "NHHH", NH3)
    basin = minimum("m_nh3", "nh3_up", members=("nh3_down", "nh3_copy"))
    inversion = ReactionInput(id="inv", reactant="nh3_up", product="nh3_down")
    same = ReactionInput(id="same", reactant="nh3_up", product="nh3_copy")

    inv, rel = select(basins(basin), [up, down, copy], [], [inversion, same], load)
    assert inv.degenerate and inv.minima == ("m_nh3", "m_nh3")
    assert inv.endpoints == ("nh3_up", "nh3_down") and inv.torsional
    assert not rel.degenerate and rel.minima == ("m_nh3", "m_nh3")
    assert (inv.outcome, rel.outcome, rel.reasons) == (None, r.CaseOutcome.SAME_BASIN,
                                                       ("same_basin",))


def test_rough_inputs_are_judged_on_the_basin_structure():
    """P5a (DME): a methyl rotation drawn by hand stays degenerate; the inputs alone are 0.06 A
    apart after relabelling, which failed the 0.05 A test on input coordinates."""
    noise = np.random.default_rng(0).normal(0.0, 0.05, (2, 6, 3))
    a = species("rotor_a", "COHHHH", methanol() + noise[0])
    b = species("rotor_b", "COHHHH", methanol(120.0) + noise[1])
    assert not identity.mapped_equivalent(a.geometry.symbols, STORE["rotor_a"].coords,
                                          STORE["rotor_b"].coords)
    basin = minimum("m_rotor", "rotor_a", members=("rotor_b",))
    rotation = ReactionInput(id="rot", reactant="rotor_a", product="rotor_b")
    [rec] = select([(basin, species("rotor_opt", "COHHHH", methanol()).geometry)], [a, b], [],
                   [rotation], load)
    assert rec.degenerate and rec.torsional


def test_a_declared_enantiomerization_is_degenerate_not_same_basin():
    r_form = species("r_form", ["C", "H", "F", "Cl", "Br"], CHFCLBR)
    s_form = species("s_form", ["C", "H", "F", "Cl", "Br"], CHFCLBR * [-1.0, 1.0, 1.0])
    basin = minimum("m_chfclbr", "r_form", members=("s_form",))
    [rec] = select(basins(basin), [r_form, s_form], [],
                   [ReactionInput(id="rs", reactant="r_form", product="s_form")], load)
    assert rec.degenerate and rec.minima == ("m_chfclbr", "m_chfclbr")


def test_a_bond_exchanging_discovery_within_one_basin_is_a_degenerate_hypothesis():
    """U4-P5 / S10: the NH3·HF double H exchange found by explore ends in the source's basin with
    H3 and H4 swapped, which exchanges bonded partners; it is evaluated as degenerate like a
    declared one. Found by mode-follow at the DFT tier it carries its verified saddle instead."""
    source = species("xch_src", NH3_HF_SYMBOLS, NH3_HF)
    swapped = species("xch_prod", NH3_HF_SYMBOLS, NH3_HF_EXCHANGED)
    screen = minimum("s_xch", "xch_src").model_copy(update={"tier": "screen"})
    basin = minimum("d_xch", "xch_src", members=("xch_prod",))
    ts = species("xch_ts", NH3_HF_SYMBOLS, (NH3_HF + NH3_HF_EXCHANGED) / 2).geometry
    found = [product("disc_xch", "s_xch", "xch_prod", ts=ts)]

    [rec] = select(basins(screen, basin), [source, swapped], found, [], load)
    assert (rec.source, rec.minima, rec.endpoints) == ("discovery", ("d_xch", "d_xch"),
                                                       ("xch_src", "xch_prod"))
    assert rec.degenerate and not rec.torsional
    assert (rec.low_level_ts, rec.ts_calc) == ((ts,), None)

    [dft] = select(basins(basin), [source, swapped], [saddle("mf", "d_xch", "xch_prod")], [], load)
    assert (dft.source, dft.low_level_ts, dft.ts_calc) == ("mode_follow", (), "calc_saddle")
    [xtb] = select(basins(screen, basin), [source, swapped],
                   [saddle("mf", "s_xch", "xch_prod", ts_calc=None)], [], load)
    assert (xtb.low_level_ts, xtb.ts_calc) == ((swapped.geometry,), None)  # an xTB saddle


def test_an_undeclared_torsion_or_enantiomerization_is_no_hypothesis():
    """S4: the mode-follow saddle of CH2OH. (Cs) joins two twisted mirror images of one basin,
    and a rotor turned by 120 degrees is a relabelling; no bond changes, so neither is a
    hypothesis unless declared (declared ones: the rotor and enantiomerization tests above)."""
    rot_a, rot_b = species("rot_a", "COHHHH", methanol()), species("rot_b", "COHHHH", methanol(120))
    mirror_r = species("mirror_r", ["C", "H", "F", "Cl", "Br"], CHFCLBR)
    mirror_s = species("mirror_s", ["C", "H", "F", "Cl", "Br"], CHFCLBR * [-1.0, 1.0, 1.0])
    rotor = minimum("m_rot", "rot_a", members=("rot_b",))
    chiral = minimum("m_mirror", "mirror_r", members=("mirror_s",))
    found = [saddle("mf_rot", "m_rot", "rot_b"), saddle("mf_mirror", "m_mirror", "mirror_s")]
    assert select(basins(rotor, chiral), [rot_a, rot_b, mirror_r, mirror_s], found, [], load) == []


def test_a_declared_inversion_borrows_the_verified_saddle_of_its_basin():
    """R3 with a planar seed: the D3h saddle's mode-follow sides fall into the NH3 basin. The
    inversion changes no bond, so the discovery is no hypothesis itself, but the declared case
    borrows its DFT saddle by basin, also after an explore TS was lent (the saddle is validated
    first)."""
    up, down = species("inv_up", "NHHH", NH3), species("inv_down", "NHHH", NH3 * [1, 1, -1])
    basin = basins(minimum("m_inv", "inv_up", members=("inv_down",)))
    inversion = ReactionInput(id="inv", reactant="inv_up", product="inv_down")
    follow = saddle("mf_inv", "m_inv", "inv_down")
    assert select(basin, [up, down], [follow], [], load) == []
    [rec] = select(basin, [up, down], [follow], [inversion], load)
    assert (rec.reaction_id, rec.ts_calc, rec.low_level_ts) == ("inv", "calc_saddle", ())
    explored = product("nt2_inv", "m_inv", "inv_down", ts=up.geometry)
    [rec] = select(basin, [up, down], [follow, explored], [inversion], load)
    assert (rec.ts_calc, rec.low_level_ts) == ("calc_saddle", (up.geometry,))


def test_a_discovery_mapped_onto_its_source_is_no_hypothesis():
    source = species("same_src", NH3_HF_SYMBOLS, NH3_HF)
    again = species("same_prod", NH3_HF_SYMBOLS, NH3_HF + 0.01)  # identity mapping
    screen = minimum("s_same", "same_src").model_copy(update={"tier": "screen"})
    basin = minimum("d_same", "same_src", members=("same_prod",))
    found = [product("disc_same", "s_same", "same_prod")]
    assert select(basins(screen, basin), [source, again], found, [], load) == []


def test_declared_endpoint_without_minimum_is_kept_for_blocking():
    hcn = species("hcn2", "HCN", [[0.0, 0, 0], [1.06, 0, 0], [2.22, 0, 0]])
    [rec] = select([], [hcn], [], [ReactionInput(id="r", reactant="hcn2", product="gone")], load)
    assert rec.minima == ("", "") and rec.products == ()
    assert (rec.outcome, rec.reasons) == (r.CaseOutcome.BLOCKED, ("endpoint_without_dft_minimum",))
    screen = minimum("s_hcn", "hcn2").model_copy(update={"tier": "screen"})
    [only] = select(basins(screen), [hcn], [], [ReactionInput(id="r", reactant="hcn2",
                                                               product="hcn2")], load)
    assert only.outcome is r.CaseOutcome.BLOCKED  # a screen minimum is no DFT minimum


def test_declared_endpoint_collapsed_at_screen_takes_the_dft_basin():
    """W7 TMA·(HF)2: an endpoint merged into a screen basin was BLOCKED (endpoints_not_on_one_pes)."""
    a = species("pt_a", "NHHH", NH3)
    b = species("pt_b", "NHHH", NH3 + [[0.0, 0.0, 0.2], [0, 0, 0], [0, 0, 0], [0, 0, 0]])
    screen = minimum("s_a", "pt_a", members=("pt_b",)).model_copy(update={"tier": "screen"})
    declared = ReactionInput(id="pt", reactant="pt_a", product="pt_b")
    [rec] = select(basins(screen, minimum("d_a", "pt_a")), [a, b], [], [declared], load)
    assert rec.minima == ("d_a", "d_a") and not rec.degenerate
    assert rec.outcome is r.CaseOutcome.SAME_BASIN


def test_conformers_of_one_state_are_no_hypothesis_but_a_declared_torsion_is():
    """U5-P7: three DFT minima of one state_label give no hypothesis (Curtin-Hammett); the
    declared trans -> cis rotation of HONO (S3) stays."""
    confs = [species(f"hooh{d}", "OOHH", h2o2(d)) for d in (180, 140, 100)]
    minima = [minimum(f"m{i}", s.species_id, -1.0 + 1e-4 * i) for i, s in enumerate(confs)]
    assert select(basins(*minima), confs, [], [], load) == []

    trans = species("hono_t", "HONO", [[-0.1975, -0.9436, 0.0], [0.0248, -0.0016, 0.0],
                                        [1.4070, 0.0083, 0.0], [1.8276, 1.0971, 0.0]])
    cis = species("hono_c", "HONO", [[-0.2479, 0.9422, 0.0], [0.0449, 0.0070, 0.0],
                                      [1.3960, 0.0144, 0.0], [1.8689, 1.0940, 0.0]])
    dihedral = r.CoordinateTerm(kind="dihedral", atoms=(0, 1, 2, 3))
    rotation = ReactionInput(id="trans_to_cis", reactant="hono_t", product="hono_c",
                             coordinate=[dihedral], torsional=True)
    hono = basins(minimum("m_t", "hono_t", -205.3455), minimum("m_c", "hono_c", -205.3453))
    [rec] = select(hono, [trans, cis], [], [rotation], load)
    assert (rec.reaction_id, rec.source, rec.minima) == ("trans_to_cis", "declared",
                                                         ("m_t", "m_c"))
    assert rec.torsional and not rec.degenerate


def test_an_undeclared_pair_needs_a_bond_change():
    """CH-07 on bonds: the proton moved across F-H...N (F-H broken) is a hypothesis; the proton
    0.05 A further out, or the H...N contact 0.3 A longer (0.2 A sufficed before), changes no
    bond and is none. F-H stretched to 0.05 A past r_thr (1.28 A) splits the bond graph inside
    the band (G8-P5): no bond change, so no hypothesis, and a declared one is torsional."""
    fhn = species("fhn", "FHN", [[0.0, 0, 0], [1.04, 0, 0], [2.28, 0, 0]])
    near = species("near", "FHN", [[0.0, 0, 0], [1.09, 0, 0], [2.28, 0, 0]])
    far = species("far", "FHN", [[0.0, 0, 0], [1.40, 0, 0], [2.28, 0, 0]])
    edge = species("edge", "FHN", [[0.0, 0, 0], [1.33, 0, 0], [2.28, 0, 0]])
    contact = species("contact", "FHN", [[0.0, 0, 0], [0.93, 0, 0], [2.83, 0, 0]])
    loose = species("loose", "FHN", [[0.0, 0, 0], [0.93, 0, 0], [3.13, 0, 0]])
    names = ("fhn", "near", "far", "edge", "contact", "loose")
    minima = basins(*(minimum(f"m_{n}", n) for n in names))
    found = [product("d_near", "m_fhn", "near"), product("d_far", "m_fhn", "far"),
             product("d_edge", "m_fhn", "edge"), product("d_loose", "m_contact", "loose")]
    everything = [fhn, near, far, edge, contact, loose]

    records = select(minima, everything, found, [], load)
    assert [(r.minima, r.torsional) for r in records] == [(("m_fhn", "m_far"), False)]
    assert topology.bonds(edge.geometry.symbols, STORE["edge"].coords) != topology.bonds(
        fhn.geometry.symbols, STORE["fhn"].coords)
    [stretch] = select(minima, everything, [], [ReactionInput(id="s", reactant="fhn",
                                                              product="edge")], load)
    assert stretch.torsional


def test_one_hypothesis_per_state_pair_holds_every_distinct_ts_of_it():
    """G8-P7 (W3/VAL7 S6): two discoveries reach two basins of one product state (H on N) from
    H on F. They make one hypothesis, the first's, which holds both distinct TSs in priority
    order; a third TS in one basin with the first adds none. Declared reactions of that state
    pair keep their own records and each takes every TS."""
    fh = species("g_fh", "FHN", [[0.0, 0, 0], [1.04, 0, 0], [2.28, 0, 0]])
    hn = [species("g_hn1", "FHN", [[0.0, 0, 0], [1.40, 0, 0], [2.28, 0, 0]]),
          species("g_hn2", "FHN", [[0.0, 0, 0], [1.50, 0.2, 0], [2.45, 0, 0]])]
    ts = [species(f"g_ts{k}", "FHN", x).geometry for k, x in enumerate((
        [[0.0, 0, 0], [1.22, 0, 0], [2.28, 0, 0]], [[0.0, 0, 0], [1.25, 0.4, 0], [2.40, 0, 0]],
        [[0.0, 0, 0], [1.23, 0.01, 0], [2.28, 0, 0]]))]
    minima = basins(minimum("g_mfh", "g_fh").model_copy(update={"state_label": "FH+N"}),
                    *(minimum(f"g_m{s.species_id[2:]}", s.species_id).model_copy(
                        update={"state_label": "HN+F"}) for s in hn))
    found = [product(f"g_d{k}", "g_mfh", sid, ts=t)
             for k, (sid, t) in enumerate(zip(("g_hn1", "g_hn2", "g_hn2"), ts, strict=True))]

    [rec] = select(minima, [fh, *hn], found, [], load)
    assert (rec.minima, rec.low_level_ts) == (("g_mfh", "g_mhn1"), (ts[0], ts[1]))
    declared = [ReactionInput(id="fwd", reactant="g_fh", product="g_hn1"),
                ReactionInput(id="back", reactant="g_hn2", product="g_fh")]
    out = select(minima, [fh, *hn], found, declared, load)
    assert [(r.reaction_id, r.low_level_ts) for r in out] == [("fwd", (ts[0], ts[1])),
                                                              ("back", (ts[0], ts[1]))]


def test_an_edge_without_a_ts_runs_while_its_source_keeps_its_state_at_dft():
    """R6 (S6): the seed (H at F) collapsed at screen into H at N and represents that basin. The
    edge runs from the unrelaxed seed as a species of its own, whose DFT job kept the seed state,
    to the collapse basin: one hypothesis, with no TS. When DFT collapsed it too, it joined the
    collapse basin: no hypothesis."""
    seed = species("seed", "FHN", [[0.0, 0, 0], [0.93, 0, 0], [2.83, 0, 0]]).model_copy(
        update={"state_label": "FH+N"})
    own = seed.model_copy(update={"species_id": "spc_relax_seed_0"})
    collapsed = species("collapsed", "FHN", [[0.0, 0, 0], [1.40, 0, 0], [2.28, 0, 0]]).geometry
    lost = minimum("s_lost", "seed").model_copy(update={"tier": "screen", "state_label": "HN+F"})
    basin = minimum("d_lost", "seed", -1.01).model_copy(update={"state_label": "HN+F"})
    kept = minimum("d_seed", "spc_relax_seed_0").model_copy(update={"state_label": "FH+N"})
    relax = r.DiscoveryRecord(discovery_id="relax_seed", mechanism="relaxation",
                              outcome="product", source_species="spc_relax_seed_0",
                              product_species="seed", generation=1)
    minima = [(lost, collapsed), (basin, collapsed), (kept, seed.geometry)]

    [rec] = select(minima, [seed, own], [relax], [], load)
    assert (rec.source, rec.minima, rec.endpoints) == ("discovery", ("d_seed", "d_lost"),
                                                       ("spc_relax_seed_0", "seed"))
    assert not rec.torsional and (rec.low_level_ts, rec.ts_calc) == ((), None)
    joined = basin.model_copy(update={"members": ("seed", "spc_relax_seed_0")})
    assert select([(lost, collapsed), (joined, collapsed)], [seed, own], [relax], [], load) == []


NH4_F = np.vstack([NH3_HF[:4], [[1.0, -0.2, 0.0]], NH3_HF[5:]])  # H4 moved from F to N


def change(a: str, b: str) -> tuple[frozenset, frozenset]:
    """The labelled bond change between two species' own structures."""
    return topology.bond_changes(STORE[a].symbols, STORE[a].coords, STORE[b].coords)


def test_a_discovery_keeps_its_own_ends_under_a_relabelled_representative():
    """X2 (G8-3): the DFT basins are represented by relabelled species (H3 and H4 swapped). The
    hypothesis runs between the discovery's own ends, so the case changes the bonds the discovery
    changed: in one basin the double H exchange stays a degenerate rearrangement (the
    representative's labelling made it a torsion, no hypothesis), and across two basins the
    proton moves from F5 to N0 as H4 (not as H3). A discovery without its own ends gives none."""
    source, swapped = species("own_src", NH3_HF_SYMBOLS, NH3_HF), species(
        "own_rep", NH3_HF_SYMBOLS, NH3_HF_EXCHANGED)
    exchanged = species("own_xch", NH3_HF_SYMBOLS, NH3_HF_EXCHANGED + 0.01)
    moved = species("own_pt", NH3_HF_SYMBOLS, NH4_F)
    screen = minimum("s_own", "own_src").model_copy(update={"tier": "screen"})
    basin = minimum("d_own", "own_rep", members=("own_src", "own_xch"))
    pair = minimum("d_pt", "own_pt", -0.99)
    everything = [source, swapped, exchanged, moved]

    [rec] = select(basins(screen, basin), everything, [product("xch", "s_own", "own_xch")], [],
                   load)
    assert (rec.minima, rec.endpoints) == (("d_own", "d_own"), ("own_src", "own_xch"))
    assert rec.degenerate and not rec.torsional

    [rec] = select(basins(screen, basin, pair), everything, [product("pt", "s_own", "own_pt")],
                   [], load)
    assert (rec.minima, rec.endpoints) == (("d_own", "d_pt"), ("own_src", "own_pt"))
    assert change("own_src", "own_pt") == ({(0, 4)}, {(4, 5)})  # the discovery's own change
    ends = [identity.basin_coords(NH3_HF_SYMBOLS, STORE["own_rep"].coords, STORE["own_src"].coords),
            STORE["own_pt"].coords]
    assert topology.bond_changes(NH3_HF_SYMBOLS, *ends) == change("own_src", "own_pt")

    anonymous = product("pt", "s_own", "own_pt").model_copy(update={"source_species": None})
    assert select(basins(screen, basin, pair), everything, [anonymous], [], load) == []


def test_a_mode_follow_hypothesis_runs_between_its_side_species():
    """The DFT saddle's sides (mf1, mf2) are the discovery's ends even when side 1 joined a basin
    whose representative is labelled otherwise; its verified saddle comes along."""
    rep = species("mf_rep", NH3_HF_SYMBOLS, NH3_HF_EXCHANGED)
    one, two = (species("t_mf1", NH3_HF_SYMBOLS, NH3_HF + 0.01),
                species("t_mf2", NH3_HF_SYMBOLS, NH4_F))
    minima = basins(minimum("m_one", "mf_rep", members=("t_mf1",)), minimum("m_two", "t_mf2"))
    [rec] = select(minima, [rep, one, two], [saddle("mf", "m_one", "t_mf2", side="t_mf1")], [],
                   load)
    assert (rec.source, rec.minima, rec.endpoints) == ("mode_follow", ("m_one", "m_two"),
                                                       ("t_mf1", "t_mf2"))
    assert (rec.ts_calc, rec.low_level_ts, rec.torsional) == ("calc_saddle", (), False)


# CH3OO (C O O H H H) and CH3·O2, the O2 moved 1.5 A out along C-O (C...O 2.95 A)
CH3OO = np.array([[0.0, 0, 0], [1.45, 0, 0], [1.905, 1.25, 0], [-0.363, 1.028, 0],
                  [-0.363, -0.514, 0.89], [-0.363, -0.514, -0.89]])
CH3_O2 = CH3OO + np.outer([0, 1, 1, 0, 0, 0], [1.5, 0, 0])
CH2_HOO = CH3_O2 + np.outer([0, 0, 0, 1, 0, 0], [3.768, 1.192, 0])  # H3 moved from C to O2
# 1,2-dioxetane (C C O O H H H H) and C2H4·O2, the O2 moved 1.5 A out
DIOXETANE = np.array([[0.0, 0, 0], [1.54, 0, 0], [0, 1.45, 0], [1.54, 1.45, 0],
                      [-0.35, -0.55, 0.85], [-0.35, -0.55, -0.85], [1.89, -0.55, 0.85],
                      [1.89, -0.55, -0.85]])
C2H4_O2 = DIOXETANE + np.outer([0, 0, 1, 1, 0, 0, 0, 0], [0, 1.5, 0])
# ethane (C C H H H H H H) and CH3·CH3, the second methyl moved 1.5 A out along C-C
ETHANE = np.array([[0.0, 0, 0], [1.53, 0, 0], [-0.36, 1.03, 0], [-0.36, -0.51, 0.89],
                   [-0.36, -0.51, -0.89], [1.89, -1.03, 0], [1.89, 0.51, 0.89],
                   [1.89, 0.51, -0.89]])
CH3_CH3 = ETHANE + np.outer([0, 1, 0, 0, 0, 1, 1, 1], [1.5, 0, 0])
# monomer states as thermo.separated_states reads a complex's fragments (topology.state_label)
CH3, O2, C2H4 = ((f, topology.state_label(list(s), x)) for f, s, x in (
    ("CH3", "CHHH", CH3OO[[0, 3, 4, 5]]), ("O2", "OO", CH3OO[[1, 2]]),
    ("C2H4", "CCHHHH", DIOXETANE[[0, 1, 4, 5, 6, 7]])))
MONOMERS = {("CH3O2", 0): [(CH3, O2)], ("C2H4O2", 0): [(C2H4, O2)],  # thermo.declared_monomers
            ("C2H6", 0): [(CH3, CH3)]}


def level(multiplicity: int, basis: str = "def2-svpd") -> Level:
    return Level(program="nwchem", version="7.2.3", method="pbe0", basis=basis, charge=0,
                 multiplicity=multiplicity)


def monomer(mid: str, state: tuple[str, str], energy: float) -> r.MinimumRecord:
    """A DFT minimum of a monomer ``state`` (its composition and state label)."""
    species(mid, "H", [[0.0, 0.0, 0.0]])  # its structure is never read
    return minimum(mid, mid, energy).model_copy(update={"composition_id": state[0],
                                                        "state_label": state[1]})


def test_a_one_bond_formation_between_fragments_is_an_association():
    """X4 (G2-P4): CH3·O2 -> CH3OO forms one C-O bond between the fragments and breaks none: its
    reactant side is the separated CH3 and O2, each the lowest DFT minimum of its state on the
    complex's level (charge and multiplicity aside: a lower CH3 on another basis is none), the
    complex staying minima[0]; CH3·CH3 -> C2H6 takes one CH3 twice. The H abstraction CH3·O2 ->
    CH2·HOO (a bond broken), the cycloaddition C2H4·O2 -> dioxetane (two bonds formed) and an
    association without a monomer minimum on its level stay ordinary hypotheses."""
    ends = [species(s, symbols, x) for s, symbols, x in (
        ("cpx", "COOHHH", CH3_O2), ("add", "COOHHH", CH3OO), ("abst", "COOHHH", CH2_HOO),
        ("c2o2", "CCOOHHHH", C2H4_O2), ("ring", "CCOOHHHH", DIOXETANE),
        ("ch3ch3", "CCHHHHHH", CH3_CH3), ("c2h6", "CCHHHHHH", ETHANE))]
    pairs = [minimum(f"m_{e.species_id}", e.species_id, -1.0 - 0.01 * k)
             for k, e in enumerate(ends)]
    parts = [monomer("m_ch3", CH3, -0.5), monomer("m_ch3_tz", CH3, -0.6),
             monomer("m_o2", O2, -0.4), monomer("m_c2h4", C2H4, -0.3)]
    levels = {m.minimum_id: level(2) for m in pairs} | {
        "m_ch3": level(2), "m_ch3_tz": level(2, "def2-tzvpd"), "m_o2": level(3),
        "m_c2h4": level(1)}
    declared = [ReactionInput(id=i, reactant=a, product=b) for i, a, b in (
        ("assoc", "cpx", "add"), ("abst", "cpx", "abst"), ("ring", "c2o2", "ring"),
        ("dimer", "ch3ch3", "c2h6"))]

    def run(**kw) -> dict[str, r.ReactionRecord]:
        found = select(basins(*pairs, *parts), ends, [], declared, load,
                       **({"monomers": MONOMERS, "levels": levels} | kw))
        return {rec.reaction_id: rec for rec in found}

    out = run()
    assoc, dimer = out["assoc"], out["dimer"]
    assert (assoc.monomers, assoc.minima) == (("m_ch3", "m_o2"), ("m_cpx", "m_add"))
    assert [(t.composition_id, t.coefficient) for t in assoc.reactants] == [("CH3", 1),
                                                                         ("O2", 1)]
    assert assoc.products[0].composition_id == "CH3O2_q0_m1" and not assoc.torsional
    assert dimer.monomers == ("m_ch3", "m_ch3")
    assert [(t.composition_id, t.coefficient) for t in dimer.reactants] == [("CH3", 2)]
    assert change("cpx", "abst") == ({(2, 3)}, {(0, 3)})
    assert change("c2o2", "ring") == ({(0, 2), (1, 3)}, set())
    for ordinary in (out["abst"], out["ring"]):
        assert ordinary.monomers == () and ordinary.reactants == ordinary.products
    elsewhere = run(levels=levels | {"m_o2": level(3, "def2-tzvpd")})["assoc"]
    assert elsewhere.monomers == () and elsewhere.reactants == elsewhere.products
    assert run(monomers={})["assoc"].monomers == ()


def test_a_complex_that_relaxed_into_its_adduct_is_judged_on_its_own_structure():
    """X4: a declared complex without a DFT minimum of its own (BH3 + NH3 relaxes into the
    adduct) is still an association, judged on its input structure; minima[0] is then the
    adduct's basin, and the record is neither degenerate nor torsional. Without monomer states
    it stays an ordinary hypothesis."""
    ends = [species("cpx_in", "COOHHH", CH3_O2), species("add_in", "COOHHH", CH3OO)]
    adduct = minimum("m_add_in", "add_in", -1.1, members=("cpx_in",))
    parts = [monomer("m_ch3", CH3, -0.5), monomer("m_o2", O2, -0.4)]
    levels = {"m_add_in": level(2), "m_ch3": level(2), "m_o2": level(3)}
    declared = [ReactionInput(id="assoc", reactant="cpx_in", product="add_in")]
    [rec] = select(basins(adduct, *parts), ends, [], declared, load, monomers=MONOMERS,
                   levels=levels)
    assert (rec.monomers, rec.minima) == (("m_ch3", "m_o2"), ("m_add_in", "m_add_in"))
    assert not rec.degenerate and not rec.torsional
    [plain] = select(basins(adduct, *parts), ends, [], declared, load)
    assert plain.monomers == () and plain.reactants == plain.products
    assert rec.outcome is None and plain.outcome is r.CaseOutcome.SAME_BASIN


def test_the_window_is_measured_from_the_reactant_asymptote():
    """U5-P5: one window function: an association's product lies 25 kcal/mol above its
    separated monomers (inside the 40 kcal/mol window) though 56 above its complex; without the
    monomers the complex is the reactant and the hypothesis closes out of the window."""
    ends = [species("cpx_w", "COOHHH", CH3_O2), species("add_w", "COOHHH", CH3OO)]
    k = 1.0 / HARTREE_TO_KCAL_MOL
    pairs = [minimum("m_cpx_w", "cpx_w", -0.9 - 31 * k), minimum("m_add_w", "add_w", -0.9 + 25 * k)]
    parts = [monomer("m_ch3", CH3, -0.5), monomer("m_o2", O2, -0.4)]
    levels = {"m_cpx_w": level(2), "m_add_w": level(2), "m_ch3": level(2), "m_o2": level(3)}
    declared = [ReactionInput(id="assoc", reactant="cpx_w", product="add_w")]
    [assoc] = select(basins(*pairs, *parts), ends, [], declared, load, monomers=MONOMERS,
                     levels=levels)
    assert assoc.monomers == ("m_ch3", "m_o2") and assoc.outcome is None
    [plain] = select(basins(*pairs, *parts), ends, [], declared, load)
    assert plain.outcome is r.CaseOutcome.OUT_OF_WINDOW

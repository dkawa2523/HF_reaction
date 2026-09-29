"""Reaction hypotheses (design §8.2): priority, degeneracy on basin structures (declared and
discovered), undeclared pairs only with a bond change (CH-07), TS lending, chem 13, R6."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from fakes import NH3_HF, NH3_HF_EXCHANGED, NH3_HF_SYMBOLS

from hfauto.chemistry import identity
from hfauto.chemistry.hypotheses import pick_endpoints, select
from hfauto.chemistry.xyz import XYZ
from hfauto.core import records as r
from hfauto.core.evidence import FileRef, Geometry
from hfauto.core.system import ReactionInput

STORE: dict[str, XYZ] = {}
GEOMETRY: dict[str, Geometry] = {}  # species id -> its geometry, the basin structure in tests
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
    return r.DiscoveryRecord(discovery_id=did, source_minimum=source, mechanism="nt2",
                             outcome="product", product_species=sid, ts=ts)


def saddle(did: str, source: str, sid: str, ts_calc: str | None = "calc_saddle"
           ) -> r.DiscoveryRecord:
    """A mode-follow discovery: ts_calc is its verified saddle at the DFT tier, None at screen."""
    return r.DiscoveryRecord(discovery_id=did, source_minimum=source, mechanism="mode_follow",
                             outcome="product", product_species=sid, ts=GEOMETRY[sid],
                             ts_calc=ts_calc)


def h2o2(dihedral_deg: float) -> np.ndarray:
    t = np.radians(dihedral_deg)
    return np.array([[0.0, 0.0, 0.0], [1.45, 0.0, 0.0], [-0.2, 0.95, 0.0],
                     [1.65, 0.95 * np.cos(t), 0.95 * np.sin(t)]])


def test_declared_reaction_comes_first_and_borrows_the_discovery_ts():
    hcn = species("hcn", "HCN", [[0.0, 0, 0], [1.06, 0, 0], [2.22, 0, 0]])
    hnc = species("hnc", "HCN", [[3.22, 0, 0], [1.06, 0, 0], [2.22, 0, 0]])
    minima = [minimum("m_hcn", "hcn"), minimum("m_hnc", "hnc", -0.98)]
    negative = r.DiscoveryRecord(discovery_id="d2", source_minimum="m_hcn", mechanism="nt2",
                                 outcome="negative", reason="no_nt2_maximum")
    found = [product("d1", "m_hcn", "hnc", ts=hnc.geometry), negative]
    term = r.CoordinateTerm(kind="distance", atoms=(0, 2))
    iso = ReactionInput(id="iso", reactant="hcn", product="hnc", coordinate=[term])

    [rec] = select(basins(*minima), [hcn, hnc], found, [iso], load)
    assert (rec.reaction_id, rec.source) == ("iso", "declared")
    assert rec.minima == ("m_hcn", "m_hnc") and rec.endpoints == ("hcn", "hnc")
    assert rec.reactants[0].composition_id == rec.products[0].composition_id == "CHN_q0_m1"
    assert rec.coordinate == (term,) and not rec.torsional and not rec.degenerate
    assert rec.low_level_ts == hnc.geometry  # d2 vetoes nothing

    [auto] = select(basins(*minima), [hcn, hnc], found, [], load)
    assert auto.source == "discovery" and auto.reaction_id.startswith("rxn_discovery_")
    assert auto.low_level_ts == hnc.geometry
    assert select(basins(*minima), [hcn, hnc], found, [], load, window_kcal=10.0) == []


def test_ammonia_inversion_is_a_degenerate_reaction_in_one_basin():
    up, down = species("nh3_up", "NHHH", NH3), species("nh3_down", "NHHH", NH3 * [1, 1, -1])
    copy = species("nh3_copy", "NHHH", NH3)
    basin = minimum("m_nh3", "nh3_up", members=("nh3_down", "nh3_copy"))
    inversion = ReactionInput(id="inv", reactant="nh3_up", product="nh3_down")
    same = ReactionInput(id="same", reactant="nh3_up", product="nh3_copy")

    inv, rel = select(basins(basin), [up, down, copy], [], [inversion, same], load)
    assert inv.degenerate and inv.minima == ("m_nh3", "m_nh3")
    assert inv.endpoints == ("nh3_up", "nh3_down") and inv.torsional
    assert not rel.degenerate and rel.minima == ("m_nh3", "m_nh3")  # decide(): SAME_BASIN


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
    assert (rec.low_level_ts, rec.ts_calc) == (ts, None)

    [dft] = select(basins(basin), [source, swapped], [saddle("mf", "d_xch", "xch_prod")], [], load)
    assert (dft.source, dft.low_level_ts, dft.ts_calc) == ("mode_follow", None, "calc_saddle")
    [xtb] = select(basins(screen, basin), [source, swapped],
                   [saddle("mf", "s_xch", "xch_prod", ts_calc=None)], [], load)
    assert (xtb.low_level_ts, xtb.ts_calc) == (swapped.geometry, None)  # an xTB saddle


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
    assert (rec.reaction_id, rec.ts_calc, rec.low_level_ts) == ("inv", "calc_saddle", None)
    explored = product("nt2_inv", "m_inv", "inv_down", ts=up.geometry)
    [rec] = select(basin, [up, down], [follow, explored], [inversion], load)
    assert (rec.ts_calc, rec.low_level_ts) == ("calc_saddle", up.geometry)


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


def test_declared_endpoint_collapsed_at_screen_takes_the_dft_basin():
    """W7 TMA·(HF)2: an endpoint merged into a screen basin was BLOCKED (endpoints_not_on_one_pes)."""
    a = species("pt_a", "NHHH", NH3)
    b = species("pt_b", "NHHH", NH3 + [[0.0, 0.0, 0.2], [0, 0, 0], [0, 0, 0], [0, 0, 0]])
    screen = minimum("s_a", "pt_a", members=("pt_b",)).model_copy(update={"tier": "screen"})
    declared = ReactionInput(id="pt", reactant="pt_a", product="pt_b")
    [rec] = select(basins(screen, minimum("d_a", "pt_a")), [a, b], [], [declared], load)
    assert rec.minima == ("d_a", "d_a") and not rec.degenerate  # decide(): SAME_BASIN


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
    bond and is none."""
    fhn = species("fhn", "FHN", [[0.0, 0, 0], [1.04, 0, 0], [2.28, 0, 0]])
    near = species("near", "FHN", [[0.0, 0, 0], [1.09, 0, 0], [2.28, 0, 0]])
    far = species("far", "FHN", [[0.0, 0, 0], [1.40, 0, 0], [2.28, 0, 0]])
    contact = species("contact", "FHN", [[0.0, 0, 0], [0.93, 0, 0], [2.83, 0, 0]])
    loose = species("loose", "FHN", [[0.0, 0, 0], [0.93, 0, 0], [3.13, 0, 0]])
    names = ("fhn", "near", "far", "contact", "loose")
    found = [product("d_near", "m_fhn", "near"), product("d_far", "m_fhn", "far"),
             product("d_loose", "m_contact", "loose")]

    records = select(basins(*(minimum(f"m_{n}", n) for n in names)),
                     [fhn, near, far, contact, loose], found, [], load)
    assert [(r.minima, r.torsional) for r in records] == [(("m_fhn", "m_far"), False)]


def test_a_relaxation_seed_kept_at_dft_runs_to_its_collapse_basin():
    """R6 (S6): the seed (H at F) collapsed at screen into H at N and represents that basin. Its
    own DFT job (seed_species_id) kept the seed state: one hypothesis, seed -> collapse, with no
    TS. When DFT collapsed it too, it joined the collapse basin: no hypothesis."""
    seed = species("seed", "FHN", [[0.0, 0, 0], [0.93, 0, 0], [2.83, 0, 0]]).model_copy(
        update={"state_label": "FH+N"})
    own = seed.model_copy(update={"species_id": "spc_relax_seed"})
    collapsed = species("collapsed", "FHN", [[0.0, 0, 0], [1.40, 0, 0], [2.28, 0, 0]]).geometry
    lost = minimum("s_lost", "seed").model_copy(update={"tier": "screen", "state_label": "HN+F"})
    basin = minimum("d_lost", "seed", -1.01).model_copy(update={"state_label": "HN+F"})
    kept = minimum("d_seed", "spc_relax_seed").model_copy(update={"state_label": "FH+N"})
    relax = r.DiscoveryRecord(discovery_id="relax_seed", source_minimum="s_lost",
                              mechanism="relaxation", outcome="product", product_species="seed")
    minima = [(lost, collapsed), (basin, collapsed), (kept, seed.geometry)]

    [rec] = select(minima, [seed, own], [relax], [], load)
    assert (rec.source, rec.minima, rec.endpoints) == ("discovery", ("d_seed", "d_lost"),
                                                       ("spc_relax_seed", "seed"))
    assert not rec.torsional and (rec.low_level_ts, rec.ts_calc) == (None, None)
    assert select(minima, [seed, own], [relax], [], load, max_per_composition=0) == []
    joined = basin.model_copy(update={"members": ("seed", "spc_relax_seed")})
    assert select([(lost, collapsed), (joined, collapsed)], [seed, own], [relax], [], load) == []


def test_pick_endpoints_avoids_a_permuted_representative():
    relabelled = species("rep", "NHHH", NH3[[0, 2, 1, 3]])
    natural = species("nat", "NHHH", NH3)
    target = species("tgt", "NHHH", NH3 + [[0.0, 0.0, 0.05], [0, 0, 0], [0, 0, 0], [0, 0, 0]])
    assert pick_endpoints([relabelled, natural], [target], load) == ("nat", "tgt")
    assert pick_endpoints([relabelled], [species("hf", "HF", [[0, 0, 0], [0.92, 0, 0]])],
                          load) is None

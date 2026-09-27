"""Reaction hypotheses (design §8.2): priority, degeneracy on basin structures (declared and
discovered), no conformer pairs, CH-07, chem 13."""

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


def h2o2(dihedral_deg: float) -> np.ndarray:
    t = np.radians(dihedral_deg)
    return np.array([[0.0, 0.0, 0.0], [1.45, 0.0, 0.0], [-0.2, 0.95, 0.0],
                     [1.65, 0.95 * np.cos(t), 0.95 * np.sin(t)]])


def test_declared_reaction_comes_first_and_borrows_the_discovery_ts():
    hcn = species("hcn", "HCN", [[0.0, 0, 0], [1.06, 0, 0], [2.22, 0, 0]])
    hnc = species("hnc", "HCN", [[3.22, 0, 0], [1.06, 0, 0], [2.22, 0, 0]])
    minima = [minimum("m_hcn", "hcn"), minimum("m_hnc", "hnc", -0.98)]
    negative = r.DiscoveryRecord(discovery_id="d2", source_minimum="m_hcn", mechanism="afir",
                                 outcome="negative", reason="monotonic_uphill")
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


def test_a_discovery_within_one_basin_is_a_degenerate_hypothesis():
    """U4-P5: the NH3·HF double H exchange found by explore ends in the source's basin with H3
    and H4 swapped; it is evaluated as degenerate like a declared one."""
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
    assert rec.low_level_ts == ts


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


def test_a_0p05_angstrom_change_is_not_a_reaction():
    """CH-07: the legacy connect-minima test accepted this pair; the expectation is inverted."""
    fhn = species("fhn", "FHN", [[0.0, 0, 0], [1.04, 0, 0], [2.28, 0, 0]])
    near = species("near", "FHN", [[0.0, 0, 0], [1.09, 0, 0], [2.28, 0, 0]])
    far = species("far", "FHN", [[0.0, 0, 0], [1.40, 0, 0], [2.28, 0, 0]])
    minima = [minimum("m_fhn", "fhn"), minimum("m_near", "near"), minimum("m_far", "far")]
    found = [product("d_near", "m_fhn", "near"), product("d_far", "m_fhn", "far")]

    records = select(basins(*minima), [fhn, near, far], found, [], load)
    assert [r.minima for r in records] == [("m_fhn", "m_far")]


def test_pick_endpoints_avoids_a_permuted_representative():
    relabelled = species("rep", "NHHH", NH3[[0, 2, 1, 3]])
    natural = species("nat", "NHHH", NH3)
    target = species("tgt", "NHHH", NH3 + [[0.0, 0.0, 0.05], [0, 0, 0], [0, 0, 0], [0, 0, 0]])
    assert pick_endpoints([relabelled, natural], [target], load) == ("nat", "tgt")
    assert pick_endpoints([relabelled], [species("hf", "HF", [[0, 0, 0], [0.92, 0, 0]])],
                          load) is None

"""Reaction hypotheses (design §8.2): priority, degeneracy, conformer pairs, CH-07, chem 13."""

from __future__ import annotations

import numpy as np

from hfauto.chemistry.hypotheses import pick_endpoints, select
from hfauto.chemistry.xyz import XYZ
from hfauto.core import records as r
from hfauto.core.evidence import FileRef, Geometry
from hfauto.core.system import ReactionInput

STORE: dict[str, XYZ] = {}
NH3 = np.array([[0.0, 0.0, 0.38], [0.94, 0.0, 0.0], [-0.47, 0.814, 0.0], [-0.47, -0.814, 0.0]])


def load(geometry: Geometry) -> XYZ:
    return STORE[geometry.fingerprint]


def species(sid: str, symbols: str, coords) -> r.SpeciesRecord:
    STORE[sid] = XYZ(list(symbols), np.asarray(coords, dtype=float))
    geo = Geometry(file=FileRef(path=f"{sid}.xyz", sha256="0" * 64), fingerprint=sid,
                   symbols=tuple(symbols))
    return r.SpeciesRecord(species_id=sid, composition_id="c", formula="f", charge=0,
                           multiplicity=1, geometry=geo, source="input", state_label="x")


def minimum(mid: str, sid: str, energy: float = -1.0, members: tuple[str, ...] = ()
            ) -> r.MinimumRecord:
    return r.MinimumRecord(minimum_id=mid, basin_id="b_" + mid, composition_id="c", species_id=sid,
                           tier="dft", level_key="L", opt_calc="o", freq_calc="f",
                           energy_hartree=energy, state_label="x", n_fragments=1, members=members)


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
    trial = r.ReactionTrial(trial_id="t", source_minimum="m_hcn", kind="h_shift", mechanism="afir",
                            associations=((2, 0),), dissociations=((0, 1),))
    negative = r.DiscoveryRecord(discovery_id="d2", source_minimum="m_hcn", mechanism="afir",
                                 outcome="negative", reason="monotonic_uphill", trial=trial)
    found = [product("d1", "m_hcn", "hnc", ts=hnc.geometry), negative]
    term = r.CoordinateTerm(kind="distance", atoms=(0, 2))
    iso = ReactionInput(id="iso", reactant="hcn", product="hnc", coordinate=[term])

    [rec] = select(minima, [hcn, hnc], found, [iso], load)
    assert (rec.reaction_id, rec.source) == ("iso", "declared")
    assert rec.minima == ("m_hcn", "m_hnc") and rec.endpoints == ("hcn", "hnc")
    assert rec.reactants[0].composition_id == rec.products[0].composition_id == "CHN_q0_m1"
    assert rec.coordinate == (term,) and not rec.torsional and not rec.degenerate
    assert rec.n_h_transferred == 1 and rec.low_level_ts == hnc.geometry
    assert rec.negative_evidence == ("d2:monotonic_uphill",)

    [auto] = select(minima, [hcn, hnc], found, [], load)
    assert auto.source == "discovery" and auto.reaction_id.startswith("rxn_discovery_")
    assert auto.low_level_ts == hnc.geometry
    assert select(minima, [hcn, hnc], found, [], load, window_kcal=10.0) == []


def test_ammonia_inversion_is_a_degenerate_reaction_in_one_basin():
    up, down = species("nh3_up", "NHHH", NH3), species("nh3_down", "NHHH", NH3 * [1, 1, -1])
    copy = species("nh3_copy", "NHHH", NH3)
    basin = minimum("m_nh3", "nh3_up", members=("nh3_down", "nh3_copy"))
    inversion = ReactionInput(id="inv", reactant="nh3_up", product="nh3_down")
    same = ReactionInput(id="same", reactant="nh3_up", product="nh3_copy")

    inv, rel = select([basin], [up, down, copy], [], [inversion, same], load)
    assert inv.degenerate and inv.minima == ("m_nh3", "m_nh3")
    assert inv.endpoints == ("nh3_up", "nh3_down") and inv.torsional
    assert not rel.degenerate and rel.minima == ("m_nh3", "m_nh3")  # decide(): SAME_BASIN


def test_declared_endpoint_without_minimum_is_kept_for_blocking():
    hcn = species("hcn2", "HCN", [[0.0, 0, 0], [1.06, 0, 0], [2.22, 0, 0]])
    [rec] = select([], [hcn], [], [ReactionInput(id="r", reactant="hcn2", product="gone")], load)
    assert rec.minima == ("", "") and rec.products == ()


def test_conformer_pairs_need_a_30_degree_twist():
    confs = [species(f"hooh{d}", "OOHH", h2o2(d)) for d in (180, 170, 110)]
    minima = [minimum(f"m{i}", s.species_id, -1.0 + 1e-4 * i) for i, s in enumerate(confs)]

    records = select(minima, confs, [], [], load)
    assert [r.minima for r in records] == [("m0", "m2"), ("m1", "m2")]  # 10° apart is dropped
    assert all(r.source == "conformer" and r.torsional and r.n_h_transferred == 0
               for r in records)


def test_a_0p05_angstrom_change_is_not_a_reaction():
    """CH-07: the legacy connect-minima test accepted this pair; the expectation is inverted."""
    fhn = species("fhn", "FHN", [[0.0, 0, 0], [1.04, 0, 0], [2.28, 0, 0]])
    near = species("near", "FHN", [[0.0, 0, 0], [1.09, 0, 0], [2.28, 0, 0]])
    far = species("far", "FHN", [[0.0, 0, 0], [1.40, 0, 0], [2.28, 0, 0]])
    minima = [minimum("m_fhn", "fhn"), minimum("m_near", "near"), minimum("m_far", "far")]
    found = [product("d_near", "m_fhn", "near"), product("d_far", "m_fhn", "far")]

    records = select(minima, [fhn, near, far], found, [], load)
    assert [r.minima for r in records] == [("m_fhn", "m_far")]


def test_pick_endpoints_avoids_a_permuted_representative():
    relabelled = species("rep", "NHHH", NH3[[0, 2, 1, 3]])
    natural = species("nat", "NHHH", NH3)
    target = species("tgt", "NHHH", NH3 + [[0.0, 0.0, 0.05], [0, 0, 0], [0, 0, 0], [0, 0, 0]])
    assert pick_endpoints([relabelled, natural], [target], load) == ("nat", "tgt")
    assert pick_endpoints([relabelled], [species("hf", "HF", [[0, 0, 0], [0.92, 0, 0]])],
                          load) is None

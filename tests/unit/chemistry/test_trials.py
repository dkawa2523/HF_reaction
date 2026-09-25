"""Reaction trials (design §8.2 explore): tier order, relays, H shifts, linear perturbation."""

import numpy as np

from hfauto.chemistry import trials
from hfauto.chemistry.vibrations import external_basis

HCN = (["H", "C", "N"], np.array([[0, 0, -1.066], [0, 0, 0], [0, 0, 1.156]]))
H2CO_WATER = (["C", "O", "H", "H", "O", "H", "H"],  # H-O-H ··· O=CH2 hydrogen bond
              np.array([[0, 0, 0], [0, 0, 1.21], [0, 0.94, -0.54], [0, -0.94, -0.54],
                        [0, 0, 4.1], [0, 0, 3.14], [0.93, 0, 4.34]]))
CHAIN = (["N", "H", "H", "H", "H", "F", "H", "F"],  # H3N ··· H-F ··· H-F (legacy relay case)
         np.array([[0, 0, 0], [-0.33, 0.94, 0], [-0.33, -0.47, 0.81], [-0.33, -0.47, -0.81],
                   [1.7, 0, 0], [2.62, 0, 0], [4.3, 0, 0], [5.22, 0, 0]]))
TIERS = ("polar_h", "h_shift", "heavy_bond", "association")


def _drives(ts):
    return [(t.kind, set(t.associations), set(t.dissociations)) for t in ts]


def test_four_tiers_in_priority_order_and_capped():
    _, ts = trials.generate("m", *H2CO_WATER, max_trials=50)
    kinds = [t.kind for t in ts]
    assert set(kinds) == set(TIERS) and kinds == sorted(kinds, key=TIERS.index)
    assert ("h_shift", {(1, 2)}, {(0, 2)}) in _drives(ts)  # formaldehyde 1,2-H shift C -> O
    assert ("association", {(1, 4)}, set()) in _drives(ts)  # closest heavy pair O···O
    assert len({t.trial_id for t in ts}) == len(ts) > 3 == len(trials.generate(
        "m", *H2CO_WATER, max_trials=3)[1])


def test_relay_comes_before_single_transfers():
    _, ts = trials.generate("m", *CHAIN)
    assert _drives(ts)[0] == ("polar_h", {(0, 4), (5, 6)}, {(4, 5), (6, 7)})
    assert all(t.kind == "polar_h" for t in ts[:4]) and not ts[0].perturbed


def test_linear_hcn_is_bent_and_gets_the_ch_shift():
    start, ts = trials.generate("m", *HCN)
    assert ("h_shift", {(0, 2)}, {(0, 1)}) in _drives(ts) and all(t.perturbed for t in ts)
    assert external_basis(HCN[0], start).shape[1] == 6
    v1, v2 = start[0] - start[1], start[2] - start[1]
    bend = 180 - np.degrees(np.arccos(v1 @ v2 / np.linalg.norm(v1) / np.linalg.norm(v2)))
    assert 5 < bend < 15 and np.abs(start - HCN[1]).max() < 0.3
    assert np.array_equal(start, trials.generate("m", *HCN)[0])  # seeded


def test_aromatic_ring_bonds_are_never_broken():
    ring = np.array([(np.cos(a), np.sin(a), 0.0) for a in np.radians(np.arange(0, 360, 60))])
    _, ts = trials.generate("m", ["C"] * 6 + ["H"] * 6, np.vstack([1.39 * ring, 2.47 * ring]),
                            max_trials=100)
    assert all(not t.dissociations for t in ts if t.kind == "heavy_bond")


def test_product_verdict_window_and_ts_requirement():
    def verdict(mechanism, barrier, reaction, ts=True):
        return trials.product_verdict(mechanism, ts_validated=ts, barrier_kj=barrier,
                                      reaction_kj=reaction)

    assert verdict("nt2", 10, 5, ts=False) == "ts_not_validated"
    assert verdict("nt2", 151, 5) == verdict("afir", None, 101) == "out_of_window"
    assert verdict("afir", None, None) == "out_of_window" and verdict("afir", None, 99) is None


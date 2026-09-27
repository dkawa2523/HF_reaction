"""Reaction trials (design §8.2 explore): the four bond-change templates, the valence rule,
drive equivalence classes and the linear perturbation."""

import numpy as np

from hfauto.chemistry import trials
from hfauto.chemistry.vibrations import external_basis

HCN = (["H", "C", "N"], np.array([[0, 0, -1.066], [0, 0, 0], [0, 0, 1.156]]))
_CH4 = [[0, 0, 0], [0, 0, 1.09], [1.028, 0, -0.363], [-0.514, 0.89, -0.363],
        [-0.514, -0.89, -0.363]]
SN2 = (["C", "H", "H", "H", "Cl", "Cl"],  # Cl⁻ ··· CH3Cl ion-dipole complex (backside)
       np.array([[0, 0, 0], [1.034, 0, -0.346], [-0.517, 0.8955, -0.346],
                 [-0.517, -0.8955, -0.346], [0, 0, 1.8], [0, 0, -3.2]]))
CH4_OH = (["C", "H", "H", "H", "H", "O", "H"], np.array(_CH4 + [[0, 0, 3.39], [0.94, 0, 3.63]]))
H_CH4 = (["C", "H", "H", "H", "H", "H"], np.array(_CH4 + [[0, 0, 2.89]]))
CH3O = (["C", "O", "H", "H", "H"], np.array([[0, 0, 0], [0, 0, 1.37], [1.05, 0, -0.39],
                                             [-0.517, 0.8955, -0.376], [-0.517, -0.8955, -0.376]]))
HONO = (["H", "O", "N", "O"], np.array([[-0.2017, 0.9487, 0], [0, 0, 0], [1.43, 0, 0],
                                        [1.8336, 1.1089, 0]]))
_TMA = [[0.682, -1.219, -0.115], [0.107, 0, 0.45], [-1.352, 0.017, 0.332], [0.711, 1.202, -0.122],
        [1.766, -1.223, 0.051], [0.249, -2.091, 0.389], [0.482, -1.29, -1.197],
        [-1.748, 0.908, 0.832], [-1.769, -0.861, 0.839], [-1.662, 0.017, -0.726],
        [0.298, 2.087, 0.377], [0.512, 1.272, -1.204], [1.795, 1.182, 0.045]]
TMA = (["C", "N", "C", "C"] + ["H"] * 9, np.array(_TMA))
TMA_HF2 = (TMA[0] + ["H", "F", "H", "F"],  # N ··· H13–F14 ··· H15–F16, hydrogen-bonded chain
           np.array(_TMA + [[0.466, 0, 2.06], [0.669, 0, 2.968], [-0.679, 0, 4.005],
                            [-1.416, 0, 4.572]]))
MALONALDEHYDE = (["O", "C", "C", "C", "O", "H", "H", "H", "H"],  # enol, O0–H5 ··· O4
                 np.array([[1.217, 0.791, -0.974], [0.949, -0.409, -0.520], [-0.127, -0.713, 0.249],
                           [-1.057, 0.318, 0.612], [-0.954, 1.490, 0.271], [0.496, 1.402, -0.649],
                           [1.674, -1.160, -0.815], [-0.292, -1.719, 0.590],
                           [-1.906, 0.000, 1.237]]))
NH3_HF = (["N", "H", "H", "H", "H", "F"],  # H3N ··· H–F, CREST seed of amine_pilot2 (F···H 3.0 Å)
          np.array([[-1.3137, 0, 0], [-1.6606, -0.1372, 0.9409], [-1.6644, 0.8821, -0.3518],
                    [-1.6609, -0.7473, -0.588], [0.2409, 0.0054, -0.0015], [1.2203, 0, 0]]))
NH3_2HF = (["N", "H", "H", "H", "H", "F", "H", "F"],  # two HF donate to N at 1.70 and 2.40 Å
           np.array([[0, 0, 0], [0.94, 0, -0.38], [-0.47, 0.814, -0.38], [-0.47, -0.814, -0.38],
                     [0, 0, 1.7], [0, 0, 2.62], [1.039, 1.8, 1.2], [1.438, 2.49, 1.66]]))


def _drives(system, charge=0, multiplicity=1):
    _, ts = trials.generate("m", *system, charge=charge, multiplicity=multiplicity,
                            max_trials=100)
    return [(t.kind, set(t.associations), set(t.dissociations)) for t in ts]


def test_backside_sn2_first_and_no_halogen_h_shift():
    drives = _drives(SN2, charge=-1)
    assert drives[0] == ("transfer", {(0, 5)}, {(0, 4)})  # Cl'–C–Cl closest to linear
    assert ("dissociation", set(), {(0, 4)}) in drives  # charged: heterolysis is a drive
    symbols = SN2[0]
    assert not any({symbols[i] for p in form for i in p} == {"Cl", "H"} for _, form, _ in drives)


def test_h_abstraction_by_oh_and_by_an_h_atom():
    assert ("transfer", {(1, 5)}, {(0, 1)}) in _drives(CH4_OH, multiplicity=2)
    assert ("transfer", {(1, 5)}, {(0, 1)}) in _drives(H_CH4, multiplicity=2)


def test_dissociation_only_for_open_shell_or_charged_sources():
    drives = _drives(CH3O, multiplicity=2)
    assert ("transfer", {(1, 3)}, {(0, 3)}) in drives  # 1,2-H shift C -> O
    cuts = [cut for kind, _, cut in drives if kind == "dissociation"]
    assert cuts == [{(0, 2)}, {(0, 1)}]  # one C–H class (most stretched first), then C–O
    assert ("transfer", {(0, 3)}, {(0, 1)}) in _drives(HONO)  # 1,3-H shift O -> O
    assert all(kind != "dissociation" for kind, _, _ in _drives(HONO))


def test_tma_hf2_starts_with_the_n_side_transfer_and_the_hf_relay():
    _, ts = trials.generate("m", *TMA_HF2)
    assert [(t.kind, set(t.associations), set(t.dissociations)) for t in ts[:2]] == [
        ("transfer", {(1, 13)}, {(13, 14)}),
        ("relay", {(1, 13), (14, 15)}, {(13, 14), (15, 16)})]
    assert len(ts) <= 10 and len({t.trial_id for t in ts}) == len(ts)


def test_through_space_transfer_comes_before_12_shifts():
    """Geminal pairs are close by the bond angle: the skeletal 1,2-shifts of the enol (r/Σr_cov
    about 1.6) rank after its O–H···O proton transfer (1.7)."""
    drives = _drives(MALONALDEHYDE)
    assert drives[0] == ("transfer", {(4, 5)}, {(0, 5)})
    assert ("transfer", {(2, 4)}, {(3, 4)}) in drives[1:]  # O4 1,2-shift C3 -> C2 stays, later


def test_cyclic_relay_is_the_nh3_hf_double_h_exchange():
    """The amine_pilot2 degenerate exchange (−1265i): H4 F -> N, an NH3 H N -> F. The returning
    H is 3.0 Å from F (beyond Σr_vdW); the first transfer holds N and F together."""
    relays = [(form, cut) for kind, form, cut in _drives(NH3_HF) if kind == "relay"]
    assert len(relays) == 1  # the three N–H are one class
    [(form, cut)] = relays
    [h] = {i for p in form for i in p} - {0, 4, 5}
    assert (form, cut) == ({(0, 4), (h, 5)}, {(0, h), (4, 5)})


def test_equivalent_drives_are_one_trial_and_distinct_sites_stay():
    drives = _drives(TMA)  # nine methyl 1,2-H shifts, sixty methyl relays and exchanges
    assert len(drives) == 7 and [k for k, _, _ in drives].count("transfer") == 1
    n_to_h = [form for kind, form, _ in _drives(NH3_2HF) if kind == "transfer"
              and {NH3_2HF[0][i] for p in form for i in p} == {"N", "H"}]
    assert n_to_h == [{(0, 4)}, {(0, 6)}]  # same atom classes, N···H 1.70 vs 2.40 Å


def test_linear_hcn_is_bent_and_gets_the_12_shift():
    start, ts = trials.generate("m", *HCN)
    assert _drives(HCN) == [("transfer", {(0, 2)}, {(0, 1)})] and all(t.perturbed for t in ts)
    assert external_basis(HCN[0], start).shape[1] == 6
    v1, v2 = start[0] - start[1], start[2] - start[1]
    bend = 180 - np.degrees(np.arccos(v1 @ v2 / np.linalg.norm(v1) / np.linalg.norm(v2)))
    assert 5 < bend < 15 and np.abs(start - HCN[1]).max() < 0.3
    assert np.array_equal(start, trials.generate("m", *HCN)[0])  # seeded


def test_product_verdict_window_and_ts_requirement():
    def verdict(mechanism, barrier, reaction, ts=True):
        return trials.product_verdict(mechanism, ts_validated=ts, barrier_kj=barrier,
                                      reaction_kj=reaction)

    assert verdict("nt2", 10, 5, ts=False) == "ts_not_validated"
    assert verdict("nt2", 151, 5) == verdict("afir", None, 101) == "out_of_window"
    assert verdict("afir", None, None) == "out_of_window" and verdict("afir", None, 99) is None

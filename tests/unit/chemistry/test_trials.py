"""Reaction trials (design §8.2 explore): the four bond-change templates, the valence rule,
drive equivalence classes, their independence of the atom numbering (G5-P4), the cap and the
linear perturbation."""

import numpy as np
import pytest

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
WATER_HF_HCL = (["O", "H", "H", "H", "F", "H", "Cl"],  # HF and HCl donate to O, both at 2 Å
                np.array([[0, 0, 0], [0.76, 0, -0.59], [-0.76, 0, -0.59], [0, 0, 2.0],
                          [0, 0, 2.92], [0, 2.0, 0], [0, 3.28, 0]]))
NH3_2HF = (["N", "H", "H", "H", "H", "F", "H", "F"],  # two HF donate to N at 1.70 and 2.40 Å
           np.array([[0, 0, 0], [0.94, 0, -0.38], [-0.47, 0.814, -0.38], [-0.47, -0.814, -0.38],
                     [0, 0, 1.7], [0, 0, 2.62], [1.039, 1.8, 1.2], [1.438, 2.49, 1.66]]))


def _drives(system, charge=0, multiplicity=1, max_trials=100, order=None):
    """(kind, formed, broken) of each trial in the source's own numbering; ``order`` numbers
    the atoms of the generated source (new index k is atom order[k])."""
    symbols, x = system
    order = np.arange(len(symbols)) if order is None else order
    _, ts = trials.generate("m", [symbols[i] for i in order], x[order], charge=charge,
                            multiplicity=multiplicity, max_trials=max_trials)
    return [(t.kind, *({tuple(sorted(int(order[i]) for i in p)) for p in bonds}
                       for bonds in (t.associations, t.dissociations))) for t in ts]


def _ids(system, max_trials=100, **kw):
    return [t.trial_id for t in trials.generate("m", *system, max_trials=max_trials, **kw)[1]]


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


@pytest.mark.parametrize(("system", "charge", "multiplicity", "mirror"), [
    (TMA_HF2, 0, 1, False), (MALONALDEHYDE, 0, 1, False), (NH3_2HF, 0, 1, False),
    (SN2, -1, 1, True), (CH3O, 0, 2, True), (WATER_HF_HCL, 0, 1, True)])
@pytest.mark.parametrize("max_trials", [10, 100])
def test_trials_depend_on_the_structure_not_on_the_atom_numbering(
        system, charge, multiplicity, mirror, max_trials):
    """G5-P4: under atom permutations the trial ids (class description and source) come in the
    same order and, carried back, the bonds are the same (up to a mirror image of the source:
    which of two mirror atoms a trial names is not a structural fact)."""
    symbols, x = system
    kw = {"charge": charge, "multiplicity": multiplicity, "max_trials": max_trials}
    ids, drives = _ids(system, **kw), _drives(system, **kw)
    for seed in range(5):
        order = np.random.default_rng(seed).permutation(len(symbols))
        assert _ids(([symbols[i] for i in order], x[order]), **kw) == ids
        assert mirror or _drives(system, order=order, **kw) == drives


def test_classes_tied_at_the_cap_are_all_in_or_all_out():
    """HF and HCl donate to O at exactly 2 Å, both linear: their two O–H transfers tie in value,
    so a cap between them takes neither (which comes first is only their class description)."""
    full = trials.generate("m", *WATER_HF_HCL, max_trials=100)[1]
    tied = {t.trial_id for t in full
            if t.kind == "transfer" and {(0, 3), (0, 5)} & set(t.associations)}
    assert len(tied) == 2 and len(full) == 5
    for cap in range(len(full) + 1):
        kept = _ids(WATER_HF_HCL, max_trials=cap)
        assert len(kept) <= cap and kept == [t.trial_id for t in full if t.trial_id in kept]
        assert len(tied & set(kept)) in (0, 2)
    assert [kind for kind, _, _ in _drives(WATER_HF_HCL, max_trials=2)] == ["relay"]


def test_linear_hcn_is_bent_and_gets_the_12_shift():
    start, _ = trials.generate("m", *HCN)
    assert _drives(HCN) == [("transfer", {(0, 2)}, {(0, 1)})]
    assert external_basis(HCN[0], start).shape[1] == 6
    v1, v2 = start[0] - start[1], start[2] - start[1]
    bend = 180 - np.degrees(np.arccos(v1 @ v2 / np.linalg.norm(v1) / np.linalg.norm(v2)))
    assert 5 < bend < 15 and np.abs(start - HCN[1]).max() < 0.3
    assert np.array_equal(start, trials.generate("m", *HCN)[0])  # seeded

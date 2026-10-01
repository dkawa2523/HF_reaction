"""chemistry.thermo: textbook values and the fixed rule for negative modes."""

import importlib.metadata
import math

import numpy as np
import pytest

from hfauto.chemistry import thermo as th
from hfauto.chemistry.symmetry import Symmetry
from hfauto.chemistry.vibrations import external_basis
from hfauto.core.constants import HARTREE_TO_KCAL_MOL as H2K
from hfauto.core.constants import R_KCAL_MOL_K

RT = R_KCAL_MOL_K * 298.15


def test_standard_states_and_association():
    assert th.standard_state_shift(1, 298.15, "1M") == pytest.approx(1.894, abs=1e-3)
    bound, ts = -3.0 - 10 / H2K, -3.0 + 5 / H2K  # 10 kcal/mol below the separated monomers
    assert th.association(bound, ts, [-1.0, -2.0], [1, 1], 298.15, "1atm") == pytest.approx(
        (-10.0, 5.0))
    assert th.association(bound, None, [-1.5], [2], 298.15, "1M")[0] == pytest.approx(
        -10.0 - 1.894, abs=1e-3)  # dn = 1 - 2


@pytest.mark.parametrize(("G_ts", "G_P", "expected"), [
    (12.0, 3.0, 12.0),  # normal: the TS is the highest point
    (-0.5, 3.0, 3.0),  # dG_act < 0 (a submerged TS): the product state is the bottleneck
    (2.0, 5.0, 5.0),  # dG_act < dG_rxn: the TS lies below the product
    (None, 4.0, 4.0),  # barrierless (or submerged): max(dG_rxn, 0)
    (None, -6.0, 0.0),
])
def test_effective_barrier_is_the_highest_point_above_the_reactant(G_ts, G_P, expected):
    assert th.effective_barrier(G_ts, 0.0, G_P) == pytest.approx(expected)
    assert th.effective_barrier(None if G_ts is None else G_ts - 7.0, -7.0, G_P - 7.0
                                ) == pytest.approx(expected)  # relative to G_R


def test_a_chiral_structure_gains_minus_rt_ln2():  # m = 2 (HONO TS: 12.23 -> 11.82)
    assert th.chiral_G(298.15) * H2K == pytest.approx(-0.4107, abs=1e-4)
    assert th.chiral_G(298.15) == pytest.approx(-RT * math.log(2) / H2K)


def test_minimum_keeps_every_negative_mode_as_its_magnitude():  # noise (-4) and soft (-30)
    assert th.thermo_frequencies((500.0, -4.0, -30.0), saddle=False) == (4.0, 30.0, 500.0)


def test_saddle_drops_the_lowest_mode_and_flips_the_others():  # incl. an unannotated -8
    assert th.thermo_frequencies((600.0, -8.0, -700.0, -30.0), saddle=True) == (
        8.0, 30.0, 600.0)


def test_positive_or_empty_modes_pass_unchanged():
    assert th.thermo_frequencies((900.0, 300.0), saddle=False) == (300.0, 900.0)
    assert th.thermo_frequencies((), saddle=False) == th.thermo_frequencies((), saddle=True) == ()


def test_thermal_modes_follow_the_point_group_not_the_rank():
    """S18: H...H2 0.006 Å off the axis has rank 6 (3 modes, 2.27 kcal/mol in G); as C∞v at its
    symmetrized structure it has 3N - 5 = 4 modes, a saddle 3."""
    x = np.array([[0.0, 0.0, -1.8], [0.0, 0.006, 0.0], [0.0, 0.0, 0.74]])
    line = x * [0.0, 0.0, 1.0]
    h = 0.3 * np.eye(9)
    assert external_basis(["H"] * 3, x).shape[1] == 6
    linear, bent = (Symmetry(g, 1, lin, 1, c) for g, lin, c in (("Cinfv", True, line),
                                                                  ("C1", False, x)))
    assert len(th.thermal_modes(["H"] * 3, linear, h, saddle=False)) == 4
    assert len(th.thermal_modes(["H"] * 3, linear, h, saddle=True)) == 3
    assert len(th.thermal_modes(["H"] * 3, bent, h, saddle=False)) == 3
    assert th.thermal_modes(["H"], Symmetry("Kh", 1, False, 1, np.zeros((1, 3))),
                            np.zeros((3, 3)), saddle=False) == ()


def test_pymsym_is_the_pinned_version():  # R14: another libmsym build may find another sigma
    pytest.importorskip("pymsym")
    assert importlib.metadata.version("pymsym") == "0.3.5"

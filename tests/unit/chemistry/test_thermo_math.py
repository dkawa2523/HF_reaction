"""chemistry.thermo: textbook values and the fixed rule for negative modes."""

import math

import pytest

from hfauto.chemistry import thermo as th
from hfauto.core.constants import HARTREE_TO_KCAL_MOL as H2K
from hfauto.core.constants import R_KCAL_MOL_K

RT = R_KCAL_MOL_K * 298.15


def test_standard_states_populations_and_association():  # + port of test_basin_populations
    assert th.standard_state_shift(1, 298.15, "1M") == pytest.approx(1.894, abs=1e-3)
    weights, w = th.boltzmann_populations([-10.0, -10.0 + 1 / H2K], 298.15), math.exp(-1 / RT)
    assert sum(weights) == pytest.approx(1.0) and weights[1] == pytest.approx(w / (1 + w))
    assert th.ensemble_G([-10.0, -10.0], 298.15) == pytest.approx(-10.0 - RT * math.log(2) / H2K)
    bound, ts = -3.0 - 10 / H2K, -3.0 + 5 / H2K  # 10 kcal/mol below the separated monomers
    assert th.association(bound, ts, [-1.0, -2.0], [1, 1], 298.15, "1atm") == pytest.approx(
        (-10.0, 5.0))
    assert th.association(bound, None, [-1.5], [2], 298.15, "1M")[0] == pytest.approx(
        -10.0 - 1.894, abs=1e-3)  # dn = 1 - 2
    assert th.composite(-2.0, -0.95, -1.0) == pytest.approx(-1.95)


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

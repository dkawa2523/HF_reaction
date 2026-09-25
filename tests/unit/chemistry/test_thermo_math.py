"""chemistry.thermo: textbook values and the scale-factor lookup order."""

import math

import pytest

from hfauto.chemistry import thermo as th
from hfauto.core.constants import HARTREE_TO_KCAL_MOL as H2K
from hfauto.core.constants import R_KCAL_MOL_K
from hfauto.core.evidence import Level
from hfauto.core.method import ThermoSettings

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


@pytest.mark.parametrize(("method", "basis", "expected"), [
    ("wb97x-d", "def2-TZVP", (0.989, 0.975, None)),
    ("pbe0", "def2-svpd", (0.989, 0.975, "scale_factor_from:PBE0/MG3S")),
    ("wb97x-d3", "def2-tzvpd", (0.989, 0.975, "scale_factor_from:wB97XD/def2TZVP")),
    ("gfn2", None, (1.0, 1.0, "scale_factor_unverified")),
])
def test_scale_factor_lookup_order(method, basis, expected):
    level = Level(program="p", version="1", method=method, basis=basis, charge=0, multiplicity=1)
    assert th.scale_factors(level) == expected
    assert th.resolve_scales(ThermoSettings(vib_scale=0.97), level)[:2] == (0.97, expected[1])

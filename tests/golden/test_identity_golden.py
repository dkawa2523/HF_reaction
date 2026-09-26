"""G21: two water optimizations land in the same minimum (design §10.2)."""

import re

import pytest

from hfauto.chemistry.identity import permutation_invariant_rmsd, same_minimum
from hfauto.chemistry.xyz import read_xyz

pytestmark = pytest.mark.golden


def test_g21_water_structures_are_same(golden):
    ref = read_xyz(golden.path("nwchem/G21/water_reference_final.xyz"))
    other = read_xyz(golden.path("nwchem/G21/water_distorted_final.xyz"))
    outputs = [golden.text(f"nwchem/G21/water_{n}.out") for n in ("reference", "distorted")]
    energies = [float(re.findall(r"Total DFT energy =\s+(\S+)", out)[-1]) for out in outputs]
    rmsd, _ = permutation_invariant_rmsd(ref.symbols, ref.coords, other.coords)
    assert rmsd == pytest.approx(8e-6, abs=2e-6)
    assert same_minimum(ref.symbols, ref.coords, other.coords, *energies)

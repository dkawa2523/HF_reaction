"""G13 / G14 / G15: string convergence, stagnation and initial paths (design §10.2)."""

import re

import numpy as np
import pytest

from hfauto.chemistry.interpolation import idpp, min_interatomic_distance
from hfauto.chemistry.profile import stagnated, string_converged
from hfauto.chemistry.xyz_trajectory import read_xyz_trajectory

pytestmark = pytest.mark.golden


def test_g13_worsening_gmax_is_not_converged(golden):
    text = golden.text("nwchem/G13/nwchem_string.out")
    gmax = [float(v) for v in re.findall(r"string: gmax,grms,xrms,xmax=\s+(\S+)", text)]
    assert gmax[0] < 1e-3 < gmax[-1]
    assert not string_converged(gmax)


def test_g14_zero_steps_are_stagnation(golden):
    text = golden.text("nwchem/G14/hfauto_neb.neb_epath")
    gmax = [float(v) for v in re.findall(r"# Gmax\s+=\s+(\S+)", text)]
    xmax = [float(v) for v in re.findall(r"# Xmax\s+=\s+(\S+)", text)]
    assert stagnated(gmax, xmax)
    assert stagnated((), xmax)


def test_g15_collided_path_versus_idpp(golden):
    frames = read_xyz_trajectory(golden.path("nwchem/G15/hfauto_neb.neb_final.xyz"))
    # The fixture is the final NEB path (0.6015 Å); the review's 0.529 Å is its initial path.
    shortest = min(min_interatomic_distance(f.coords) for f in frames)
    assert shortest == pytest.approx(0.6015, abs=1e-3)
    path = idpp(frames[0].symbols, frames[0].coords, frames[-1].coords, 11)
    assert len(path) == 11 and np.array_equal(path[0], frames[0].coords)
    assert min(min_interatomic_distance(c) for c in path) >= 0.7

"""G13 / G15: settled string energies and initial paths (design §10.2; G14 is retired)."""

import re

import numpy as np
import pytest

from hfauto.chemistry.interpolation import idpp, min_interatomic_distance
from hfauto.chemistry.profile import energies_settled
from hfauto.chemistry.xyz_trajectory import read_xyz_trajectory
from hfauto.core.constants import HARTREE_TO_KCAL_MOL

pytestmark = pytest.mark.golden


def test_g13_bead_energies_settle_while_gmax_worsens(golden):
    text = golden.text("nwchem/G13/nwchem_string.out")
    gmax = [float(v) for v in re.findall(r"string: gmax,grms,xrms,xmax=\s+(\S+)", text)]
    blocks = re.findall(r"Path Energy #.*\n((?: string:[ \t]+\d+[ \t]+\S+[ \t]*\n)+)", text)
    history = [[float(row.split()[2]) for row in block.splitlines()] for block in blocks]
    assert gmax[0] < 1e-3 < gmax[-1] and len(history) == 3
    assert energies_settled(history, 0.1 / HARTREE_TO_KCAL_MOL)  # largest step 0.095 kcal/mol
    assert not energies_settled(history, 0.05 / HARTREE_TO_KCAL_MOL)


def test_g15_collided_path_versus_idpp(golden):
    frames = read_xyz_trajectory(golden.path("nwchem/G15/hfauto_neb.neb_final.xyz"))
    # The fixture is the final NEB path (0.6015 Å); the review's 0.529 Å is its initial path.
    shortest = min(min_interatomic_distance(f.coords) for f in frames)
    assert shortest == pytest.approx(0.6015, abs=1e-3)
    path = idpp(frames[0].symbols, frames[0].coords, frames[-1].coords, 11)
    assert len(path) == 11 and np.array_equal(path[0], frames[0].coords)
    assert min(min_interatomic_distance(c) for c in path) >= 0.7

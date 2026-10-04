"""G13 / G15: a string classified as it stands, and initial paths (design §10.2; G14 is
retired)."""

import re

import numpy as np
import pytest

from hfauto.chemistry.gates import Policy
from hfauto.chemistry.interpolation import idpp, min_interatomic_distance
from hfauto.chemistry.profile import classify
from hfauto.chemistry.xyz import read_xyz_trajectory
from hfauto.core.constants import HARTREE_TO_KCAL_MOL

pytestmark = pytest.mark.golden


def g13_final_beads(golden) -> list[float]:
    """Bead energies of the last iteration of the TMA·(HF)2 proton reorganization string."""
    text = golden.text("nwchem/G13/nwchem_string.out")
    blocks = re.findall(r"Path Energy #.*\n((?: string:[ \t]+\d+[ \t]+\S+[ \t]*\n)+)", text)
    return [float(row.split()[2]) for row in blocks[-1].splitlines()]


def test_g13_string_is_barrierless_at_the_resolution(golden):
    """U5-P1 / U5-P4: the chunk is classified as it stands (its gmax rose 9.1e-4 -> 9.2e-3):
    the beads fall 0.93 kcal/mol with a 0.01 kcal/mol ripple, no hill and no well."""
    resolution = Policy().resolution_kcal / HARTREE_TO_KCAL_MOL
    assert classify(g13_final_beads(golden), resolution) == "barrierless"


def test_g15_collided_path_versus_idpp(golden):
    frames = read_xyz_trajectory(golden.path("nwchem/G15/hfauto_neb.neb_final.xyz"))
    # The fixture is the final NEB path (0.6015 Å); the review's 0.529 Å is its initial path.
    shortest = min(min_interatomic_distance(f.coords) for f in frames)
    assert shortest == pytest.approx(0.6015, abs=1e-3)
    path = idpp(frames[0].symbols, frames[0].coords, frames[-1].coords, 11)
    assert len(path) == 11 and np.array_equal(path[0], frames[0].coords)
    assert min(min_interatomic_distance(c) for c in path) >= 0.7

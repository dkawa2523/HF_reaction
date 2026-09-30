"""Mapped alignment, resampling and IDPP (design §5.5, CH-09)."""

from itertools import pairwise

import numpy as np
import pytest

from hfauto.chemistry import interpolation as itp
from hfauto.chemistry.geometry import declared_coordinate
from hfauto.core.records import CoordinateTerm

SYMBOLS = ["N", "H", "H", "H"]
NH3 = np.array([[0.0, 0.0, 0.38], [0.94, 0.0, 0.0], [-0.47, 0.814, 0.0], [-0.47, -0.814, 0.0]])
# PBE0-D3BJ/def2-SVPD minima of HONO (H O N O), flattened to z = 0
CIS_HONO = np.array([[-0.32417822, -0.09440887, 0.0], [0.39510676, -0.75997798, 0.0],
                     [1.56130940, -0.07839078, 0.0], [1.42966204, 1.09297763, 0.0]])
TRANS_HONO = np.array([[-0.19774781, -0.94387557, 0.0], [0.02491128, -0.00145366, 0.0],
                       [1.40706589, 0.00832925, 0.0], [1.82767063, 1.09719998, 0.0]])


def largest_step(frames):
    return max(float(np.linalg.norm(b - a)) for a, b in pairwise(frames))


def test_align_mapped_undoes_rigid_motion():
    c, s = np.cos(0.7), np.sin(0.7)
    moved = NH3 @ np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]).T + [3.0, 1.0, -2.0]
    assert itp.align_mapped(NH3, moved) == pytest.approx(NH3, abs=1e-10)


def test_idpp_inversion_path_is_smooth_with_fixed_ends():
    inverted = NH3 * [1.0, 1.0, -1.0]
    path = itp.idpp(SYMBOLS, NH3, inverted, 11)
    assert len(path) == 11
    assert np.array_equal(path[0], NH3)
    assert np.array_equal(path[-1], itp.align_mapped(NH3, inverted))
    assert largest_step(path) < 0.1
    assert min(itp.min_interatomic_distance(x) for x in path) > 0.9


def test_resample_spaces_points_evenly_and_keeps_the_ends():
    frames = [NH3, NH3 + [0.1, 0.0, 0.0], NH3 + [0.4, 0.0, 0.0]]  # arc lengths 0.2 and 0.6
    out = itp.resample(frames, 5)
    steps = [np.linalg.norm(b - a) for a, b in pairwise(out)]
    assert len(out) == 5 and np.allclose(steps, steps[0])
    assert np.array_equal(out[0], NH3) and np.array_equal(out[-1], frames[-1])
    with pytest.raises(ValueError, match="zero-length"):
        itp.resample([NH3, NH3], 5)


def test_idpp_rejects_short_contacts():
    squeezed = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
    with pytest.raises(ValueError, match="0.7"):
        itp.idpp(["H", "H"], squeezed, squeezed + [0.0, 0.1, 0.0], 5)


def test_align_sequential_removes_a_rigid_jump_inside_a_path():
    inverted = NH3 * [1.0, 1.0, -1.0]
    path = [NH3 + t * (inverted - NH3) for t in np.linspace(0.0, 1.0, 7)]
    path = [x - x.mean(axis=0) for x in path]  # each frame already aligned onto its neighbours
    rotated = list(path)
    rotated[3] = path[3] @ np.array([[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]]).T
    assert largest_step(rotated) > 1.0
    aligned = itp.align_sequential(rotated)
    assert np.array_equal(aligned[0], rotated[0])
    assert largest_step(aligned) == pytest.approx(largest_step(path), abs=1e-10)


def test_idpp_leaves_the_plane_between_planar_cis_trans_hono():
    # S18: a converged IDPP twists out of the plane (an unconverged one inverts H-O-N in-plane)
    path = itp.idpp(["H", "O", "N", "O"], CIS_HONO, TRANS_HONO, 11)
    hon = CoordinateTerm(kind="angle", atoms=(0, 1, 2))
    hono = CoordinateTerm(kind="dihedral", atoms=(0, 1, 2, 3))
    angle, dihedral = ([abs(declared_coordinate([t], x)) for x in path] for t in (hon, hono))
    assert max(angle[1:-1]) < 120.0
    assert dihedral[0] < 1.0 and dihedral[-1] > 179.0 and np.all(np.diff(dihedral) > 0.0)
    assert any(30.0 < d < 150.0 for d in dihedral[1:-1])

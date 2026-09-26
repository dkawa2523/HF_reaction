"""Mapped alignment and IDPP (design §5.5, CH-09)."""

import numpy as np
import pytest

from hfauto.chemistry import interpolation as itp
from hfauto.chemistry.profile import max_node_spacing

SYMBOLS = ["N", "H", "H", "H"]
NH3 = np.array([[0.0, 0.0, 0.38], [0.94, 0.0, 0.0], [-0.47, 0.814, 0.0], [-0.47, -0.814, 0.0]])


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
    assert max_node_spacing(path) < 0.1
    assert min(itp.min_interatomic_distance(x) for x in path) > 0.9


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
    assert max_node_spacing(rotated) > 1.0
    aligned = itp.align_sequential(rotated)
    assert np.array_equal(aligned[0], rotated[0])
    assert max_node_spacing(aligned) == pytest.approx(max_node_spacing(path), abs=1e-10)

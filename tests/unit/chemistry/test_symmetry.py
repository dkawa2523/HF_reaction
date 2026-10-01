"""chemistry.symmetry: one point group, accepted by the basin criterion (analysis 2026-09-30
X3-2). Model Hessians: k·I (Eh/bohr²), so a symmetrization rises ½ k |δx|²."""

import math

import numpy as np
import pytest

from hfauto.chemistry.identity import BASIN_DE_HARTREE
from hfauto.chemistry.symmetry import analyze
from hfauto.core.constants import BOHR_TO_ANGSTROM

pytest.importorskip("pymsym")
NH3_SYMBOLS = ["N", "H", "H", "H"]
NH3 = np.array([[0.0, 0.0, 0.38], [0.94, 0.0, 0.0], [-0.47, 0.814, 0.0], [-0.47, -0.814, 0.0]])
CH4 = np.array([[0, 0, 0], [0.63, 0.63, 0.63], [-0.63, -0.63, 0.63], [-0.63, 0.63, -0.63],
                [0.63, -0.63, -0.63]])
H2O2 = np.array([[0, 0.7, 0], [0, -0.7, 0], [0.9, 1, 0.3], [-0.9, -1, 0.3]])  # gauche, C2
# I-...CH3I as NWChem left it in a validation run, the I- 0.03° off the C3 axis
SN2_I_SYMBOLS = ["C", "H", "H", "H", "I", "I"]
SN2_I = np.array([[-0.00040698, -1.12e-06, 0.04484225], [1.04164864, -8.4e-07, 0.37030592],
                  [-0.5205438, 0.90190957, 0.37320459], [-0.52054263, -0.90191372, 0.37320135],
                  [0.00429119, 3.12e-06, 3.5046527], [-0.00444641, 3e-06, -2.1462068]])


def _iso(x: np.ndarray, k: float = 0.5) -> np.ndarray:
    return k * np.eye(np.size(x))


def _key(symbols, x, hessian=None):
    s = analyze(symbols, x, _iso(x) if hessian is None else hessian)
    return s.point_group, s.sigma, s.linear, s.m


def _hcn(angle_deg: float) -> np.ndarray:
    t = math.radians(180.0 - angle_deg)
    h = [-1.066 * math.cos(t), 1.066 * math.sin(t), 0.0]
    return np.array([[0.0, 0.0, 0.0], [1.153, 0.0, 0.0], h])


@pytest.mark.parametrize(("symbols", "x", "expected"), [
    (NH3_SYMBOLS, NH3, ("C3v", 3, False, 1)),
    (["C", "H", "H", "H", "H"], CH4, ("Td", 12, False, 1)),
    (["O", "O", "H", "H"], H2O2, ("C2", 2, False, 2)),  # chiral: no improper operation
    (["C", "H", "F", "Cl", "Br"], CH4, ("C1", 1, False, 2)),
    (SN2_I_SYMBOLS, SN2_I, ("C3v", 3, False, 1)),
    (["C", "N", "H"], _hcn(180.0), ("Cinfv", 1, True, 1)),
    (["O", "O"], np.array([[0, 0, 0], [0, 0, 1.2]]), ("Dinfh", 2, True, 1)),
    (["O", "C", "O"], np.array([[-1.160, 0, 0], [0, 0, 0], [1.161, 0, 0]]), ("Dinfh", 2, True, 1)),
    (["Ar"], np.zeros((1, 3)), ("Kh", 1, False, 1)),
])
def test_known_point_groups(symbols, x, expected):
    assert _key(symbols, x) == expected


def test_libmsym_default_splits_the_complex_the_accepted_group_does_not():
    import pymsym

    assert pymsym.get_symmetry_number([6, 1, 1, 1, 53, 53], SN2_I.tolist()) == 1  # Cs
    sym = analyze(SN2_I_SYMBOLS, SN2_I, _iso(SN2_I))
    assert sym.point_group == "C3v" and sym.coords.shape == SN2_I.shape


def test_a_symmetrized_linear_structure_is_exactly_linear_in_the_input_frame():
    """S18: H...H2 0.006 Å off the axis is one basin with its C∞v projection."""
    x = np.array([[0.0, 0.0, -1.8], [0.0, 0.006, 0.0], [0.0, 0.0, 0.74]])
    sym = analyze(["H", "H", "H"], x, _iso(x))
    centred = sym.coords - sym.coords.mean(axis=0)
    assert (sym.point_group, sym.linear) == ("Cinfv", True)
    assert np.linalg.matrix_rank(centred, tol=1e-10) == 1
    assert np.sqrt(np.mean(np.sum((sym.coords - x) ** 2, axis=1))) < 0.005  # not re-oriented


def test_a_bent_molecule_is_not_projected_onto_a_line():
    assert _key(["C", "N", "H"], _hcn(170.0)) == ("Cs", 1, False, 1)


def test_a_group_counts_only_within_the_basin_energy():
    """NH3 with one H 0.02 Å out of its C3v place, within a mirror plane: the C3v structure is
    ½ k |δx|² away, inside the basin on a soft surface and outside it on a stiff one."""
    x = NH3.copy()
    x[1] += [0.02, 0.0, 0.0]
    soft, stiff = _key(NH3_SYMBOLS, x, _iso(x, 0.01)), _key(NH3_SYMBOLS, x, _iso(x, 0.5))
    assert soft == ("C3v", 3, False, 1) and stiff == ("Cs", 1, False, 1)
    sym = analyze(NH3_SYMBOLS, x, _iso(x, 0.01))
    d = (sym.coords - x).ravel() / BOHR_TO_ANGSTROM
    assert 0.5 * 0.5 * d @ d > BASIN_DE_HARTREE >= 0.5 * 0.01 * d @ d  # the stiff rise is too high


def test_a_group_counts_only_within_the_basin_rmsd():
    """A flat surface (H = 0) accepts any rise; the RMSD to the group's structure still has to
    stay within BASIN_A."""
    x = NH3.copy()
    x[1] += [0.15, 0.0, 0.0]  # C3v lies ~0.06 Å RMSD away, Cs holds
    assert _key(NH3_SYMBOLS, x, np.zeros((12, 12))) == ("Cs", 1, False, 1)
    x[1] -= [0.1, 0.0, 0.0]
    assert _key(NH3_SYMBOLS, x, np.zeros((12, 12))) == ("C3v", 3, False, 1)


def test_a_saddle_counts_its_negative_curvature_by_magnitude():
    x = NH3.copy()
    x[1] += [0.02, 0.0, 0.0]
    assert _key(NH3_SYMBOLS, x, -_iso(x, 0.5)) == _key(NH3_SYMBOLS, x, _iso(x, 0.5))


@pytest.mark.parametrize(("symbols", "x"), [
    (NH3_SYMBOLS, NH3), (["C", "H", "H", "H", "H"], CH4), (["O", "O", "H", "H"], H2O2),
    (SN2_I_SYMBOLS, SN2_I), (["C", "N", "H"], _hcn(180.0))])
def test_the_group_survives_noise(symbols, x):
    rng, ref = np.random.default_rng(5), _key(symbols, x)
    for _ in range(20):
        noisy = x + rng.normal(size=x.shape) * 1e-3 / np.sqrt(3)
        assert _key(symbols, noisy, _iso(x)) == ref


def test_the_hessian_must_match_the_structure():
    with pytest.raises(ValueError, match="Hessian"):
        analyze(NH3_SYMBOLS, NH3, np.eye(9))

"""Permutation-invariant identity: mirror images are one basin, labels stay proper (§5.5)."""

import numpy as np
import pytest

from hfauto.chemistry import identity as idn
from hfauto.chemistry import topology

NH3_SYMBOLS = ["N", "H", "H", "H"]
NH3 = np.array([[0.0, 0.0, 0.38], [0.94, 0.0, 0.0], [-0.47, 0.814, 0.0], [-0.47, -0.814, 0.0]])
INVERTED = NH3 * [1.0, 1.0, -1.0]
CHFCLBR_SYMBOLS = ["C", "H", "F", "Cl", "Br"]
CHFCLBR = np.array(
    [[0.0, 0.0, 0.0], [0.63, 0.63, 0.63], [-0.8, -0.8, 0.8], [-1.0, 1.0, -1.0], [1.1, -1.1, -1.1]]
)
MIRROR = CHFCLBR * [-1.0, 1.0, 1.0]
CH4 = np.array([[0, 0, 0], [0.63, 0.63, 0.63], [-0.63, -0.63, 0.63], [-0.63, 0.63, -0.63],
                [0.63, -0.63, -0.63]])
# QRC ± starts displaced from the DFT TS (VAL runs, Å): images at the D3h NH3 and SN2 TSs,
# not at the HCN -> HNC TS.
QRC_STARTS = {  # name: (symbols, plus, minus)
    "nh3": (NH3_SYMBOLS,
            np.array([[0.0, 0.0, 0.01079], [1.00157, -2e-05, -0.05019],
                      [-0.5008, -0.86742, -0.04958], [-0.5008, 0.86738, -0.05023]]),
            np.array([[0.0, 0.0, -0.0108], [1.00159, 2e-05, 0.04981],
                      [-0.50078, -0.86738, 0.05042], [-0.50078, 0.86741, 0.04977]])),
    "sn2_cl": (["C", "H", "H", "H", "Cl", "Cl"],
               np.array([[1e-05, 1e-05, 0.05], [-0.54054, 0.93624, 0.00099],
                         [1.08107, 0.0, 0.0005], [-0.54054, -0.93624, 0.00115],
                         [-0.00083, -0.0002, -2.32979], [0.00084, 0.00019, 2.31256]]),
               np.array([[-3e-05, 0.0, -0.05], [-0.54054, 0.93624, -0.00076],
                         [1.08107, 0.0, -0.00127], [-0.54054, -0.93624, -0.00061],
                         [-0.00082, -0.00019, -2.31256], [0.00085, 0.00019, 2.32979]])),
    "hcn": (["C", "N", "H"],
            np.array([[-0.02629, -0.46541, -0.42539], [0.72056, 0.21233, 0.19408],
                      [-0.73139, 0.22832, 0.20868]]),
            np.array([[-0.02406, -0.47242, -0.4318], [0.71292, 0.21442, 0.19598],
                      [-0.65175, 0.28276, 0.25844]])),
}


def _proper_rotation(seed: int) -> np.ndarray:
    q, _ = np.linalg.qr(np.random.default_rng(seed).standard_normal((3, 3)))
    return q * np.sign(np.linalg.det(q))


def test_ammonia_inversion_is_one_basin_and_a_degenerate_pair():
    assert idn.assign(NH3_SYMBOLS, NH3, -56.5, {"inverted": (INVERTED, -56.5)}) == "inverted"
    assert idn.mapped_equivalent(NH3_SYMBOLS, NH3, INVERTED)  # labels differ: proper only
    assert not idn.same_as_labelled(NH3, INVERTED, -56.5, -56.5)  # two structures as labelled
    assert not idn.mapped_equivalent(NH3_SYMBOLS, NH3, NH3)
    assert not idn.is_chiral(NH3_SYMBOLS, NH3)


def test_same_as_labelled_is_the_identity_mapping_within_the_basin_tolerances():
    moved = NH3 @ _proper_rotation(4).T + [0.5, 0.0, -1.0]
    assert idn.same_as_labelled(NH3, moved, 0.0, 3e-5)
    assert not idn.same_as_labelled(NH3, moved, 0.0, 6e-5)  # 5e-5 Eh
    assert not idn.same_as_labelled(NH3, NH3[[0, 2, 1, 3]] @ _proper_rotation(4).T, 0.0, 0.0)


def test_enantiomers_are_one_chiral_basin():
    assert idn.permutation_invariant_rmsd(CHFCLBR_SYMBOLS, CHFCLBR, MIRROR)[0] > 1.0  # proper
    moved = MIRROR @ _proper_rotation(5).T + [0.3, 0.0, -1.0]
    assert idn.assign(CHFCLBR_SYMBOLS, CHFCLBR, -1.0, {"s": (moved, -1.0)}) == "s"
    assert idn.is_chiral(CHFCLBR_SYMBOLS, CHFCLBR) and idn.is_chiral(CHFCLBR_SYMBOLS, moved)
    assert idn.mapped_equivalent(CHFCLBR_SYMBOLS, CHFCLBR, MIRROR)  # R -> S is degenerate


def test_methane_is_achiral():
    assert not idn.is_chiral(["C", "H", "H", "H", "H"], CH4)


def test_basin_coords_take_the_atom_order_and_handedness_of_the_member():
    member = MIRROR @ _proper_rotation(7).T + 0.02  # the S member of an R basin
    x = idn.basin_coords(CHFCLBR_SYMBOLS, CHFCLBR, member)
    assert idn.mapped_rmsd(x, member) < 0.05 < idn.mapped_rmsd(x, CHFCLBR)
    order = [0, 2, 3, 1]  # an achiral basin is only relabelled
    x = idn.basin_coords(NH3_SYMBOLS, NH3, NH3[order] @ _proper_rotation(3).T + 0.03)
    assert np.allclose(x, NH3[order])


def test_member_coords_return_the_representatives_input_bit_for_bit():
    """G7-P4: the representative takes the basin's coordinates unchanged (an alignment could
    move the job keys downstream); any other member, basin_coords."""
    basin = NH3 @ _proper_rotation(2).T + 0.1234567891
    assert idn.member_coords(NH3_SYMBOLS, basin, "rep", "rep", INVERTED).tobytes() == (
        basin.tobytes())
    member = NH3[[0, 2, 3, 1]] @ _proper_rotation(3).T + 0.03
    assert np.array_equal(idn.member_coords(NH3_SYMBOLS, NH3, "rep", "m", member),
                          idn.basin_coords(NH3_SYMBOLS, NH3, member))


def test_basin_coords_keep_the_bonds_of_a_member_far_from_the_basin():
    """S6 (W5, VAL7): an xTB product CH3OH + H, 0.96 A from the PBE0 basin that holds it. As
    elements alone the nearest match moves H1 onto O and H6 onto C (a fake H scramble); within
    atom classes the member's atom-indexed bonds are kept."""
    symbols = ["C", "H", "H", "H", "H", "O", "H"]
    basin = np.array([[-0.3753, -0.5165, 0.3011], [-0.8441, -1.5066, 0.3711],
                      [0.6085, 0.4446, -1.0785], [0.3826, -0.4401, 1.1004],
                      [-1.1552, 0.2442, 0.4804], [0.188, -0.413, -0.982],
                      [1.0871, 2.5978, 0.1369]])
    own = np.array([[-0.0454, 0.2054, -0.3942], [-0.6118, -0.64, 0.1764],
                    [-0.0044, -0.0841, -1.4423], [-0.6312, 1.1203, -0.2668],
                    [-1.2368, -1.598, 0.9061], [1.2626, 0.3256, 0.0606],
                    [1.267, 0.6708, 0.9603]])
    by_element = basin[idn.permutation_invariant_rmsd(symbols, own, basin)[1]]
    assert topology.bonds(symbols, by_element) != topology.bonds(symbols, own)
    x = idn.basin_coords(symbols, basin, own)
    assert topology.bonds(symbols, x) == topology.bonds(symbols, own)


def test_rotated_permuted_copy_is_recovered():
    order = [0, 3, 1, 2]
    moved = (NH3 @ _proper_rotation(3).T + [1.0, -2.0, 0.5])[order]
    rmsd, perm = idn.permutation_invariant_rmsd(NH3_SYMBOLS, NH3, moved)
    assert rmsd < 1e-6
    assert [order[i] for i in perm] == [0, 1, 2, 3]


def test_one_criterion_energy_first_then_the_best_structure():
    squeezed = NH3 * [1.0, 1.0, 0.5]
    assert idn.assign(NH3_SYMBOLS, NH3, 3e-5, {"up": (NH3, 0.0)}) == "up"  # 5e-5 Eh
    assert idn.assign(NH3_SYMBOLS, NH3, 6e-5, {"up": (NH3, 0.0)}) is None
    assert idn.assign(NH3_SYMBOLS, squeezed, 0.0, {"up": (NH3, 0.0)}) is None  # 0.05 A
    candidates = {"up": (NH3, 0.0), "flat": (squeezed, 0.0)}
    assert idn.assign(NH3_SYMBOLS, NH3 + 1e-4, 1e-6, candidates) == "up"
    near = NH3 * [1.0, 1.0, 0.98]  # both within 0.05 A: the closer one joins (U3-P7)
    assert idn.assign(NH3_SYMBOLS, NH3, 0.0, {"a": (near, 0.0), "b": (NH3, 0.0)}) == "b"
    twins = {"b": (NH3, 0.0), "a": (NH3 + 0.004, 0.0)}  # one structure: the lowest id
    assert idn.assign(NH3_SYMBOLS, NH3 + 0.002, 0.0, twins) == "a"
    twins["a"] = (NH3 + 0.004, 1e-3)  # another energy never competes
    assert idn.assign(NH3_SYMBOLS, NH3 + 0.002, 0.0, twins) == "b"


@pytest.mark.parametrize("name", ["nh3", "sn2_cl"])
def test_the_qrc_starts_of_a_symmetric_ts_are_exact_images(name):
    symbols, plus, minus = QRC_STARTS[name]
    rmsd, carried = idn.carry(symbols, minus, plus, plus)
    assert rmsd < idn.IMAGE_A and idn.is_image(symbols, minus, plus)
    assert np.sqrt(np.mean(np.sum((carried - minus) ** 2, axis=1))) == pytest.approx(rmsd)


def test_the_qrc_starts_of_an_asymmetric_ts_are_not_images():
    symbols, plus, minus = QRC_STARTS["hcn"]
    assert idn.carry(symbols, minus, plus, plus)[0] == pytest.approx(0.042, abs=5e-4)
    assert not idn.is_image(symbols, minus, plus)
    symbols, plus, minus = QRC_STARTS["sn2_cl"]
    bent = plus.copy()
    bent[[1, 3], 2] += 0.055  # two H atoms off the image: 0.025 A, like acac's ± starts
    rmsd = idn.carry(symbols, minus, bent, bent)[0]
    assert idn.IMAGE_A < rmsd == pytest.approx(0.025, abs=1e-3)
    assert not idn.is_image(symbols, minus, bent)


def test_carry_takes_other_into_the_atom_order_and_frame_of_ref():
    rotation, shift = _proper_rotation(3).T, np.array([1.0, -2.0, 0.5])
    asym = NH3 + [[0.0, 0.0, 0.0], [0.05, 0.0, 0.0], [0.0, 0.03, 0.0], [0.0, 0.0, -0.04]]
    order = [0, 3, 1, 2]  # a relabelled copy of asym, and its inversion in the same frame
    x, other = ((y @ rotation + shift)[order] for y in (asym, asym * [1.0, 1.0, -1.0]))
    rmsd, carried = idn.carry(NH3_SYMBOLS, asym, x, other)
    assert rmsd < 1e-6 and np.allclose(carried, asym * [1.0, 1.0, -1.0])
    x, other = (y * [-1.0, 1.0, 1.0] @ rotation + shift for y in (CHFCLBR, 1.1 * CHFCLBR))
    rmsd, carried = idn.carry(CHFCLBR_SYMBOLS, CHFCLBR, x, other)  # through the mirror image
    assert rmsd < 1e-6 and np.allclose(carried, 1.1 * CHFCLBR)


def test_periodic_nearest_wraps():
    assert idn.periodic_nearest(355.0, [0.0, 120.0, 240.0]) == 0
    assert idn.periodic_nearest(170.0, [-170.0, 90.0]) == 0
    with pytest.raises(ValueError):
        idn.periodic_nearest(0.0, [])


def test_the_radial_prefilter_is_a_lower_bound_of_the_basin_rmsd():
    """W4 (S8: thousands of TSs compared): the sorted centroid distances never differ by more
    than the permutation-invariant RMSD, so assign skips only candidates that cannot match."""
    from hfauto.chemistry.identity import _basin_match, _radii, _rms
    rng = np.random.default_rng(7)
    symbols = ["C", "H", "H", "H", "O", "H", "H"]
    x = rng.normal(size=(7, 3))
    for scale in (0.01, 0.03, 0.1, 0.5):
        y = (x + rng.normal(scale=scale, size=(7, 3)))[[0, 2, 1, 3, 4, 6, 5]]
        y = y @ np.linalg.qr(rng.normal(size=(3, 3)))[0]  # rotated (or mirrored)
        sym = [symbols[i] for i in [0, 2, 1, 3, 4, 6, 5]]
        assert _rms(_radii(symbols, x) - _radii(sym, y)) <= _basin_match(symbols, x, y)[0] + 1e-12

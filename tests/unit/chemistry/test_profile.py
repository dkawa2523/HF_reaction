"""Path class, peak interpolation and the one profile judgement (design §5.5, CH-08, U5-P1)."""

import numpy as np
import pytest

from hfauto.chemistry import profile as prof
from hfauto.chemistry.gates import Policy
from hfauto.core.constants import HARTREE_TO_KCAL_MOL


@pytest.mark.parametrize(
    ("energies", "expected"),
    [
        ([0.0, 0.4, 0.2, 0.9, 2.0, 3.0], "barrierless"),  # wiggle below the resolution
        ([0.0, 1.0, 2.0, 3.0], "barrierless"),
        ([0.0, 2.0, 5.0, 3.0, 1.0], "single"),
        ([0.0, 10.0, 9.99, 10.0, 0.0], "single"),  # shallow dip counts once
        ([0.0, 5.0, 1.0], "single"),  # three points: a low-level TS between two minima
        ([0.0, 0.5, 1.0], "barrierless"),
        ([0.0, 5.0, 2.0, 5.0, 0.0], "intermediate"),  # two peaks always have a well between
        ([0.0, -1.5, -3.0, -1.0, 1.0, 3.0, 5.0], "intermediate"),  # well below both ends
        ([0.0, 2.0, 0.5, 3.0, 6.0, 9.0, 10.0], "intermediate"),  # uphill two-step path
    ],
)
def test_classify_counts_hills_and_wells_alike(energies, expected):
    assert prof.classify(energies, resolution=1.0) == expected


def test_hei_interpolates_parabola_and_coordinates_at_the_given_peak():
    s = np.arange(6.0)
    energies = -((s - 2.3) ** 2)
    frames = [np.array([[x, 0.0, 0.0]]) for x in s]
    index, energy, coords = prof.hei(frames, energies, 2)
    assert index == pytest.approx(2.3) and energy == pytest.approx(0.0, abs=1e-12)
    assert coords == pytest.approx(np.array([[2.3, 0.0, 0.0]]))
    two = [0.0, 3.0, 1.0, 5.0, 0.0]  # the detected peak, not only the global maximum
    assert prof.hei(frames[:5], two, 1)[0] == pytest.approx(1.1)  # toward the higher side
    with pytest.raises(ValueError, match="interior"):
        prof.hei(frames, energies, 0)


K = 1.0 / HARTREE_TO_KCAL_MOL
POLICY = Policy()  # resolution 1 kcal/mol, spin_tol 0.1


def line(n):
    return [np.array([[float(x), 0.0, 0.0]]) for x in range(n)]


def test_judge_densifies_a_barrierless_profile_once_beside_its_highest_node():
    """U5-P1: new nodes at the midpoints of the two segments beside the highest interior node,
    asked by the node before each; a hill they find makes it a single step, judged once more
    and never densified again; a failed one leaves it unavailable."""
    path = prof.Profile(line(5), tuple(e * K for e in (0.0, 0.5, 0.8, 0.3, -1.0)), "string")
    asked = []

    def sample(new):
        asked.extend((i, float(x[0, 0])) for i, x in new)
        return [(0.9 * K, None), (2.0 * K, None)]

    verdict, judged = prof.judge(path, POLICY, sample)
    assert asked == [(1, 1.5), (2, 2.5)] and verdict.verdict == "single"
    assert [float(x[0, 0]) for x in judged.frames] == [0, 1, 1.5, 2, 2.5, 3, 4]
    assert judged.energies == tuple(e * K for e in (0.0, 0.5, 0.9, 0.8, 2.0, 0.3, -1.0))
    assert prof.judge(path, POLICY)[0].verdict == "barrierless"  # without an energy function
    failed, kept = prof.judge(path, POLICY, lambda new: [None] * len(new))
    assert (failed.verdict, failed.reasons, kept) == ("unavailable", ("midpoint_single_point",),
                                                      path)


def test_judge_puts_no_node_in_a_segment_of_zero_length():
    """An association's scan: the monomers' sum stands on the first scan frame, so only the
    segment after its highest (first) point gets a node."""
    frames = line(4)
    path = prof.Profile([frames[0], *frames], tuple(e * K for e in (0, -0.2, -0.4, -2, -9)),
                        "scan", (None, 1.7, 1.6, 1.2, 0.76))
    verdict, judged = prof.judge(path, POLICY, lambda new: [(-0.3 * K, 1.65)] * len(new))
    assert verdict.verdict == "barrierless" and len(judged.energies) == 6
    assert judged.s2 == (None, 1.7, 1.65, 1.6, 1.2, 0.76)


def test_judge_leaves_a_maximum_off_its_scf_branch_and_a_short_profile_unavailable():
    s2 = (1.7114, 0.7604, 0.7651, 0.7547)  # VAL7 S5's maximum beside the minimum's branch
    path = prof.Profile(line(4), tuple(e * K for e in (0.0, 26.0, 18.0, -35.0)), "screen", s2)
    assert prof.judge(path, POLICY)[0].reasons == ("scf_branch_jump",)
    assert prof.judge(path._replace(s2=()), POLICY)[0].verdict == "single"  # not observed
    short = prof.Profile(line(2), (0.0, -1.0), "string")
    assert prof.judge(short, POLICY)[0].reasons == ("too_few_points",)

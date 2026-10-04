"""DFT path profiles: one type and one judgement for every source (§5.5).

A path is classified as it stands, converged or not: its maximum between two minima bounds the
saddle from above, and its peak is only a seed.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Literal, NamedTuple

import numpy as np

from hfauto.chemistry.gates import Policy
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.records import BarrierVerdict, ProfileSource


def interior_maxima(energies: Sequence[float], resolution: float) -> tuple[int, ...]:
    """Indices of interior maxima that rise and then fall by at least `resolution`.

    Wiggles smaller than the resolution are ignored, so two peaks separated by a shallow
    dip count once.
    """

    peaks: list[int] = []
    low = float(energies[0])
    candidate: tuple[float, int] | None = None
    for index, value in enumerate(map(float, energies)):
        if candidate is None:
            if value < low:
                low = value
            elif value - low >= resolution:
                candidate = (value, index)
        elif value > candidate[0]:
            candidate = (value, index)
        elif candidate[0] - value >= resolution:
            peaks.append(candidate[1])
            low, candidate = value, None
    return tuple(peaks)


PathClass = Literal["barrierless", "single", "intermediate"]


def classify(energies: Sequence[float], resolution: float) -> PathClass:
    """Hills and wells counted by the same rule: a well (an interior minimum a resolution deep
    on both sides) marks an intermediate, and two peaks always have one between them; without
    a well one peak is a single step, and without either the path is barrierless."""

    if interior_maxima([-float(e) for e in energies], resolution):
        return "intermediate"
    return "single" if interior_maxima(energies, resolution) else "barrierless"


def _frames(coords_list: Sequence[np.ndarray]) -> np.ndarray:
    frames = np.asarray([np.asarray(c, dtype=float).reshape(-1, 3) for c in coords_list])
    if len(frames) < 2:
        raise ValueError("a path needs at least two images")
    return frames


def hei(
    coords_list: Sequence[np.ndarray], energies: Sequence[float], k: int
) -> tuple[float, float, np.ndarray]:
    """Interior image ``k`` (a detected peak) refined by a parabola through it and its
    neighbours.

    Returns (fractional index, interpolated energy, coordinates interpolated linearly at
    the fractional index).
    """

    frames, e = _frames(coords_list), np.asarray(energies, dtype=float)
    if len(e) != len(frames) or not 0 < k < len(e) - 1:
        raise ValueError("need matching coordinates and energies and an interior image")
    left, mid, right = e[k - 1], e[k], e[k + 1]
    curvature = left - 2.0 * mid + right
    offset = 0.0
    if curvature < 0:
        offset = float(np.clip(0.5 * (left - right) / curvature, -0.5, 0.5))
    energy = mid + 0.5 * (right - left) * offset + 0.5 * curvature * offset**2
    neighbour = k + 1 if offset >= 0 else k - 1
    coords = frames[k] + abs(offset) * (frames[neighbour] - frames[k])
    return k + offset, float(energy), coords


class Profile(NamedTuple):
    """A DFT profile: the nodes of a continuous path from end to end, their energies and ⟨S²⟩
    (empty when not observed, as a string's beads), and what made it. The ends are the DFT
    minima; an association's scan starts at its separated monomers' energy sum, which stands
    on the first scan point's frame."""

    frames: list[np.ndarray]
    energies: tuple[float, ...]
    source: ProfileSource
    s2: tuple[float | None, ...] = ()


Point = tuple[float, float | None]  # a new node's (energy, ⟨S²⟩)
# a profile's energy function at new nodes: (index of the node before it, structure) -> its
# Point, None when its SP failed
Sample = Callable[[list[tuple[int, np.ndarray]]], Sequence[Point | None]]


def branch_jump(energies: Sequence[float], s2: Sequence[float | None], tol: float) -> bool:
    """Whether a profile's highest interior point lies on another SCF branch than its
    neighbours: its ⟨S²⟩ outside theirs and more than ``tol`` (Policy.spin_tol) from one of
    them. At a genuine radical TS the UKS contamination peaks smoothly, at most 0.022 above a
    neighbour (H + H2 0.7674 beside 0.7565 and 0.7570, the OH + CH4 shortcut 0.7722 beside
    0.7544 and 0.7500), and through the Coulson–Fischer region it changes monotonically (CH3 +
    O2 1.7115 → 0.7545, its 0.80 kcal/mol hill at 1.5438 between 1.6261 and 1.4145); an SCF
    branch jump moves it 0.77–0.95 (CH3 + O2 0.7604 beside 1.7114; from an atomic guess 0.7591
    beside 1.7115, from mixed guesses 0.7697 beside 1.5438). An unobserved ⟨S²⟩ is no evidence."""
    if not s2:
        return False
    k = 1 + int(np.argmax(energies[1:-1]))
    before, at, after = s2[k - 1:k + 2]
    if before is None or at is None or after is None:
        return False
    outside = not min(before, after) <= at <= max(before, after)
    return outside and max(abs(at - before), abs(at - after)) > tol


def _verdict(path: Profile, policy: Policy) -> BarrierVerdict:
    def unavailable(reason: str) -> BarrierVerdict:
        return BarrierVerdict(verdict="unavailable", source=path.source, reasons=(reason,))

    if len(path.energies) < 3:
        return unavailable("too_few_points")
    if branch_jump(path.energies, path.s2, policy.spin_tol):
        return unavailable("scf_branch_jump")
    resolution = policy.resolution_kcal / HARTREE_TO_KCAL_MOL
    return BarrierVerdict(verdict=classify(path.energies, resolution), source=path.source)


def judge(path: Profile, policy: Policy, sample: Sample | None = None
          ) -> tuple[BarrierVerdict, Profile]:
    """The class of ``path`` at resolution_kcal and the profile it was judged on. A maximum off
    its SCF branch (``branch_jump``) leaves it unavailable. A barrierless class is accepted
    only after densifying with ``sample``: new nodes at the midpoints of the two segments beside
    the highest interior node (none in a segment of zero length: the monomers' sum standing on
    the first scan frame), then judged again; a failed SP leaves it unavailable. A harmonic top
    hidden between nodes h apart (measured 0.04-0.30 Å) rises at most h²|E''|/8 above the higher
    one, and it lies beside the highest node only on a unimodal profile: a hill hidden
    elsewhere is not looked for."""
    verdict = _verdict(path, policy)
    if verdict.verdict != "barrierless" or sample is None:
        return verdict, path
    k = 1 + int(np.argmax(path.energies[1:-1]))
    frames, energies, s2 = list(path.frames), list(path.energies), list(path.s2)
    new = [(i, 0.5 * (frames[i] + frames[i + 1])) for i in (k - 1, k)
           if not np.array_equal(frames[i], frames[i + 1])]
    points = [p for p in sample(new) if p is not None]
    if len(points) < len(new):
        return BarrierVerdict(verdict="unavailable", source=path.source,
                              reasons=("midpoint_single_point",)), path
    for (i, x), (energy, spin) in reversed(list(zip(new, points, strict=True))):
        frames.insert(i + 1, x)
        energies.insert(i + 1, energy)
        if s2:
            s2.insert(i + 1, spin)
    path = path._replace(frames=frames, energies=tuple(energies), s2=tuple(s2))
    return _verdict(path, policy), path

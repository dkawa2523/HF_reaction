"""Path profiles: class and peak interpolation (§5.5).

A path is classified as it stands, converged or not: its maximum between two minima bounds the
saddle from above, and its peak is only a seed.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Literal

import numpy as np


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

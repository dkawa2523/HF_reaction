"""Path profiles: shape, highest energy image, convergence and stagnation (§5.5).

Convergence is judged from the gradient history only; the program's own "converged"
message is never used (CH-08).
"""

from __future__ import annotations

from collections.abc import Sequence
from itertools import pairwise
from typing import Literal

import numpy as np

# Xmax below this (the engine's step unit) counts as a zero step; exact zeros are rare.
_XMAX_ZERO = 1.0e-6


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


def shape(
    energies: Sequence[float], resolution: float
) -> Literal["monotonic", "single_max", "multi_max"]:
    count = len(interior_maxima(energies, resolution))
    return "monotonic" if count == 0 else "single_max" if count == 1 else "multi_max"


def _frames(coords_list: Sequence[np.ndarray]) -> np.ndarray:
    frames = np.asarray([np.asarray(c, dtype=float).reshape(-1, 3) for c in coords_list])
    if len(frames) < 2:
        raise ValueError("a path needs at least two images")
    return frames


def hei(
    coords_list: Sequence[np.ndarray], energies: Sequence[float]
) -> tuple[float, float, np.ndarray]:
    """Highest interior image refined by a parabola through it and its neighbours.

    Returns (fractional index, interpolated energy, coordinates interpolated linearly at
    the fractional index).
    """

    frames, e = _frames(coords_list), np.asarray(energies, dtype=float)
    if len(e) != len(frames) or len(e) < 3:
        raise ValueError("need matching coordinates and energies for at least three images")
    k = 1 + int(np.argmax(e[1:-1]))
    left, mid, right = e[k - 1], e[k], e[k + 1]
    curvature = left - 2.0 * mid + right
    offset = 0.0
    if curvature < 0:
        offset = float(np.clip(0.5 * (left - right) / curvature, -0.5, 0.5))
    energy = mid + 0.5 * (right - left) * offset + 0.5 * curvature * offset**2
    neighbour = k + 1 if offset >= 0 else k - 1
    coords = frames[k] + abs(offset) * (frames[neighbour] - frames[k])
    return k + offset, float(energy), coords


def string_converged(
    gmax_history: Sequence[float], gmax_tol: float = 1.0e-3, non_worsening: int = 3
) -> bool:
    """Final gmax ≤ tol and no increase over the last `non_worsening` iterations."""

    if not gmax_history or gmax_history[-1] > gmax_tol:
        return False
    tail = list(gmax_history[-(non_worsening + 1) :])
    return all(later <= earlier for earlier, later in pairwise(tail))


def stagnated(
    gmax_history: Sequence[float],
    xmax_history: Sequence[float] = (),
    window: int = 10,
    min_improvement: float = 0.05,
    xmax_zero_run: int = 5,
) -> bool:
    """True when the last `window` iterations improved the best gmax by less than
    `min_improvement` (relative), or the last `xmax_zero_run` steps were zero."""

    if len(gmax_history) > window:
        best_before = min(gmax_history[:-window])
        if min(gmax_history[-window:]) > (1.0 - min_improvement) * best_before:
            return True
    tail = list(xmax_history[-xmax_zero_run:])
    return len(tail) == xmax_zero_run and all(abs(x) < _XMAX_ZERO for x in tail)


def max_node_spacing(coords_list: Sequence[np.ndarray]) -> float:
    """Largest Cartesian distance (Å, all atoms) between neighbouring images."""

    frames = _frames(coords_list)
    steps = (frames[1:] - frames[:-1]).reshape(len(frames) - 1, -1)
    return float(np.linalg.norm(steps, axis=1).max())


def tangent(coords_list: Sequence[np.ndarray], index: int) -> np.ndarray:
    """Unit (N, 3) tangent: central difference inside the path, one-sided at its ends."""

    frames = _frames(coords_list)
    lo, hi = max(index - 1, 0), min(index + 1, len(frames) - 1)
    step = frames[hi] - frames[lo]
    norm = float(np.linalg.norm(step))
    if norm == 0.0:
        raise ValueError("coincident neighbouring images have no tangent")
    return step / norm

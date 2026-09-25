"""Displacements along normal modes, mode overlaps and mode-follow outcomes (§5.5)."""

from __future__ import annotations

from typing import Literal

import numpy as np

from hfauto.core.constants import BOHR_TO_ANGSTROM


def _unit_max_displacement(mode: np.ndarray) -> np.ndarray:
    """Mode as (N, 3) scaled so that the largest atomic displacement is 1 Å."""

    u = np.asarray(mode, dtype=float).reshape(-1, 3)
    largest = float(np.linalg.norm(u, axis=1).max())
    if largest == 0.0 or not np.isfinite(largest):
        raise ValueError("mode has no finite displacement")
    return u / largest


def displace(
    coords: np.ndarray, mode: np.ndarray, amplitude_A: float
) -> tuple[np.ndarray, np.ndarray]:
    """coords ± mode, with the largest atomic displacement equal to amplitude_A."""

    x = np.asarray(coords, dtype=float).reshape(-1, 3)
    step = amplitude_A * _unit_max_displacement(mode)
    if step.shape != x.shape:
        raise ValueError("mode and coordinates differ in atom count")
    return x + step, x - step


def overlap(a: np.ndarray, b: np.ndarray) -> float:
    """|cos| between two displacement vectors (mode signs are arbitrary)."""

    va, vb = np.asarray(a, dtype=float).ravel(), np.asarray(b, dtype=float).ravel()
    return float(abs(va @ vb) / (np.linalg.norm(va) * np.linalg.norm(vb)))


def qrc_amplitude(
    hessian: np.ndarray,
    mode: np.ndarray,
    *,
    target_hartree: float,
    bounds_A: tuple[float, float] = (0.05, 0.4),
) -> float:
    """Largest atomic displacement s with ½|κ|s² = target, clipped to bounds_A.

    κ = uᵀHu (Eh/Å²) for the mode scaled to a 1 Å largest atomic displacement; the
    Hessian is in Eh/bohr².
    """

    u = _unit_max_displacement(mode).ravel()
    h = np.asarray(hessian, dtype=float)
    if h.shape != (u.size, u.size):
        raise ValueError("Hessian and mode differ in dimension")
    curvature = abs(float(u @ h @ u)) / BOHR_TO_ANGSTROM**2
    low, high = bounds_A
    if curvature == 0.0:
        return high
    return float(np.clip(np.sqrt(2.0 * target_hartree / curvature), low, high))


def classify_mode_follow(
    source: str | None, plus: str | None, minus: str | None
) -> Literal["replace", "ts_candidate", "one_side", "same_as_source"]:
    """Classify where the ± displaced relaxations ended up.

    Arguments are the identities (basin ids or labels) of the source structure and of the
    minima reached from each side; None means that side did not converge.
    """

    reached = [side for side in (plus, minus) if side is not None and side != source]
    if not reached:
        return "same_as_source"
    if len(reached) == 1:
        return "one_side"
    return "replace" if reached[0] == reached[1] else "ts_candidate"

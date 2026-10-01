"""Displacements along normal modes, mode overlaps and mode-follow outcomes (§5.5).

A displacement off a stationary point is ``qrc_step`` (along one imaginary mode, at one energy
target) or ``off_saddle`` (along the sum of the imaginary modes below a threshold, each by its
``qrc_step``); both cap the largest atomic displacement at BOUNDS_A[1]. The ± sides of QRC from
a TS and of the mode-follow of a first-order saddle are therefore one displacement.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Literal

import numpy as np

from hfauto.chemistry.elements import mass
from hfauto.chemistry.gates import qrc_drop
from hfauto.core.constants import AMU_TO_ME, BOHR_TO_ANGSTROM, CM1_TO_HARTREE
from hfauto.core.evidence import Evidence

TARGET_HARTREE = 3.0e-4  # expected energy change of one displacement (QRC and mode-follow)
BOUNDS_A = (0.05, 0.4)  # largest atomic displacement (Å); an opt reuses a Hessian within 0.5 Å


def _unit_max_displacement(mode: np.ndarray) -> np.ndarray:
    """Mode as (N, 3) scaled so that the largest atomic displacement is 1 Å."""

    u = np.asarray(mode, dtype=float).reshape(-1, 3)
    largest = float(np.linalg.norm(u, axis=1).max())
    if largest == 0.0 or not np.isfinite(largest):
        raise ValueError("mode has no finite displacement")
    return u / largest


def displace(
    coords: np.ndarray, mode: np.ndarray, largest_A: float
) -> tuple[np.ndarray, np.ndarray]:
    """coords ± mode, with the largest atomic displacement equal to largest_A."""

    x = np.asarray(coords, dtype=float).reshape(-1, 3)
    step = largest_A * _unit_max_displacement(mode)
    if step.shape != x.shape:
        raise ValueError("mode and coordinates differ in atom count")
    return x + step, x - step


def overlap(a: np.ndarray, b: np.ndarray) -> float:
    """|cos| between two displacement vectors (mode signs are arbitrary)."""

    va, vb = np.asarray(a, dtype=float).ravel(), np.asarray(b, dtype=float).ravel()
    return float(abs(va @ vb) / (np.linalg.norm(va) * np.linalg.norm(vb)))


def amplitude(
    nu_cm1: float,
    mode: np.ndarray,
    symbols: Sequence[str],
    target_hartree: float = TARGET_HARTREE,
    bounds_A: tuple[float, float] = BOUNDS_A,
) -> float:
    """Largest atomic displacement s (Å) with ½κs² = target along ``mode``, clipped to bounds_A.

    κ = ω²·Σ m_i|u_i|² for the mode u scaled to a 1 Å largest atomic displacement: the harmonic
    energy E = ½ω²Q², equal to uᵀHu for an exact normal mode (QRC, Goodman & Silva 2003).
    """

    u = _unit_max_displacement(mode)
    masses = np.array([mass(symbol) for symbol in symbols])
    if len(masses) != len(u):
        raise ValueError("mode and symbols differ in atom count")
    omega2 = (nu_cm1 * CM1_TO_HARTREE / BOHR_TO_ANGSTROM) ** 2 * AMU_TO_ME  # Eh/(amu Å²)
    curvature = omega2 * float(masses @ (u * u).sum(axis=1))
    low, high = bounds_A
    if curvature == 0.0:
        return high
    return float(np.clip(np.sqrt(2.0 * target_hartree / curvature), low, high))


def capped(step: np.ndarray) -> np.ndarray:
    """``step`` (N, 3) with its largest atomic displacement capped at BOUNDS_A[1]."""

    largest = float(np.linalg.norm(step, axis=1).max())
    return step * min(1.0, BOUNDS_A[1] / largest) if largest > 0.0 else step


def qrc_step(freq: Evidence, index: int, symbols: Sequence[str], factor: float = 1.0
             ) -> np.ndarray:
    """The displacement (N, 3) along imaginary mode ``index`` of ``freq`` (QRC, Goodman & Silva
    2003): its ``amplitude`` for the energy target max(3 x qrc_drop, TARGET_HARTREE), as a QRC
    side must end a qrc_drop below the TS, times ``factor``, capped at BOUNDS_A[1]. With the
    bundled methods (scf_energy_tol 1e-7) 3 x qrc_drop is 3e-5 Eh: the target is TARGET_HARTREE."""

    mode = np.asarray(freq.imaginary_modes[index])
    target = max(3.0 * qrc_drop(freq.level), TARGET_HARTREE)
    s = amplitude((freq.frequencies_cm1 or ())[index], mode, symbols, target_hartree=target)
    return min(s * factor, BOUNDS_A[1]) * _unit_max_displacement(mode)


def off_saddle(freq: Evidence, symbols: Sequence[str], *, below_cm1: float,
               keep: int | None = None) -> np.ndarray:
    """One side off a stationary point: the sum of its imaginary modes below -below_cm1 but
    ``keep`` (frequencies ascending, one mode per negative one), each by its ``qrc_step``,
    capped at BOUNDS_A[1]; with one such mode, its qrc_step."""

    picked = [i for i, nu in enumerate(freq.frequencies_cm1 or ()) if nu < -below_cm1 and i != keep]
    return capped(sum((qrc_step(freq, i, symbols) for i in picked), np.zeros((len(symbols), 3))))


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

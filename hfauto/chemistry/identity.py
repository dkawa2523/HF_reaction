"""Permutation-invariant structure identity (§5.5).

Same-element atoms are matched with scipy's Hungarian solver, alternating with a proper Kabsch
fit from several starting orientations (the alternation only finds a local optimum). A basin
(``same_basin``, ``assign``) admits proper and improper rotations: mirror images are one minimum,
a chiral one with m = 2 (``is_chiral``). Labelled comparisons (``mapped_rmsd``) stay proper, so
that a degenerate rearrangement such as the NH3 inversion differs from the identity.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

import numpy as np
from scipy.optimize import linear_sum_assignment

_RMSD_A = 0.05  # one basin: permutation-invariant RMSD (mirror image included) ...
_DE_HARTREE = 5.0e-5  # ... and |dE|, with the best match clearly ahead of the runner-up
_RUNNER_UP_RATIO, _RUNNER_UP_GAP_A = 3.0, 0.1
_MAX_REFINE = 5
_DEGENERATE_REL = 0.05  # principal moments this close (relative) count as degenerate
_IN_PLANE_STEP_DEG = 30
# Proper sign flips of the principal axes.
_SIGN_FLIPS = tuple(np.diag(s) for s in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)))


def _centered(coords: np.ndarray) -> np.ndarray:
    x = np.asarray(coords, dtype=float).reshape(-1, 3)
    return x - x.mean(axis=0)


def _mirror(coords: np.ndarray) -> np.ndarray:
    return np.asarray(coords, dtype=float).reshape(-1, 3) * [-1.0, 1.0, 1.0]


def _rotation(mobile: np.ndarray, fixed: np.ndarray) -> np.ndarray:
    """Proper rotation R minimizing |mobile @ R − fixed| (both centered)."""

    u, _, vt = np.linalg.svd(mobile.T @ fixed)
    handedness = 1.0 if np.linalg.det(u @ vt) >= 0 else -1.0
    return u @ np.diag([1.0, 1.0, handedness]) @ vt


def _rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def _principal_frame(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Ascending moments of the unit-mass inertia tensor and a right-handed axis matrix."""

    moments, axes = np.linalg.eigh(np.trace(x.T @ x) * np.eye(3) - x.T @ x)
    if np.linalg.det(axes) < 0:
        axes[:, 2] *= -1.0
    return moments, axes


def _spins(moments: np.ndarray) -> list[np.ndarray]:
    """In-plane rotations (principal frame) for each degenerate pair of moments."""

    spins = [np.eye(3)]
    scale = max(float(moments.max()), 1e-12)
    for i, j in ((1, 2), (0, 2), (0, 1)):
        if abs(moments[i] - moments[j]) > _DEGENERATE_REL * scale:
            continue
        for step in range(1, 360 // _IN_PLANE_STEP_DEG):
            angle = math.radians(step * _IN_PLANE_STEP_DEG)
            spin = np.eye(3)
            spin[i, i] = spin[j, j] = math.cos(angle)
            spin[i, j], spin[j, i] = -math.sin(angle), math.sin(angle)
            spins.append(spin)
    return spins


def _assign(symbols: Sequence[str], a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """perm with b[perm[i]] matched to a[i], same elements only."""

    perm = np.arange(len(symbols))
    elements = np.asarray(symbols)
    for element in set(symbols):
        idx = np.flatnonzero(elements == element)
        cost = np.sum((a[idx, None, :] - b[None, idx, :]) ** 2, axis=-1)
        rows, cols = linear_sum_assignment(cost)
        perm[idx[rows]] = idx[cols]
    return perm


def _check(symbols: Sequence[str], a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    xa, xb = _centered(a), _centered(b)
    if xa.shape != xb.shape or len(xa) != len(symbols):
        raise ValueError("structures must have the same atoms as symbols")
    return xa, xb


def permutation_invariant_rmsd(
    symbols: Sequence[str], a: np.ndarray, b: np.ndarray
) -> tuple[float, np.ndarray]:
    """(rmsd, perm) such that b[perm] is the best proper-rotation match of a."""

    xa, xb = _check(symbols, a, b)
    moments, frame_a = _principal_frame(xa)
    _, frame_b = _principal_frame(xb)
    starts = [_rotation(xb, xa)]  # identity mapping first
    spins = _spins(moments)
    starts += [frame_b @ flip @ spin @ frame_a.T for flip in _SIGN_FLIPS for spin in spins]
    best = (math.inf, np.arange(len(xa)))
    for start in starts:
        perm = _assign(symbols, xa, xb @ start)
        for _ in range(_MAX_REFINE):
            new_perm = _assign(symbols, xa, xb @ _rotation(xb[perm], xa))
            if np.array_equal(new_perm, perm):
                break
            perm = new_perm
        value = _rmsd(xb[perm] @ _rotation(xb[perm], xa), xa)
        if value < best[0]:
            best = (value, perm)
    return best


def mapped_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """Kabsch RMSD with the identity atom mapping (proper rotations only)."""

    xa, xb = _centered(a), _centered(b)
    if xa.shape != xb.shape:
        raise ValueError("structures differ in atom count")
    return _rmsd(xb @ _rotation(xb, xa), xa)


def _basin_match(symbols: Sequence[str], a: np.ndarray, b: np.ndarray
                 ) -> tuple[float, np.ndarray, bool]:
    """(rmsd, perm, mirrored): the better of b and its mirror image matched onto a."""

    rmsd, perm = permutation_invariant_rmsd(symbols, a, b)
    m_rmsd, m_perm = permutation_invariant_rmsd(symbols, a, _mirror(b))
    return (m_rmsd, m_perm, True) if m_rmsd < rmsd else (rmsd, perm, False)


def is_chiral(symbols: Sequence[str], x: np.ndarray) -> bool:
    """True when x and its mirror image are distinct under proper rotations and relabelling
    (optical isomer number m = 2)."""

    return permutation_invariant_rmsd(symbols, x, _mirror(x))[0] > _RMSD_A


def assign(
    symbols: Sequence[str],
    coords: np.ndarray,
    energy: float,
    candidates: Mapping[str, tuple[np.ndarray, float]],
) -> str | None:
    """Id of the candidate (coords, energy) whose basin holds the structure, else None.

    Candidates within 5e-5 Eh come first; the closest of them (mirror image included) must lie
    within 0.05 A and clearly ahead of the runner-up."""

    scored = sorted((_basin_match(symbols, coords, xyz)[0], key)
                    for key, (xyz, e) in candidates.items() if abs(energy - e) <= _DE_HARTREE)
    if not scored or scored[0][0] > _RMSD_A:
        return None
    if len(scored) > 1:
        best, runner = scored[0][0], scored[1][0]
        separated = runner - best >= _RUNNER_UP_GAP_A or (
            runner >= _RUNNER_UP_RATIO * best and runner > best)
        if not separated:
            return None
    return scored[0][1]


def same_basin(symbols: Sequence[str], a: np.ndarray, b: np.ndarray, ea: float, eb: float
               ) -> bool:
    """Structures a and b (energies ea, eb) lie in one basin by assign's criterion."""

    return assign(symbols, a, ea, {"b": (b, eb)}) is not None


def basin_coords(symbols: Sequence[str], basin: np.ndarray, own: np.ndarray) -> np.ndarray:
    """The basin's structure in the atom order of ``own`` (a structure of that basin), mirrored
    when the basin is chiral and ``own`` has the other handedness."""

    x = np.asarray(basin, dtype=float).reshape(-1, 3)
    if not is_chiral(symbols, x):
        return x[permutation_invariant_rmsd(symbols, own, x)[1]]
    _, perm, mirrored = _basin_match(symbols, own, x)
    return (_mirror(x) if mirrored else x)[perm]


def mapped_equivalent(symbols: Sequence[str], a: np.ndarray, b: np.ndarray) -> bool:
    """True when a and b differ as labelled (proper rotations) but are one basin structure
    (mirror image included): a degenerate rearrangement such as the NH3 inversion or the
    enantiomerization of a chiral minimum."""

    return mapped_rmsd(a, b) > _RMSD_A and _basin_match(symbols, a, b)[0] <= _RMSD_A


def periodic_nearest(value_deg: float, targets_deg: Sequence[float]) -> int:
    """Index of the target closest to value on a 360° circle (first one on ties)."""

    if not targets_deg:
        raise ValueError("no targets")
    distances = [abs((value_deg - t + 180.0) % 360.0 - 180.0) for t in targets_deg]
    return int(np.argmin(distances))

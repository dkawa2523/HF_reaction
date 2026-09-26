"""Mapped alignment of endpoints and paths, and IDPP initial paths (§5.5, CH-09).

IDPP (Smidstrup et al., JCP 140, 214106 (2014)): every image is pulled toward the linearly
interpolated interatomic distances, S = Σ_{i<j} d⁻⁴ (d_target − d)², with nudged-elastic-band
projection so images stay spread along the path.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from hfauto.chemistry.geometry import align_coordinates

MIN_DISTANCE_A = 0.7
_MAX_ITERATIONS = 5000
_FORCE_TOL = 1.0e-3
_SPRING = 1.0
_STEP = 0.1
_MAX_STEP_A = 0.05
# Coherent sin(πt) kick that lets exactly collinear endpoints leave the line.
_SYMMETRY_BREAK_A = 0.01


def align_mapped(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """b rigidly aligned onto a with the identity atom mapping (proper rotation)."""

    return align_coordinates(np.asarray(a, dtype=float), np.asarray(b, dtype=float))


def align_sequential(frames: Sequence[np.ndarray]) -> list[np.ndarray]:
    """Each frame aligned onto the already aligned previous one (the first is kept), so a
    path from internal coordinates carries no rigid jumps between neighbouring images."""

    aligned = [np.asarray(frames[0], dtype=float)]
    for frame in frames[1:]:
        aligned.append(align_mapped(aligned[-1], frame))
    return aligned


def _distances(x: np.ndarray) -> np.ndarray:
    return np.linalg.norm(x[:, None] - x[None], axis=-1)


def min_interatomic_distance(coords: np.ndarray) -> float:
    d = _distances(np.asarray(coords, dtype=float).reshape(-1, 3))
    return float(d[np.triu_indices(len(d), 1)].min())


def _idpp_gradient(x: np.ndarray, target: np.ndarray) -> np.ndarray:
    diff = x[:, None] - x[None]
    d = np.maximum(np.linalg.norm(diff, axis=-1), 1.0e-3)
    residual = target - d
    ds_dd = -4.0 * residual**2 / d**5 - 2.0 * residual / d**4
    np.fill_diagonal(ds_dd, 0.0)
    return np.einsum("ij,ijk->ik", ds_dd / d, diff)


def _neb_forces(images: np.ndarray, targets: np.ndarray) -> np.ndarray:
    forces = np.zeros_like(images)
    for k in range(1, len(images) - 1):
        tau = (images[k + 1] - images[k - 1]).ravel()
        tau /= max(float(np.linalg.norm(tau)), 1e-12)
        grad = _idpp_gradient(images[k], targets[k]).ravel()
        spring = _SPRING * (
            np.linalg.norm(images[k + 1] - images[k]) - np.linalg.norm(images[k] - images[k - 1])
        )
        forces[k] = (-(grad - (grad @ tau) * tau) + spring * tau).reshape(-1, 3)
    return forces


def _check_contacts(images: Sequence[np.ndarray]) -> None:
    for index, image in enumerate(images):
        shortest = min_interatomic_distance(image)
        if shortest < MIN_DISTANCE_A:
            raise ValueError(
                f"IDPP image {index} has an interatomic distance of {shortest:.3f} Å "
                f"(< {MIN_DISTANCE_A} Å)"
            )


def idpp(symbols: Sequence[str], a: np.ndarray, b: np.ndarray, n_images: int) -> list[np.ndarray]:
    """n_images geometries from a to b (aligned onto a first); the endpoints are fixed.

    Raises ValueError if any image has an interatomic distance below 0.7 Å.
    """

    start = np.asarray(a, dtype=float).reshape(-1, 3)
    end = align_mapped(start, np.asarray(b, dtype=float).reshape(-1, 3))
    if len(start) != len(symbols) or n_images < 2:
        raise ValueError("idpp needs matching symbols and at least two images")
    _check_contacts([start, end])
    t = np.linspace(0.0, 1.0, n_images)[:, None, None]
    kick = np.random.default_rng(0).standard_normal(start.shape)
    kick -= kick.mean(axis=0)
    kick *= _SYMMETRY_BREAK_A / np.linalg.norm(kick, axis=1).max()
    images = start + t * (end - start) + np.sin(np.pi * t) * kick
    images[0], images[-1] = start, end
    d0, d1 = _distances(start), _distances(end)
    targets = d0 + t * (d1 - d0)
    for _ in range(_MAX_ITERATIONS):
        forces = _neb_forces(images, targets)
        if np.linalg.norm(forces, axis=2).max() < _FORCE_TOL:
            break
        step = _STEP * forces
        largest = np.linalg.norm(step, axis=2).max(axis=1)
        images += step * np.minimum(1.0, _MAX_STEP_A / np.maximum(largest, 1e-300))[:, None, None]
    _check_contacts(images)
    return [image.copy() for image in images]

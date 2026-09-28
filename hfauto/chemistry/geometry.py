"""Small geometry primitives shared by chemistry workflows, and the value and gradient of a
declared reaction coordinate."""

from __future__ import annotations

from collections.abc import Collection, Sequence

import numpy as np

from hfauto.core.records import CoordinateTerm

_FD_STEP_A = 1.0e-4  # central-difference step of a declared coordinate's gradient


def align_coordinates(reference: np.ndarray, moving: np.ndarray) -> np.ndarray:
    """Rigidly align ``moving`` onto ``reference`` (same atom order, Kabsch)."""

    fixed = np.asarray(reference, dtype=float)
    mobile = np.asarray(moving, dtype=float)
    if (
        fixed.shape != mobile.shape
        or fixed.ndim != 2
        or fixed.shape[1] != 3
        or len(fixed) == 0
        or not np.isfinite(fixed).all()
        or not np.isfinite(mobile).all()
    ):
        raise ValueError("coordinate arrays must be finite, non-empty Nx3 arrays")
    fixed_center = fixed.mean(axis=0)
    mobile_center = mobile.mean(axis=0)
    try:
        u, _singular_values, vt = np.linalg.svd((mobile - mobile_center).T @ (fixed - fixed_center))
    except np.linalg.LinAlgError as exc:
        raise ValueError("Kabsch alignment failed") from exc
    rotation = u @ np.diag([1.0, 1.0, np.linalg.det(u @ vt)]) @ vt
    return (mobile - mobile_center) @ rotation + fixed_center


def dihedral_deg(coords: np.ndarray, atoms: Sequence[int]) -> float:
    """Dihedral angle of the atoms i-j-k-l in degrees, in (-180, 180]."""

    p = np.asarray(coords, dtype=float).reshape(-1, 3)[list(atoms)]
    b0, b1, b2 = p[0] - p[1], p[2] - p[1], p[3] - p[2]
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - (b0 @ b1) * b1, b2 - (b2 @ b1) * b1
    return float(np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w)))


def declared_coordinate(terms: Sequence[CoordinateTerm], x: np.ndarray) -> float:
    """Σ coefficient × (distance Å | angle ° | dihedral °)."""
    total = 0.0
    for term in terms:
        p = np.asarray(x, dtype=float).reshape(-1, 3)[list(term.atoms)]
        if term.kind == "distance":
            value = float(np.linalg.norm(p[1] - p[0]))
        elif term.kind == "angle":
            u, v = p[0] - p[1], p[2] - p[1]
            value = float(np.degrees(np.arccos(np.clip(
                u @ v / (np.linalg.norm(u) * np.linalg.norm(v)), -1.0, 1.0))))
        else:
            value = dihedral_deg(x, term.atoms)
        total += term.coefficient * value
    return total


def declared_coordinate_gradient(terms: Sequence[CoordinateTerm], x: np.ndarray) -> np.ndarray:
    """Central-difference gradient of ``declared_coordinate``; dihedral steps wrap at ±180°."""
    flat, grad = np.asarray(x, dtype=float).ravel(), np.zeros(np.size(x))
    for i in range(flat.size):
        step = np.zeros_like(flat)
        step[i] = _FD_STEP_A
        delta = declared_coordinate(terms, flat + step) - declared_coordinate(terms, flat - step)
        grad[i] = ((delta + 180.0) % 360.0 - 180.0) / (2 * _FD_STEP_A)
    return grad


def most_changed_dihedral(bonded: Collection[tuple[int, int]], a: np.ndarray, b: np.ndarray
                          ) -> tuple[CoordinateTerm, ...]:
    """The dihedral along bonds i-j, j-k, k-l whose angle differs most (periodically) between
    the structures ``a`` and ``b``, as a declared coordinate; () when the bond graph has no such
    chain."""

    neighbours: dict[int, set[int]] = {}
    for i, j in bonded:
        neighbours.setdefault(i, set()).add(j)
        neighbours.setdefault(j, set()).add(i)
    chains = [(i, j, k, m) for j, k in bonded for i in neighbours[j] - {k}
              for m in neighbours[k] - {i, j}]
    if not chains:
        return ()
    change = [abs((dihedral_deg(b, c) - dihedral_deg(a, c) + 180.0) % 360.0 - 180.0)
              for c in chains]
    return (CoordinateTerm(kind="dihedral", atoms=chains[int(np.argmax(change))]),)

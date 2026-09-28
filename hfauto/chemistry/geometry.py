"""Geometry primitives of chemistry (Kabsch rotation, rotation about an axis, dihedral angle,
bonded neighbours) and the value and gradient of a declared reaction coordinate."""

from __future__ import annotations

from collections.abc import Collection, Sequence

import numpy as np

from hfauto.core.records import CoordinateTerm

_FD_STEP_A = 1.0e-4  # central-difference step of a declared coordinate's gradient


def kabsch(mobile: np.ndarray, fixed: np.ndarray) -> np.ndarray:
    """Proper rotation R minimizing |mobile @ R − fixed| (both centered, N×3)."""

    u, _, vt = np.linalg.svd(mobile.T @ fixed)
    handedness = 1.0 if np.linalg.det(u @ vt) >= 0 else -1.0
    return u @ np.diag([1.0, 1.0, handedness]) @ vt


def rotation_about(axis: np.ndarray, degrees: float) -> np.ndarray:
    """Rodrigues rotation matrix about the unit ``axis`` (x @ R.T rotates the rows of x)."""

    x, y, z = axis
    k = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    t = np.radians(degrees)
    return np.eye(3) + np.sin(t) * k + (1.0 - np.cos(t)) * (k @ k)


def neighbours(bonded: Collection[tuple[int, int]]) -> dict[int, set[int]]:
    """Bonded partners of every atom that has one."""

    out: dict[int, set[int]] = {}
    for i, j in bonded:
        out.setdefault(i, set()).add(j)
        out.setdefault(j, set()).add(i)
    return out


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

    partners = neighbours(bonded)
    chains = [(i, j, k, m) for j, k in bonded for i in partners[j] - {k}
              for m in partners[k] - {i, j}]
    if not chains:
        return ()
    change = [abs((dihedral_deg(b, c) - dihedral_deg(a, c) + 180.0) % 360.0 - 180.0)
              for c in chains]
    return (CoordinateTerm(kind="dihedral", atoms=chains[int(np.argmax(change))]),)

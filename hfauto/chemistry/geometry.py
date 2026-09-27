"""Small geometry primitives shared by chemistry workflows."""

from __future__ import annotations

from collections.abc import Collection, Sequence

import numpy as np


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


def most_changed_dihedral(bonded: Collection[tuple[int, int]], a: np.ndarray, b: np.ndarray
                          ) -> tuple[int, int, int, int] | None:
    """The dihedral along bonds i-j, j-k, k-l whose angle differs most (periodically) between
    the structures ``a`` and ``b``; None when the bond graph has no such chain."""

    neighbours: dict[int, set[int]] = {}
    for i, j in bonded:
        neighbours.setdefault(i, set()).add(j)
        neighbours.setdefault(j, set()).add(i)
    chains = [(i, j, k, m) for j, k in bonded for i in neighbours[j] - {k}
              for m in neighbours[k] - {i, j}]
    change = [abs((dihedral_deg(b, c) - dihedral_deg(a, c) + 180.0) % 360.0 - 180.0)
              for c in chains]
    return chains[int(np.argmax(change))] if chains else None

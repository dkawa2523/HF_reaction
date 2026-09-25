"""Small geometry primitives shared by chemistry workflows."""

from __future__ import annotations

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

"""Small geometry primitives shared by chemistry workflows."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np

from hfauto.chemistry.xyz import XYZ, read_xyz


def align_coordinates(
    reference: np.ndarray,
    moving: np.ndarray,
) -> np.ndarray:
    """Rigidly align ``moving`` coordinates onto an ordered reference."""

    fixed = np.asarray(reference, dtype=float)
    return align_coordinates_by_indices(
        fixed, moving, range(len(fixed))
    )


def align_coordinates_by_indices(
    reference: np.ndarray,
    moving: np.ndarray,
    atom_indices: Sequence[int],
) -> np.ndarray:
    """Align a full geometry using a stable, atom-mapped reference subset."""

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
    indices = np.asarray(list(atom_indices), dtype=int)
    if (
        len(indices) == 0
        or np.any(indices < 0)
        or np.any(indices >= len(fixed))
        or len(set(indices.tolist())) != len(indices)
    ):
        raise ValueError("alignment atom indices are invalid")
    fixed_subset = fixed[indices]
    mobile_subset = mobile[indices]
    fixed_center = fixed_subset.mean(axis=0)
    mobile_center = mobile_subset.mean(axis=0)
    centered_fixed = fixed_subset - fixed_center
    centered_mobile = mobile_subset - mobile_center
    try:
        u, _singular_values, vt = np.linalg.svd(
            centered_mobile.T @ centered_fixed
        )
    except np.linalg.LinAlgError as exc:
        raise ValueError("Kabsch alignment failed") from exc
    rotation = u @ np.diag([1.0, 1.0, np.linalg.det(u @ vt)]) @ vt
    return (mobile - mobile_center) @ rotation + fixed_center


def kabsch_rmsd(
    coords_a: np.ndarray,
    coords_b: np.ndarray,
    atom_indices: Sequence[int] | None = None,
) -> float | None:
    """Return atom-order RMSD after optimal rigid alignment.

    ``None`` means the coordinate arrays or selected atom indices are not
    comparable.  A subset is useful for path continuity checks that should not
    be dominated by mobile hydrogens.
    """

    a = np.asarray(coords_a, dtype=float)
    b = np.asarray(coords_b, dtype=float)
    if a.shape != b.shape or a.ndim != 2 or a.shape[1] != 3 or len(a) == 0:
        return None
    if atom_indices is not None:
        indices = np.asarray(list(atom_indices), dtype=int)
        if len(indices) == 0 or np.any(indices < 0) or np.any(indices >= len(a)):
            return None
        a = a[indices]
        b = b[indices]

    try:
        aligned_a = align_coordinates(b, a)
    except ValueError:
        centered_a = a - a.mean(axis=0)
        centered_b = b - b.mean(axis=0)
        return float(
            np.sqrt(np.mean(np.sum((centered_a - centered_b) ** 2, axis=1)))
        )
    difference = aligned_a - b
    return float(np.sqrt(np.mean(np.sum(difference * difference, axis=1))))


def xyz_rmsd(
    xyz_a: XYZ,
    xyz_b: XYZ,
    atom_indices: Sequence[int] | None = None,
) -> float | None:
    """Return aligned RMSD when two XYZ objects have the same atom order."""

    if list(xyz_a.symbols) != list(xyz_b.symbols):
        return None
    return kabsch_rmsd(xyz_a.coords, xyz_b.coords, atom_indices)


def xyz_file_rmsd(
    path_a: str | Path,
    path_b: str | Path,
    atom_indices: Sequence[int] | None = None,
) -> float | None:
    """Read two XYZ files and return their aligned RMSD."""

    try:
        return xyz_rmsd(read_xyz(path_a), read_xyz(path_b), atom_indices)
    except (OSError, TypeError, ValueError):
        return None

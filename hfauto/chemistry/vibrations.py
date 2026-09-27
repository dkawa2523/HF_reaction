"""Harmonic analysis in the complement of translations and rotations (design §5.5).

Every engine's Hessian goes through here, so imaginary modes are counted one way.
Units: Hessians in Eh/bohr², coordinates in Å, masses in amu, frequencies in cm⁻¹
(imaginary values negative).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from pathlib import Path

import numpy as np

from hfauto.chemistry.elements import mass
from hfauto.core.constants import AMU_TO_ME, CM1_TO_HARTREE

# A singular value of the translation/rotation block counts when it exceeds this fraction
# of the largest one: 1e-4 Å noise on an optimized linear chain stays linear, 179° does not.
_EXTERNAL_RANK_TOL = 1.0e-3
# sqrt(Eh / (bohr² amu)) -> cm⁻¹ (CODATA 2018 through hfauto.core.constants).
_CM1_PER_SQRT_EIGENVALUE = 1.0 / (math.sqrt(AMU_TO_ME) * CM1_TO_HARTREE)
_PLANCK_J_S = 6.62607015e-34  # CODATA 2018, exact
_AMU_KG = 1.66053906660e-27  # CODATA 2018
_GHZ_AMU_A2 = _PLANCK_J_S / (8.0 * math.pi**2 * _AMU_KG * 1.0e-20) / 1.0e9


def _centered(symbols: Sequence[str], coords: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    masses = np.array([mass(symbol) for symbol in symbols])
    x = np.asarray(coords, dtype=float).reshape(-1, 3)
    if len(x) != len(masses):
        raise ValueError("symbols and coordinates differ in atom count")
    return x - masses @ x / masses.sum(), masses


def _square_hessian(hessian: np.ndarray, n_atoms: int) -> np.ndarray:
    h = np.asarray(hessian, dtype=float)
    if h.shape != (3 * n_atoms, 3 * n_atoms) or not np.isfinite(h).all():
        raise ValueError(f"Hessian must be a finite ({3 * n_atoms}, {3 * n_atoms}) array")
    return 0.5 * (h + h.T)


def _external_svd(symbols: Sequence[str], coords: np.ndarray) -> tuple[np.ndarray, int, np.ndarray]:
    """Full left singular basis of the mass-weighted translations/rotations and their rank k."""

    x, masses = _centered(symbols, coords)
    root = np.sqrt(masses)[:, None]
    axes = np.eye(3)
    vectors = [(root * axis).ravel() for axis in axes]
    vectors += [(root * np.cross(axis, x)).ravel() for axis in axes]
    u, singular, _ = np.linalg.svd(np.array(vectors).T, full_matrices=True)
    k = int(np.count_nonzero(singular > _EXTERNAL_RANK_TOL * singular.max()))
    return u, k, masses


def external_basis(symbols: Sequence[str], coords: np.ndarray) -> np.ndarray:
    """Orthonormal mass-weighted basis (3N, k) of translations and rotations; k = 5 if linear."""

    u, k, _ = _external_svd(symbols, coords)
    return u[:, :k]


def projected_frequencies(
    hessian_eh_bohr2: np.ndarray, symbols: Sequence[str], coords_A: np.ndarray
) -> tuple[np.ndarray, np.ndarray, int]:
    """Diagonalize Bᵀ H_mw B in the (3N − k)-dim complement of the external space.

    Returns ascending signed frequencies (cm⁻¹), Cartesian modes as unit rows (3N − k, 3N)
    with the mass weighting removed, and n_external = k.
    """

    u, k, masses = _external_svd(symbols, coords_A)
    h = _square_hessian(hessian_eh_bohr2, len(masses))
    inv_root = np.repeat(1.0 / np.sqrt(masses), 3)
    complement = u[:, k:]
    eigenvalues, vectors = np.linalg.eigh(
        complement.T @ (h * np.outer(inv_root, inv_root)) @ complement
    )
    freqs = np.sign(eigenvalues) * np.sqrt(np.abs(eigenvalues)) * _CM1_PER_SQRT_EIGENVALUE
    modes = (complement @ vectors).T * inv_root
    modes /= np.linalg.norm(modes, axis=1, keepdims=True)
    return freqs, modes, k


def shape_hessian(hessian_eh_bohr2: np.ndarray, coords_A: np.ndarray, reaction_mode: np.ndarray,
                  floor: float = 1.0e-3) -> np.ndarray:
    """Initial saddle-search Hessian whose only negative curvature is the reaction mode: P·H·P
    (P projects out the unweighted rigid motions) rebuilt from max(|λ_i|, floor) Eh/bohr²,
    negative for the eigenvector most parallel to ``reaction_mode``. Eigenvector following
    climbs that mode (moddir 1) and descends all others (Baker 1986); the inertia survives
    the internal-coordinate transformation (Sylvester)."""
    n = np.size(coords_A) // 3
    u, k, _ = _external_svd(["H"] * n, coords_A)  # equal masses: unweighted, about the centroid
    values, vectors = np.linalg.eigh(u[:, k:].T @ _square_hessian(hessian_eh_bohr2, n) @ u[:, k:])
    modes = u[:, k:] @ vectors
    curvature = np.maximum(np.abs(values), floor)
    curvature[np.argmax(np.abs(modes.T @ np.ravel(reaction_mode)))] *= -1.0
    return (modes * curvature) @ modes.T


def rotational_constants_ghz(symbols: Sequence[str], coords: np.ndarray) -> tuple[float, ...]:
    """(A, B, C) in GHz from ascending principal moments; 0.0 marks a vanishing moment.

    A moment vanishes by the same criterion that makes external_basis drop a rotation,
    so a linear molecule gives (0.0, B, B) as GoodVibes expects.
    """

    x, masses = _centered(symbols, coords)
    inertia = np.einsum("i,ij,ik->jk", masses, x, x)
    moments = np.linalg.eigvalsh(np.trace(inertia) * np.eye(3) - inertia)
    floor = _EXTERNAL_RANK_TOL**2 * max(masses.sum(), float(moments.max()))
    return tuple(0.0 if m <= floor else _GHZ_AMU_A2 / float(m) for m in moments)


def to_canonical_npy(hessian: np.ndarray, path: str | Path) -> Path:
    """Write the canonical form: symmetric (3N, 3N) float64 in Eh/bohr², input frame."""

    h = np.asarray(hessian, dtype=float)
    if h.ndim != 2 or h.shape[0] % 3:
        raise ValueError("Hessian must be a (3N, 3N) array")
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("wb") as handle:
        np.save(handle, _square_hessian(h, h.shape[0] // 3))
    return target

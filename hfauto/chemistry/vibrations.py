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
from hfauto.core.constants import AMU_TO_ME, BOHR_TO_ANGSTROM, CM1_TO_HARTREE

# Without a point group (an engine reading its own freq, trial directions, shaped Hessians), a
# singular value of the translation/rotation block counts when it exceeds this fraction of the
# largest one: 1e-4 Å noise on an optimized linear chain stays linear, 179° does not.
_EXTERNAL_RANK_TOL = 1.0e-3
# sqrt(Eh / (bohr² amu)) -> cm⁻¹ (CODATA 2018 through hfauto.core.constants).
_CM1_PER_SQRT_EIGENVALUE = 1.0 / (math.sqrt(AMU_TO_ME) * CM1_TO_HARTREE)
_PLANCK_J_S = 6.62607015e-34  # CODATA 2018, exact
_AMU_KG = 1.66053906660e-27  # CODATA 2018
_GHZ_AMU_A2 = _PLANCK_J_S / (8.0 * math.pi**2 * _AMU_KG * 1.0e-20) / 1.0e9
_CURVATURE_FLOOR = 1.0e-3  # Eh/bohr², the least curvature of a shaped Hessian's internal mode
KAPPA_MIN = 0.05  # Eh/bohr², the least negative curvature along a saddle search's direction


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


def _external_svd(symbols: Sequence[str], coords: np.ndarray, linear: bool | None = None
                  ) -> tuple[np.ndarray, int, np.ndarray]:
    """Full left singular basis of the mass-weighted translations/rotations and their count k:
    3 for an atom, else 5 if ``linear`` and 6 if not; without ``linear``, the numerical rank."""

    x, masses = _centered(symbols, coords)
    root = np.sqrt(masses)[:, None]
    axes = np.eye(3)
    vectors = [(root * axis).ravel() for axis in axes]
    vectors += [(root * np.cross(axis, x)).ravel() for axis in axes]
    u, singular, _ = np.linalg.svd(np.array(vectors).T, full_matrices=True)
    if linear is None:
        return u, int(np.count_nonzero(singular > _EXTERNAL_RANK_TOL * singular.max())), masses
    return u, 3 if len(masses) == 1 else 5 if linear else 6, masses


def external_basis(symbols: Sequence[str], coords: np.ndarray) -> np.ndarray:
    """Orthonormal mass-weighted basis (3N, k) of translations and rotations; k = 5 if linear."""

    u, k, _ = _external_svd(symbols, coords)
    return u[:, :k]


def _mass_weighted_internal(hessian_eh_bohr2: np.ndarray, symbols: Sequence[str],
                            coords_A: np.ndarray, linear: bool | None = None
                            ) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Eigenvalues and mass-weighted eigenvectors (3N, 3N − k) of Bᵀ H_mw B in the complement B
    of the external space, the weights 1/√m (3N) and k."""
    u, k, masses = _external_svd(symbols, coords_A, linear)
    h = _square_hessian(hessian_eh_bohr2, len(masses))
    inv_root = np.repeat(1.0 / np.sqrt(masses), 3)
    complement = u[:, k:]
    values, vectors = np.linalg.eigh(
        complement.T @ (h * np.outer(inv_root, inv_root)) @ complement
    )
    return values, complement @ vectors, inv_root, k


def projected_frequencies(
    hessian_eh_bohr2: np.ndarray, symbols: Sequence[str], coords_A: np.ndarray,
    linear: bool | None = None,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Diagonalize Bᵀ H_mw B in the (3N − k)-dim complement of the external space, k from
    ``linear`` (a point group's verdict, chemistry.symmetry) or else the numerical rank.

    Returns ascending signed frequencies (cm⁻¹), Cartesian modes as unit rows (3N − k, 3N)
    with the mass weighting removed, and n_external = k.
    """

    values, vectors, inv_root, k = _mass_weighted_internal(
        hessian_eh_bohr2, symbols, coords_A, linear)
    freqs = np.sign(values) * np.sqrt(np.abs(values)) * _CM1_PER_SQRT_EIGENVALUE
    modes = vectors.T * inv_root
    modes /= np.linalg.norm(modes, axis=1, keepdims=True)
    return freqs, modes, k


def stationarity_gap(hessian_eh_bohr2: np.ndarray, gradient_eh_bohr: Sequence[float],
                     symbols: Sequence[str], coords_A: np.ndarray) -> float:
    """ΔE_N = Σ g_i²/(2|λ_i|) (Eh) over the modes of projected_frequencies: the energy one
    Newton step of the quadratic model from the gradient and Hessian at ``coords`` moves by. A
    point is stationary when it is within an identity basin's energy (identity.
    BASIN_DE_HARTREE): an optimization converged on a flat surface can stop on a shoulder
    whose frequencies show nothing. |λ| is not floored: a floor (shape_hessian's) caps the soft
    modes' terms, which carry the shoulder."""
    values, vectors, inv_root, _ = _mass_weighted_internal(hessian_eh_bohr2, symbols, coords_A)
    g = vectors.T @ (np.ravel(gradient_eh_bohr) * inv_root)
    return float(np.sum(g**2 / (2.0 * np.abs(values))))


def _internal_eigen(hessian_eh_bohr2: np.ndarray, coords_A: np.ndarray
                    ) -> tuple[np.ndarray, np.ndarray]:
    """Eigenvalues and Cartesian eigenvectors (3N, 3N − k) of P_r·H·P_r on the internal
    motions, P_r projecting out the unweighted rigid motions (equal masses, about the
    centroid): the modes of every driver model."""
    n = np.size(coords_A) // 3
    u, k, _ = _external_svd(["H"] * n, coords_A)
    internal = u[:, k:]
    values, vectors = np.linalg.eigh(internal.T @ _square_hessian(hessian_eh_bohr2, n) @ internal)
    return values, internal @ vectors


def shape_hessian(hessian_eh_bohr2: np.ndarray, coords_A: np.ndarray,
                  direction: np.ndarray | None = None) -> np.ndarray:
    """Initial driver Hessian. H₊ = P_r·H·P_r on its eigenvectors with max(|λ|,
    _CURVATURE_FLOOR), rigid motions at 0: positive definite on the internal motions (a
    minimization from a TS's side). With a ``direction`` d (3N; d̂ its internal part,
    normalized): P·H₊·P − κ·d̂d̂ᵀ with P = I − d̂d̂ᵀ and κ = max(d̂ᵀH₊d̂, KAPPA_MIN), so d̂ is the
    only negative curvature, whatever the eigenvectors of H near it: eigenvector following
    climbs it (moddir 1) and descends all others (Baker 1986). The inertia survives the
    internal-coordinate transformation (Sylvester)."""
    values, modes = _internal_eigen(hessian_eh_bohr2, coords_A)
    positive = (modes * np.maximum(np.abs(values), _CURVATURE_FLOOR)) @ modes.T
    if direction is None:
        return positive
    d = modes @ (modes.T @ np.ravel(direction))
    norm = float(np.linalg.norm(d))
    if norm <= 1.0e-8 * float(np.linalg.norm(direction)):
        raise ValueError("the direction has no internal component")
    d /= norm
    kappa = max(float(d @ positive @ d), KAPPA_MIN)
    p = np.eye(d.size) - np.outer(d, d)
    return p @ positive @ p - kappa * np.outer(d, d)


def newton_step(hessian_eh_bohr2: np.ndarray, gradient_eh_bohr: Sequence[float],
                coords_A: np.ndarray, *, signed: bool = False) -> np.ndarray:
    """One Newton step −Σ g_i/λ'_i v_i (N, 3; Å) on shape_hessian's internal modes, |λ'| =
    max(|λ|, _CURVATURE_FLOOR): with λ' > 0, −H₊⁺g, downhill along every mode, imaginary ones
    included (a minimum's, as the optimization from H₊ would start); ``signed``, λ' with λ's
    sign, to the quadratic model's stationary point (a saddle's)."""
    values, modes = _internal_eigen(hessian_eh_bohr2, coords_A)
    curvature = np.maximum(np.abs(values), _CURVATURE_FLOOR)
    if signed:
        curvature = np.where(values < 0.0, -curvature, curvature)
    step = -modes @ ((modes.T @ np.ravel(gradient_eh_bohr)) / curvature)
    return step.reshape(-1, 3) * BOHR_TO_ANGSTROM


def rotational_constants_ghz(symbols: Sequence[str], coords: np.ndarray, linear: bool
                             ) -> tuple[float, ...]:
    """(A, B, C) in GHz from the ascending principal moments of a structure of two or more
    atoms (a point group's symmetrized structure); a linear one gives (0.0, B, B) as GoodVibes
    expects."""

    x, masses = _centered(symbols, coords)
    inertia = np.einsum("i,ij,ik->jk", masses, x, x)
    moments = np.linalg.eigvalsh(np.trace(inertia) * np.eye(3) - inertia)
    return tuple(0.0 if linear and i == 0 else _GHZ_AMU_A2 / float(m)
                 for i, m in enumerate(moments))


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

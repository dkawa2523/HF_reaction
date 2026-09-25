"""Small, auditable Cartesian mode-following geometry operations."""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from hfauto.chemistry.xyz import XYZ
from hfauto.core.hashing import fingerprint_dict

MODE_FOLLOWING_SCHEMA_VERSION = "hfauto.minimum_mode_following.v1"
DEFAULT_MAXIMUM_ATOM_DISPLACEMENT_A = 0.10
MAXIMUM_ALLOWED_ATOM_DISPLACEMENT_A = 0.25


def _finite_cartesian_array(value: Any, shape: tuple[int, int], label: str) -> np.ndarray:
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be a finite {shape[0]}x3 array") from exc
    if array.shape != shape:
        raise ValueError(f"{label} must have shape {shape}, found {array.shape}")
    if not np.isfinite(array).all():
        raise ValueError(f"{label} contains a non-finite value")
    return array


def build_cartesian_mode_following_seeds(
    source: XYZ,
    cartesian_mode: Any,
    *,
    maximum_atom_displacement_A: float = DEFAULT_MAXIMUM_ATOM_DISPLACEMENT_A,
) -> dict[str, tuple[XYZ, dict[str, Any]]]:
    """Create finite ``+`` and ``-`` seeds along one NWChem Cartesian mode.

    NWChem's projected table is explicitly in Cartesian coordinates: the
    mass-weighted Hessian eigenvector has already been transformed back by the
    program.  Consequently this helper preserves all printed component ratios
    and applies one scalar only; it never divides components by atomic masses a
    second time.  The scalar is chosen so the largest per-atom Cartesian shift
    is exactly ``maximum_atom_displacement_A``.

    The returned geometries are in Angstrom because ``XYZ`` coordinates and
    the requested displacement are both in Angstrom.  No calculation is run,
    and the two outputs remain descendants of one source seed rather than new
    independent sampling seeds.
    """

    symbols = [str(symbol) for symbol in source.symbols]
    atom_count = len(symbols)
    if atom_count <= 0:
        raise ValueError("source geometry must contain at least one atom")
    source_coords = _finite_cartesian_array(
        source.coords,
        (atom_count, 3),
        "source coordinates",
    )
    mode = _finite_cartesian_array(
        cartesian_mode,
        (atom_count, 3),
        "Cartesian mode",
    )
    if isinstance(maximum_atom_displacement_A, bool) or not isinstance(
        maximum_atom_displacement_A, (int, float)
    ):
        raise TypeError("maximum_atom_displacement_A must be a finite number")
    requested = float(maximum_atom_displacement_A)
    if not math.isfinite(requested) or not (
        0.0 < requested <= MAXIMUM_ALLOWED_ATOM_DISPLACEMENT_A
    ):
        raise ValueError(
            "maximum_atom_displacement_A must be greater than zero and no more "
            f"than {MAXIMUM_ALLOWED_ATOM_DISPLACEMENT_A:.2f} A"
        )

    atom_norms = np.linalg.norm(mode, axis=1)
    raw_maximum = float(np.max(atom_norms))
    if not math.isfinite(raw_maximum) or raw_maximum <= 1.0e-14:
        raise ValueError("Cartesian mode has zero or unusably small norm")
    scale_A_sqrt_amu = requested / raw_maximum
    displacement = mode * scale_A_sqrt_amu
    actual_atom_norms = np.linalg.norm(displacement, axis=1)
    actual_maximum = float(np.max(actual_atom_norms))
    if not np.isclose(actual_maximum, requested, rtol=1.0e-12, atol=1.0e-14):
        raise ValueError("Cartesian mode normalization did not reach the requested shift")

    common = {
        "schema_version": MODE_FOLLOWING_SCHEMA_VERSION,
        "atom_count": atom_count,
        "atom_symbols_in_order": symbols,
        "atom_order_fingerprint": fingerprint_dict({"symbols": symbols}),
        "source_coordinate_units": "angstrom",
        "raw_mode_coordinate_convention": (
            "NWChem projected normal-mode eigenvector in Cartesian coordinates"
        ),
        "raw_mode_component_units": "amu^-1/2",
        "mass_weighting_handling": (
            "NWChem already converted the mass-weighted Hessian eigenvector to "
            "Cartesian components; no second 1/sqrt(mass) factor was applied"
        ),
        "normalization": "largest per-atom Cartesian vector norm",
        "requested_maximum_atom_displacement_A": requested,
        "applied_maximum_atom_displacement_A": actual_maximum,
        "normalization_scale_A_sqrt_amu": scale_A_sqrt_amu,
        "rms_atom_displacement_A": float(
            math.sqrt(float(np.mean(np.square(actual_atom_norms))))
        ),
        "raw_mode_fingerprint": fingerprint_dict(
            {"cartesian_mode_amu^-1/2": mode.tolist()}
        ),
        "normalized_displacement_fingerprint": fingerprint_dict(
            {"cartesian_displacement_A": displacement.tolist()}
        ),
        "independent_seed_count_increased": False,
    }

    results: dict[str, tuple[XYZ, dict[str, Any]]] = {}
    for direction, sign in (("plus", 1), ("minus", -1)):
        coords = source_coords + sign * displacement
        if not np.isfinite(coords).all():
            raise ValueError(f"{direction} displaced geometry contains a non-finite value")
        geometry = XYZ(
            symbols=list(symbols),
            coords=coords,
            comment=(
                f"minimum mode-following {direction}; "
                f"max displacement={actual_maximum:.8f} A"
            ),
        )
        provenance = {
            **common,
            "direction": direction,
            "direction_sign": sign,
            "displaced_geometry_fingerprint": fingerprint_dict(
                {
                    "symbols": symbols,
                    "coords_A": coords.tolist(),
                }
            ),
        }
        results[direction] = (geometry, provenance)
    return results

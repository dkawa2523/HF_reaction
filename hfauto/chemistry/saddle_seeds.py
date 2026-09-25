"""Backend-neutral interior seeds and one-dimensional saddle brackets."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from hfauto.chemistry.geometry import align_coordinates
from hfauto.chemistry.xyz import XYZ, read_xyz, write_xyz


def endpoint_biased_fractions(
    values: list[float] | tuple[float, ...] | None = None,
) -> list[float]:
    """Return a bounded, ordered scan grid dense near the starting basin."""

    requested = values or (0.0, 0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.5, 0.75, 1.0)
    fractions = sorted({float(value) for value in requested})
    if len(fractions) < 3 or fractions[0] != 0.0 or fractions[-1] != 1.0:
        raise ValueError("saddle bracket fractions must include 0 and 1")
    if any(value < 0.0 or value > 1.0 for value in fractions):
        raise ValueError("saddle bracket fractions must stay inside [0, 1]")
    return fractions


def aligned_interpolation_xyz(
    start_xyz: str | Path,
    target_xyz: str | Path,
    out_xyz: str | Path,
    fraction: float,
) -> Path:
    """Interpolate atom-mapped structures after a proper Kabsch alignment."""

    value = float(fraction)
    if not 0.0 <= value <= 1.0:
        raise ValueError("interpolation fraction must be in [0, 1]")
    start = read_xyz(start_xyz)
    target = read_xyz(target_xyz)
    if start.symbols != target.symbols:
        raise ValueError("atom symbols/order differ between saddle endpoints")
    aligned_target = align_coordinates(start.coords, target.coords)
    coords = (1.0 - value) * start.coords + value * aligned_target
    return write_xyz(
        XYZ(
            list(start.symbols),
            coords,
            comment=f"state=saddle_seed aligned_fraction={value:.8f}",
        ),
        out_xyz,
    )


def bracket_energy_maximum(
    points: list[dict[str, Any]],
    *,
    minimum_prominence_hartree: float = 1.0e-6,
) -> dict[str, Any]:
    """Find a strict internal energy maximum; never call an endpoint a saddle."""

    ordered = sorted(points, key=lambda point: float(point["fraction"]))
    if len(ordered) < 3:
        return {
            "classification": "insufficient_profile",
            "candidate": None,
            "internal_maxima": [],
        }
    maxima: list[dict[str, Any]] = []
    threshold = float(minimum_prominence_hartree)
    for index in range(1, len(ordered) - 1):
        previous, current, following = ordered[index - 1 : index + 2]
        x0, x1, x2 = [float(point["fraction"]) for point in (previous, current, following)]
        e0, e1, e2 = [float(point["energy_hartree"]) for point in (previous, current, following)]
        prominence = min(e1 - e0, e1 - e2)
        if prominence < threshold:
            continue
        curvature = 2.0 * (
            e0 / ((x0 - x1) * (x0 - x2))
            + e1 / ((x1 - x0) * (x1 - x2))
            + e2 / ((x2 - x0) * (x2 - x1))
        )
        if curvature >= 0.0:
            continue
        maxima.append(
            {
                **current,
                "profile_index": index,
                "prominence_hartree": prominence,
                "directional_curvature_hartree": curvature,
                "bracket_fractions": [x0, x1, x2],
            }
        )
    candidate = maxima[0] if len(maxima) == 1 else None
    return {
        "classification": (
            "resolved_internal_maximum"
            if candidate
            else "multiple_internal_maxima"
            if maxima
            else "no_bracketed_maximum"
        ),
        "candidate": candidate,
        "internal_maxima": maxima,
    }

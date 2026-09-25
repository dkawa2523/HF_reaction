"""Read XYZ trajectories without depending on a calculation backend."""

from __future__ import annotations

from itertools import pairwise
from pathlib import Path

import numpy as np

from hfauto.chemistry.geometry import align_coordinates
from hfauto.chemistry.xyz import XYZ


def reaction_mode_reference(
    images: list[XYZ],
    image_index: int,
    *,
    source: str,
) -> dict[str, object]:
    """Build a rigid-motion-free local path tangent for mode following."""

    index = int(image_index)
    if len(images) < 3 or not 0 < index < len(images) - 1:
        raise ValueError("path tangent requires an internal image")
    symbols = list(images[index].symbols)
    if any(list(image.symbols) != symbols for image in images):
        raise ValueError("path tangent atom symbols/order are inconsistent")
    center = np.asarray(images[index].coords, dtype=float)
    previous = align_coordinates(center, images[index - 1].coords)
    following = align_coordinates(center, images[index + 1].coords)
    tangent = following - previous
    norm = float(np.linalg.norm(tangent))
    if not np.isfinite(norm) or norm <= 1.0e-12:
        raise ValueError("path tangent is zero or non-finite")
    return {
        "kind": "local_path_tangent",
        "source": str(source),
        "image_index": index,
        "component_units": "angstrom",
        "cartesian_displacements": (tangent / norm).tolist(),
    }


def endpoint_mode_reference(
    start: XYZ,
    end: XYZ,
    *,
    source: str = "aligned_endpoint_displacement",
) -> dict[str, object]:
    """Build a coarse path tangent when no converged local path is available."""

    if list(start.symbols) != list(end.symbols):
        raise ValueError("endpoint tangent atom symbols/order differ")
    aligned_end = align_coordinates(start.coords, end.coords)
    tangent = aligned_end - np.asarray(start.coords, dtype=float)
    norm = float(np.linalg.norm(tangent))
    if not np.isfinite(norm) or norm <= 1.0e-12:
        raise ValueError("endpoint tangent is zero or non-finite")
    return {
        "kind": "endpoint_displacement",
        "source": str(source),
        "component_units": "angstrom",
        "cartesian_displacements": (tangent / norm).tolist(),
    }


def read_xyz_trajectory(path: str | Path) -> list[XYZ]:
    """Read concatenated XYZ blocks and fail on malformed or partial data."""

    source = Path(path)
    lines = source.read_text(encoding="utf-8", errors="strict").splitlines()
    images: list[XYZ] = []
    cursor = 0
    while cursor < len(lines):
        if not lines[cursor].strip():
            cursor += 1
            continue
        try:
            atom_count = int(lines[cursor].strip())
        except ValueError as exc:
            raise ValueError(
                f"invalid XYZ atom count at line {cursor + 1} in {source}"
            ) from exc
        if atom_count <= 0 or cursor + atom_count + 2 > len(lines):
            raise ValueError(f"incomplete XYZ block at line {cursor + 1} in {source}")
        comment = lines[cursor + 1]
        symbols: list[str] = []
        coordinates: list[list[float]] = []
        for line_number, line in enumerate(
            lines[cursor + 2 : cursor + atom_count + 2], start=cursor + 3
        ):
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(f"invalid XYZ row at line {line_number} in {source}")
            try:
                coordinates.append([float(value) for value in fields[1:4]])
            except ValueError as exc:
                raise ValueError(
                    f"invalid XYZ coordinate at line {line_number} in {source}"
                ) from exc
            symbols.append(fields[0])
        images.append(
            XYZ(symbols=symbols, coords=np.asarray(coordinates, dtype=float), comment=comment)
        )
        cursor += atom_count + 2
    if not images:
        raise ValueError(f"XYZ trajectory contains no images: {source}")
    return images


def write_xyz_trajectory(images: list[XYZ], path: str | Path) -> Path:
    """Write a complete multi-XYZ path without changing image coordinates."""

    if len(images) < 2:
        raise ValueError("XYZ trajectory must contain at least two images")
    symbols = images[0].symbols
    if any(image.symbols != symbols for image in images):
        raise ValueError("XYZ trajectory atom symbols/order are inconsistent")
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    blocks: list[str] = []
    for image in images:
        blocks.extend(
            [
                str(len(image.symbols)),
                image.comment or "generated_by=hfauto",
                *(
                    f"{symbol:2s} {x: .8f} {y: .8f} {z: .8f}"
                    for symbol, (x, y, z) in zip(
                        image.symbols, image.coords
                    )
                ),
            ]
        )
    output.write_text("\n".join(blocks) + "\n", encoding="utf-8")
    return output


def normalized_xyz_arc_lengths(images: list[XYZ]) -> list[float]:
    """Return a rigid-motion-free normalized coordinate for an XYZ path.

    Different discretizations need not place image ``i`` at the same physical
    point.  Profile comparison therefore uses accumulated aligned Cartesian
    arc length instead of the image index.  Atom correspondence remains an
    explicit precondition; this helper never guesses permutations.
    """

    if len(images) < 2:
        raise ValueError("XYZ trajectory must contain at least two images")
    symbols = list(images[0].symbols)
    if any(list(image.symbols) != symbols for image in images):
        raise ValueError("XYZ trajectory atom symbols/order are inconsistent")
    segment_lengths: list[float] = []
    for previous, following in pairwise(images):
        aligned = align_coordinates(previous.coords, following.coords)
        length = float(
            np.linalg.norm(
                aligned - np.asarray(previous.coords, dtype=float)
            )
        )
        if not np.isfinite(length):
            raise ValueError("XYZ trajectory has a non-finite segment")
        segment_lengths.append(length)
    cumulative = np.concatenate(([0.0], np.cumsum(segment_lengths)))
    if cumulative[-1] <= 1.0e-12:
        raise ValueError("XYZ trajectory has zero aligned arc length")
    return [float(value / cumulative[-1]) for value in cumulative]


def resample_xyz_trajectory(images: list[XYZ], image_count: int) -> list[XYZ]:
    """Rediscretize one ordered path at uniform Cartesian arc length.

    Endpoints are preserved exactly.  This changes only path resolution; it
    does not add reaction-specific constraints or alter atom correspondence.
    """

    if len(images) < 2:
        raise ValueError("XYZ trajectory must contain at least two images")
    count = int(image_count)
    if count < 3:
        raise ValueError("resampled path requires at least three images")
    symbols = list(images[0].symbols)
    if any(list(image.symbols) != symbols for image in images):
        raise ValueError("XYZ trajectory atom symbols/order are inconsistent")
    coordinates = np.stack(
        [np.asarray(image.coords, dtype=float) for image in images]
    )
    segment_lengths = np.linalg.norm(
        np.diff(coordinates, axis=0).reshape(len(images) - 1, -1), axis=1
    )
    cumulative = np.concatenate(([0.0], np.cumsum(segment_lengths)))
    if not np.isfinite(cumulative).all() or cumulative[-1] <= 1.0e-12:
        raise ValueError("cannot resample a zero-length or non-finite path")
    keep = np.concatenate(([True], np.diff(cumulative) > 1.0e-12))
    cumulative = cumulative[keep]
    coordinates = coordinates[keep]
    targets = np.linspace(0.0, cumulative[-1], count)
    flat = coordinates.reshape(len(coordinates), -1)
    sampled = np.column_stack(
        [np.interp(targets, cumulative, flat[:, column]) for column in range(flat.shape[1])]
    ).reshape(count, len(symbols), 3)
    sampled[0] = images[0].coords
    sampled[-1] = images[-1].coords
    return [
        XYZ(
            symbols=list(symbols),
            coords=coords,
            comment=f"resampled_path image={index + 1}/{count}",
        )
        for index, coords in enumerate(sampled)
    ]


# Compatibility name used by CREST callers and public tests.

"""XYZ structures: the model, file I/O (one structure or a trajectory), fingerprints and
formulas."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from hfauto.core.evidence import FileRef, Geometry


@dataclass
class XYZ:
    symbols: list[str]
    coords: np.ndarray
    comment: str = ""


def _block(xyz: XYZ) -> list[str]:
    rows = zip(xyz.symbols, xyz.coords)
    return [str(len(xyz.symbols)), xyz.comment or "generated_by=hfauto",
            *(f"{s:2s} {x: .8f} {y: .8f} {z: .8f}" for s, (x, y, z) in rows)]


def read_xyz_trajectory(path: str | Path) -> list[XYZ]:
    """Concatenated XYZ blocks; malformed or partial data is a ValueError."""

    source = Path(path)
    lines = source.read_text(encoding="utf-8").splitlines()
    images: list[XYZ] = []
    cursor = 0
    while cursor < len(lines):
        if not lines[cursor].strip():
            cursor += 1
            continue
        try:
            atom_count = int(lines[cursor].strip())
        except ValueError as exc:
            raise ValueError(f"invalid XYZ atom count at line {cursor + 1} in {source}") from exc
        if atom_count <= 0 or cursor + atom_count + 2 > len(lines):
            raise ValueError(f"incomplete XYZ block at line {cursor + 1} in {source}")
        symbols: list[str] = []
        coordinates: list[list[float]] = []
        for number, line in enumerate(lines[cursor + 2 : cursor + atom_count + 2], cursor + 3):
            fields = line.split()
            try:
                if len(fields) < 4:
                    raise ValueError(line)
                coordinates.append([float(value) for value in fields[1:4]])
            except ValueError as exc:
                raise ValueError(f"invalid XYZ row at line {number} in {source}") from exc
            symbols.append(fields[0])
        images.append(XYZ(symbols, np.asarray(coordinates, dtype=float), lines[cursor + 1]))
        cursor += atom_count + 2
    if not images:
        raise ValueError(f"XYZ file contains no structure: {source}")
    return images


def read_xyz(path: str | Path) -> XYZ:
    """The one structure of an XYZ file."""

    images = read_xyz_trajectory(path)
    if len(images) != 1:
        raise ValueError(f"{path}: {len(images)} structures, expected one")
    return images[0]


def written_geometry(path: str | Path, file_ref: Callable[[Path], FileRef]) -> Geometry:
    """The Geometry of an xyz file: its FileRef and the fingerprint of the coordinates as the
    file holds them (8 decimals)."""

    xyz = read_xyz(path)
    return Geometry(file=file_ref(Path(path)), symbols=tuple(xyz.symbols),
                    fingerprint=geometry_fingerprint(xyz.symbols, xyz.coords))


def write_xyz(xyz: XYZ, path: str | Path) -> Path:
    return _write(_block(xyz), path)


def write_xyz_trajectory(images: list[XYZ], path: str | Path) -> Path:
    """A multi-XYZ path of at least two images with one atom list."""

    if len(images) < 2:
        raise ValueError("XYZ trajectory must contain at least two images")
    if any(image.symbols != images[0].symbols for image in images):
        raise ValueError("XYZ trajectory atom symbols/order are inconsistent")
    return _write([line for image in images for line in _block(image)], path)


def _write(lines: list[str], path: str | Path) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def geometry_fingerprint(symbols: Sequence[str], coords: np.ndarray) -> str:
    """sha256 of the element list and the coordinates rounded to 1e-6 Å (normalized JSON)."""

    rounded = np.round(np.asarray(coords, dtype=float).reshape(-1, 3), 6) + 0.0  # drops -0.0
    if len(rounded) != len(symbols):
        raise ValueError("symbols and coordinates differ in atom count")
    payload = json.dumps(
        {"symbols": list(symbols), "coords": rounded.tolist()}, separators=(",", ":")
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def hill_formula(symbols: Sequence[str]) -> str:
    """Hill-order formula: C, then H, then the rest alphabetically (all alphabetical without C)."""

    counts = Counter(symbols)
    first = ["C", "H"] if "C" in counts else []
    order = first + sorted(element for element in counts if element not in first)
    return "".join(f"{e}{counts[e] if counts[e] > 1 else ''}" for e in order if e in counts)


def composition_key(symbols: Sequence[str], charge: int, multiplicity: int) -> str:
    return f"{hill_formula(symbols)}_q{int(charge)}_m{int(multiplicity)}"


@dataclass(frozen=True, eq=False)
class Molecule:
    xyz: XYZ
    charge: int
    multiplicity: int

    def fingerprint(self) -> str:
        """Geometry fingerprint extended with charge and multiplicity."""

        payload = json.dumps(
            [geometry_fingerprint(self.xyz.symbols, self.xyz.coords), self.charge,
             self.multiplicity],
            separators=(",", ":"),
        )
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    def write(self, path: str | Path) -> Path:
        return write_xyz(self.xyz, path)

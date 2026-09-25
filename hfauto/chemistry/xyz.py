"""Minimal XYZ geometry model and file I/O.

This module has no reaction-family assumptions.  Chemistry builders consume
it, while path, conformer, and quantum backends can use it directly.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class XYZ:
    symbols: list[str]
    coords: np.ndarray
    comment: str = ""


def read_xyz(path: str | Path) -> XYZ:
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValueError(f"Empty XYZ file: {path}")
    atom_count = int(lines[0].strip())
    comment = lines[1] if len(lines) > 1 else ""
    symbols: list[str] = []
    coords: list[list[float]] = []
    for line in lines[2 : 2 + atom_count]:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Invalid XYZ line in {path}: {line!r}")
        symbols.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    if len(symbols) != atom_count:
        raise ValueError(
            f"XYZ atom count mismatch in {path}: "
            f"header={atom_count}, parsed={len(symbols)}"
        )
    return XYZ(
        symbols=symbols,
        coords=np.asarray(coords, dtype=float),
        comment=comment,
    )


def write_xyz(xyz: XYZ, path: str | Path) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    lines = [str(len(xyz.symbols)), xyz.comment or "generated_by=hfauto"]
    for symbol, (x, y, z) in zip(xyz.symbols, xyz.coords):
        lines.append(f"{symbol:2s} {x: .8f} {y: .8f} {z: .8f}")
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def xyz_files_have_same_atom_order(path_a: str | Path, path_b: str | Path) -> bool:
    """Return whether two XYZ files contain the same ordered element list."""

    xyz_a = read_xyz(path_a)
    xyz_b = read_xyz(path_b)
    return xyz_a.symbols == xyz_b.symbols

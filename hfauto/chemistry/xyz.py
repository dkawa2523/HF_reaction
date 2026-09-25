"""Minimal XYZ geometry model and file I/O.

This module has no reaction-family assumptions.  Chemistry builders consume
it, while path, conformer, and quantum backends can use it directly.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from collections.abc import Sequence
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

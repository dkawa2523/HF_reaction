"""NWChem output parsing (design §5.1, §6.3): pure functions over the output text and the
files a job leaves in its permanent directory.

The Level is observed, never assumed: version, xc functional or CCSD(T), basis (``/cart``
when cartesian), DFT-D3 variant, grid, SCF energy tolerance, charge and multiplicity. Frequencies are never taken from the text for
decisions; the ``.hess`` file is converted to the canonical ``.npy`` instead. The last DFT
gradient block of a driver job is its final frame's (487 NWChem opt / saddle Evidence of the
validation runs, within 5e-7 Å).
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path

import numpy as np

from hfauto.chemistry.vibrations import to_canonical_npy
from hfauto.chemistry.xyz import read_xyz
from hfauto.core.constants import BOHR_TO_ANGSTROM
from hfauto.core.evidence import Failure, FailureKind, Level

NORMAL_END = "Total times  cpu:"
_PROJECTED = "(Projected Frequencies expressed in cm-1)"
_GEOMETRY = re.compile(r'^\s*Geometry "[^"]*" -> "[^"]*"', re.MULTILINE)
_NUMBER = r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?"
_ATOM_ROW = re.compile(
    rf"^\s*\d+\s+([A-Za-z]+)\S*\s+{_NUMBER}\s+({_NUMBER})\s+({_NUMBER})\s+({_NUMBER})\s*$"
)
_BASIS = re.compile(
    r'^ Summary of "ao basis" -> "ao basis" \((spherical|cartesian)\)\n(?:.*\n){3}'
    r"((?:[ \t]*\S.*\n)+)",
    re.MULTILINE,
)
# fatal only: "AUTOZ failed to generate good internal coordinates. Cartesian coordinates will
# be used" is a notice, after which NWChem goes on (HCN, an HF-dimer saddle at maxiter)
_AUTOZ = re.compile("insufficient internal variables|geom_binvr: #indep variables incorrect"
                    "|regeneration of autoz failed")
_SCF = re.compile(r"Calculation failed to converge|SCF not converged")
_GEOMETRY_MAXITER = re.compile("Failed to converge in maximum number of steps")
_GRADIENT_ROW = re.compile(rf"^\s*\d+\s+[A-Za-z]\S*((?:\s+{_NUMBER}){{6}})\s*$")
_DFT_ENERGY = r"Total DFT energy =\s+(\S+)"
# the ccsd module (RHF) or the TCE (UHF; ROHF before P0d)
_CCSD_T = r"(?:Total CCSD\(T\) energy:|CCSD\(T\) total energy / hartree\s+=)\s+(\S+)"
_UHF = re.compile(r"^\s*alpha electrons\s+=\s+(\d+)\n\s*beta\s+electrons\s+=\s+(\d+)", re.MULTILINE)


def _number(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


def _last(pattern: str, text: str) -> str | None:
    """First group of the last match (the final SCF / job wins)."""
    matches = re.findall(pattern, text, re.MULTILINE)
    return matches[-1] if matches else None


def version(text: str) -> str | None:
    return _last(r"Northwest Computational Chemistry Package \(NWChem\)\s+(\S+)", text)


def basis(text: str) -> str | None:
    """The basis of every tag in the last resolved summary; several names are comma-joined."""
    matches = _BASIS.findall(text)
    if not matches:
        return None
    kind, rows = matches[-1]
    names = sorted({row.split()[1].lower() for row in rows.splitlines() if len(row.split()) > 1})
    return ",".join(names) + ("/cart" if kind == "cartesian" else "")


def _dispersion(text: str, xc: str) -> str | None:
    if xc.endswith("-d3"):  # the D3 term belongs to the functional (wb97x-d3)
        return None
    if "DFT-D3BJ Model" in text:
        return "d3bj"
    return "d3zero" if "DFT-D3 Model" in text else None


def _dft_level(text: str, observed_version: str) -> Level | None:
    xc = _last(r"^\s*(\S+)\s+Method XC (?:Functional|Potential)", text)
    charge = _last(r"^\s*Charge\s+:\s+(-?\d+)", text)
    multiplicity = _last(r"^\s*Spin multiplicity:\s+(\d+)", text)
    if xc is None or charge is None or multiplicity is None:
        return None
    tol = _last(r"Convergence on energy requested:\s+(\S+)", text)
    return Level(
        program="nwchem", version=observed_version, method=xc.lower(), basis=basis(text),
        dispersion=_dispersion(text, xc.lower()), charge=int(charge), multiplicity=int(multiplicity),
        grid=_last(r"Grid used for XC integration:\s+(\S+)", text),
        scf_tol=None if tol is None else _number(tol),
    )


def _ccsd_t_level(text: str, observed_version: str) -> Level | None:
    """The state of the SCF reference: ``open shells`` of an RHF (or ROHF) one, alpha - beta
    electrons of a UHF one."""
    charge = _last(r"^\s*charge\s+=\s+(\S+)", text)
    unpaired = _last(r"^\s*open shells\s+=\s+(\d+)", text)
    if uhf := _UHF.findall(text):
        unpaired = str(int(uhf[-1][0]) - int(uhf[-1][1]))
    if charge is None or unpaired is None:
        return None
    return Level(program="nwchem", version=observed_version, method="ccsd(t)", basis=basis(text),
                 charge=round(_number(charge)), multiplicity=int(unpaired) + 1)


def observe_level(text: str) -> Level | None:
    """The Level the output reports, or None when a required field is missing."""
    observed_version = version(text)
    if observed_version is None:
        return None
    if _last(_CCSD_T, text) is None:
        return _dft_level(text, observed_version)
    return _ccsd_t_level(text, observed_version)


def total_energy(text: str) -> float | None:
    """Final electronic energy: the CCSD(T) total, else the last DFT total energy."""
    value = _last(_CCSD_T, text) or _last(_DFT_ENERGY, text)
    return None if value is None else _number(value)


def dft_energies(text: str) -> tuple[float, ...]:
    """The total DFT energy of every SCF, in order."""
    return tuple(_number(v) for v in re.findall(_DFT_ENERGY, text))


def step_energy(text: str) -> float | None:
    """The energy of the driver's last step line (``@ <step> <energy> ...``), at its last
    frame."""
    value = _last(r"^@\s+\d+\s+(\S+)", text)
    return None if value is None else _number(value)


def gradient(text: str) -> tuple[np.ndarray, np.ndarray] | None:
    """(coordinates in bohr, gradient in Eh/bohr), (N, 3) each, of the last ``DFT ENERGY
    GRADIENTS`` block."""
    rows: list[list[float]] = []
    for line in text[text.rfind("DFT ENERGY GRADIENTS"):].splitlines()[1:]:
        row = _GRADIENT_ROW.match(line)
        if row is None and rows:
            break
        if row is not None:
            rows.append([_number(v) for v in row[1].split()])
    return (np.array(rows)[:, :3], np.array(rows)[:, 3:]) if rows else None


def s2(text: str) -> float | None:
    value = _last(r"<S2> =\s+(\S+)", text)
    return None if value is None else _number(value)


def geometry_block(text: str, index: int = 0) -> tuple[tuple[str, ...], np.ndarray] | None:
    """Symbols and coordinates (Å) of a printed Geometry block; 0 is the input echo."""
    heads = list(_GEOMETRY.finditer(text))
    if not -len(heads) <= index < len(heads):
        return None
    scale: float | None = None
    symbols: list[str] = []
    coords: list[list[float]] = []
    for line in text[heads[index].end():].splitlines():
        if scale is None:
            if "Output coordinates in" in line:
                scale = 1.0 if "in angstroms" in line else BOHR_TO_ANGSTROM
            continue
        row = _ATOM_ROW.match(line)
        if row is None and symbols:
            break
        if row is not None:
            symbols.append(row[1].capitalize())
            coords.append([_number(v) for v in row.groups()[1:]])
    return (tuple(symbols), np.array(coords) * scale) if symbols and scale else None


def frame_shift(text: str, symbols: Sequence[str], coords: np.ndarray) -> float | None:
    """Largest atomic distance (Å) between the echoed input geometry and ``coords``; None
    when no geometry was echoed or its atoms differ."""
    echo = geometry_block(text, 0)
    if echo is None or echo[0] != tuple(symbols):
        return None
    diff = echo[1] - np.asarray(coords, dtype=float).reshape(-1, 3)
    return float(np.linalg.norm(diff, axis=1).max())


def count_frequency_blocks(text: str) -> int:
    return text.count(_PROJECTED)


def read_hess(path: Path, n_atoms: int) -> np.ndarray:
    """``.hess``: lower triangle row by row, Fortran D exponents, Eh/bohr² (symmetrized)."""
    values = [_number(v) for v in Path(path).read_text(encoding="ascii").split()]
    n = 3 * n_atoms
    if len(values) != n * (n + 1) // 2:
        raise ValueError(f"{path}: {len(values)} values do not fill a 3N = {n} lower triangle")
    h = np.zeros((n, n))
    h[np.tril_indices(n)] = values
    return h + np.tril(h, -1).T


def hess_to_npy(path: Path, n_atoms: int, dest: Path) -> Path:
    return to_canonical_npy(read_hess(path, n_atoms), dest)


def final_xyz(workdir: Path, symbols: Sequence[str]) -> Path | None:
    """The driver's last frame: the highest-numbered ``final-NNN.xyz`` (not the newest), or
    None when there is none or it is unreadable or changed the atom order."""
    frames = [(int(m[1]), path) for path in Path(workdir).glob("final-*.xyz")
              if (m := re.fullmatch(r"final-(\d+)\.xyz", path.name))]
    if not frames:
        return None
    path = max(frames)[1]
    try:
        return path if read_xyz(path).symbols == list(symbols) else None
    except ValueError:
        return None


def string_energies(text: str) -> tuple[float, ...]:
    """Final bead energies (``@zts Bead number``), in bead order."""
    pattern = r"^@zts Bead number\s+(\d+)\s+Potential Energy =\s+(\S+)"
    found = re.findall(pattern, text, re.MULTILINE)
    beads = {int(n): _number(e) for n, e in found}
    return tuple(beads[n] for n in sorted(beads))


def classify_failure(text: str, *, returncode: int | None, timed_out: bool) -> Failure | None:
    """The Failure of an abnormal job, or None when it terminated normally."""
    if timed_out or returncode == 124:
        return Failure(kind=FailureKind.TIMEOUT, reason="timeout")
    if returncode == 0 and NORMAL_END in text:
        return None
    for pattern, kind, reason in (
        (_AUTOZ, FailureKind.INPUT_INVALID, "autoz"),
        (_SCF, FailureKind.SCF_NOT_CONVERGED, "scf"),
        (_GEOMETRY_MAXITER, FailureKind.GEOMETRY_MAXITER, "maxiter"),
    ):
        if pattern.search(text):
            return Failure(kind=kind, reason=reason)
    if returncode not in (0, None):
        return Failure(kind=FailureKind.NONZERO_EXIT, reason=f"returncode {returncode}")
    return Failure(kind=FailureKind.INCOMPLETE_OUTPUT, reason="no_normal_termination")

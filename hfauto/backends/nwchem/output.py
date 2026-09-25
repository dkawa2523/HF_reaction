"""NWChem output parsing (design §5.1, §6.3): pure functions over the output text and the
files a job leaves in its permanent directory.

The Level is observed, never assumed: version, xc functional or wave-function method,
basis (``/cart`` when cartesian), DFT-D3 variant, COSMO dielectric, grid, SCF energy
tolerance, charge and multiplicity. Frequencies are never taken from the text for
decisions; the ``.hess`` file is converted to the canonical ``.npy`` instead.
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
_AUTOZ = re.compile(r"AUTOZ failed|regeneration of autoz failed")
_SCF = re.compile(r"Calculation failed to converge|SCF not converged")
_GEOMETRY_MAXITER = re.compile("Failed to converge in maximum number of steps")
_WFT_ENERGY = {
    "ccsd(t)": r"Total CCSD\(T\) energy:\s+(\S+)",
    "mp2": r"Total MP2 energy:?\s+(\S+)",
}


def _number(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


def _last(pattern: str, text: str) -> str | None:
    """First group of the last match (the final SCF / job wins)."""
    matches = re.findall(pattern, text, re.MULTILINE)
    return matches[-1] if matches else None


def version(text: str) -> str | None:
    return _last(r"Northwest Computational Chemistry Package \(NWChem\)\s+(\S+)", text)


def wft_method(text: str) -> str | None:
    """"ccsd(t)" or "mp2" when the job ran a correlated wave-function method."""
    if "Total CCSD(T) energy" in text:
        return "ccsd(t)"
    return "mp2" if "Total MP2 energy" in text else None


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
    eps = _last(r"^\s*dielec:\s+(\S+)", text)
    return Level(
        program="nwchem", version=observed_version, method=xc.lower(), basis=basis(text),
        dispersion=_dispersion(text, xc.lower()),
        solvation=None if eps is None else f"cosmo:{_number(eps):g}",
        charge=int(charge), multiplicity=int(multiplicity),
        grid=_last(r"Grid used for XC integration:\s+(\S+)", text),
        scf_tol=None if tol is None else _number(tol),
    )


def _wft_level(text: str, observed_version: str, method: str) -> Level | None:
    charge = _last(r"^\s*charge\s+=\s+(\S+)", text)
    open_shells = _last(r"^\s*open shells\s+=\s+(\d+)", text)
    if charge is None or open_shells is None:
        return None
    return Level(program="nwchem", version=observed_version, method=method, basis=basis(text),
                 charge=round(_number(charge)), multiplicity=int(open_shells) + 1)


def observe_level(text: str) -> Level | None:
    """The Level the output reports, or None when a required field is missing."""
    observed_version, wft = version(text), wft_method(text)
    if observed_version is None:
        return None
    if wft is None:
        return _dft_level(text, observed_version)
    return _wft_level(text, observed_version, wft)


def total_energy(text: str) -> float | None:
    """Final electronic energy: CCSD(T) / MP2 total, else the last DFT total energy."""
    value = _last(_WFT_ENERGY.get(wft_method(text) or "", r"Total DFT energy =\s+(\S+)"), text)
    return None if value is None else _number(value)


def trajectory_energies(text: str) -> tuple[float, ...]:
    """Driver energies per step (``@`` lines; the repeated last step counts once)."""
    steps = {int(n): _number(e) for n, e in re.findall(r"^@\s+(\d+)\s+(\S+)", text, re.MULTILINE)}
    return tuple(steps[n] for n in sorted(steps))


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


def last_frequency_block(text: str) -> tuple[float, ...]:
    """NWChem's own projected frequencies of the last block (cross-checks only)."""
    if _PROJECTED not in text:
        return ()
    tail = text.rpartition(_PROJECTED)[2]
    rows = re.findall(r"^\s*P\.Frequency\s+(.*)$", tail, re.MULTILINE)
    return tuple(_number(v) for row in rows for v in row.split())


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


def string_history(text: str) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """(gmax, xmax) per string iteration from ``string: gmax,grms,xrms,xmax=`` lines."""
    found = re.findall(r"^\s*string: gmax,grms,xrms,xmax=(.*)$", text, re.MULTILINE)
    rows = [row.split() for row in found if len(row.split()) == 4]
    return tuple(_number(r[0]) for r in rows), tuple(_number(r[3]) for r in rows)


def string_energies(text: str) -> tuple[float, ...]:
    """Final bead energies (``@zts Bead number``), in bead order."""
    pattern = r"^@zts Bead number\s+(\d+)\s+Potential Energy =\s+(\S+)"
    found = re.findall(pattern, text, re.MULTILINE)
    beads = {int(n): _number(e) for n, e in found}
    return tuple(beads[n] for n in sorted(beads))


def string_converged_by_program(text: str) -> bool:
    return "@zts The string calculation converged" in text


def classify_failure(text: str, *, returncode: int | None, timed_out: bool,
                     stopped: str | None = None) -> Failure | None:
    """The Failure of an abnormal job, or None when it terminated normally."""
    if timed_out or returncode == 124:
        return Failure(kind=FailureKind.TIMEOUT, reason="walltime")
    if stopped:
        return Failure(kind=FailureKind.STAGNATED, reason=stopped)
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

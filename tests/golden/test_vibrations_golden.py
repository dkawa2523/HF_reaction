"""G07 / G08: projected frequencies from NWChem .hess files (design §10.2, CH-26)."""

import re

import numpy as np
import pytest

from hfauto.chemistry.vibrations import projected_frequencies
from hfauto.chemistry.xyz import read_xyz

pytestmark = pytest.mark.golden


def _hess(text: str, n_atoms: int) -> np.ndarray:
    """NWChem .hess: lower triangle, row by row, Fortran D exponents, Eh/bohr²."""
    values = [float(v.replace("D", "E")) for v in text.split()]
    h = np.zeros((3 * n_atoms, 3 * n_atoms))
    h[np.tril_indices(3 * n_atoms)] = values
    return h + np.tril(h, -1).T


def test_g07_matches_nwchem_with_isotopic_masses(golden):
    xyz = read_xyz(golden.path("nwchem/G07/final.xyz"))
    hessian = _hess(golden.text("nwchem/G07/hfauto_job.hess"), 3)
    freqs, _, k = projected_frequencies(hessian, xyz.symbols, xyz.coords)
    out = golden.text("nwchem/G07/nwchem.out")
    lines = re.findall(r"^ P\.Frequency(.*)$", out, re.MULTILINE)
    nwchem = [float(v) for line in lines for v in line.split() if abs(float(v)) > 50]
    assert k == 6
    assert freqs == pytest.approx(nwchem, abs=1.0)
    assert freqs[0] == pytest.approx(-1131.57, abs=1.0)


def test_g08_non_stationary_seed_has_four_imaginary_modes(golden):
    deck = golden.text("nwchem/G08/nwchem.nw")
    block = re.search(r"^geometry[^\n]*\n(.*?)^end", deck, re.MULTILINE | re.DOTALL)
    rows = [line.split() for line in block.group(1).splitlines()]
    symbols, coords = [r[0] for r in rows], np.array([[float(v) for v in r[1:4]] for r in rows])
    hessian = _hess(golden.text("nwchem/G08/hfauto_job.hess"), len(symbols))
    freqs, _, _ = projected_frequencies(hessian, symbols, coords)
    assert np.count_nonzero(freqs < 0) == 4

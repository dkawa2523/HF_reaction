"""NWChem parser on real outputs (design §10.2: G01, G03-G05, G07, G10, G13, G21-G24)."""

import json
import re

import numpy as np
import pytest

from hfauto.backends.nwchem import output as nw
from hfauto.chemistry.profile import energies_settled
from hfauto.chemistry.vibrations import external_basis, projected_frequencies
from hfauto.chemistry.xyz import read_xyz
from hfauto.core.constants import BOHR_TO_ANGSTROM, HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import FailureKind as Kind
from hfauto.core.method import MethodSpec, level_mismatches

pytestmark = pytest.mark.golden
XFINE = MethodSpec(id="m", kind="dft", functional="pbe0", basis="def2-svpd",
                   dispersion="d3zero", grid="xfine", scf_energy_tol=1e-8)


def test_g01_two_vibrational_blocks(golden):
    assert nw.count_frequency_blocks(golden.text("nwchem/G01/nwchem.out")) == 2


@pytest.mark.parametrize(("stem", "energy", "steps"), [
    ("G03/hono_trans", -205.3455299945, 8), ("G03/hono_cis", -205.3452826520, 7),
    ("G21/water_reference", -76.2956762656, 5), ("G21/water_distorted", -76.2956762657, 8)])
def test_energy_trajectory_level_and_final_structure(golden, stem, energy, steps):
    text = golden.text(f"nwchem/{stem}.out")
    assert level_mismatches(XFINE, nw.observe_level(text), version_pin="7.2.3") == []
    assert nw.total_energy(text) == pytest.approx(energy, abs=1e-9)
    trajectory = nw.trajectory_energies(text)
    assert len(trajectory) == steps and trajectory[-1] == pytest.approx(energy, abs=1e-8)
    final = read_xyz(golden.path(f"nwchem/{stem}_final.xyz"))
    symbols, coords = nw.geometry_block(text, -1)
    assert symbols == tuple(final.symbols) and np.allclose(coords, final.coords, atol=1e-8)


@pytest.mark.parametrize("stem", ["G04/hnc", "G05/hcn"])
def test_linear_molecules_have_five_external_modes(golden, stem):
    text = golden.text(f"nwchem/{stem}.out")
    symbols, coords = nw.geometry_block(text, -1)
    assert external_basis(symbols, coords).shape[1] == 5
    printed = [float(v) for row in re.findall(r"^\s*P\.Frequency(.*)$", text, re.MULTILINE)
               for v in row.split()]  # NWChem's own projection, one block
    assert len([f for f in printed if abs(f) > 1.0]) == 3 * 3 - 5


def test_g07_hess_to_canonical_npy(golden, tmp_path):
    xyz = read_xyz(golden.path("nwchem/G07/final.xyz"))
    hess = golden.path("nwchem/G07/hfauto_job.hess")
    npy = nw.hess_to_npy(hess, 3, tmp_path / "hessian.npy")
    freqs, _, n_external = projected_frequencies(np.load(npy), xyz.symbols, xyz.coords)
    assert n_external == 6 and freqs[0] == pytest.approx(-1131.6, abs=1.0)
    with pytest.raises(ValueError, match="3N"):  # ported K case: 3N mismatch
        nw.read_hess(hess, 4)


def test_g10_grid_mismatch_and_frame_change(golden):
    text = golden.text("nwchem/G10/irc_000.001.qce_nwchem_stdout")
    level = nw.observe_level(text)
    assert (level.grid, level.dispersion, level.scf_tol, level.basis) == (
        "medium", "d3zero", 1e-6, "def2-svpd")
    assert any(m.startswith("grid") for m in level_mismatches(XFINE, level, version_pin="7.2.3"))
    rows = [r.split() for r in text.split("geometry units bohr\n")[1].split("\nend")[0].split("\n")
            if r.strip()]
    deck = np.array([[float(v) for v in r[1:]] for r in rows]) * BOHR_TO_ANGSTROM
    assert nw.frame_shift(text, [r[0] for r in rows], deck) > 1e-4  # rigidly moved: frame lost
    assert nw.frame_shift(text, *nw.geometry_block(text, 0)) == 0.0


def test_g13_bead_energies_settle_although_gmax_rose(golden):
    text = golden.text("nwchem/G13/nwchem_string.out")  # gmax 9.1e-4 -> 9.2e-3
    history = nw.string_path_energies(text)
    assert [len(energies) for energies in history] == [11] * 3
    assert history[-1] == pytest.approx(nw.string_energies(text), abs=1e-9)
    assert energies_settled(history, 0.1 / HARTREE_TO_KCAL_MOL)
    assert nw.program_converged(text)  # recorded only


def test_g22_g23_g24_and_other_failures(golden):
    timeout = json.loads(golden.text("nwchem/G22/command_result.json"))
    killed = json.loads(golden.text("nwchem/G24/command_result.json"))
    failures = [
        nw.classify_failure(golden.text("nwchem/G22/nwchem.out"),
                            returncode=timeout["returncode"], timed_out=timeout["timed_out"]),
        nw.classify_failure(golden.text("nwchem/G23/nwchem.out"), returncode=1, timed_out=False),
        nw.classify_failure("", returncode=killed["returncode"], timed_out=killed["timed_out"]),
        nw.classify_failure(" AUTOZ failed to generate", returncode=1, timed_out=False),
        nw.classify_failure("Calculation failed to converge", returncode=1, timed_out=False),
        nw.classify_failure("", returncode=0, timed_out=False),
    ]
    assert [f.kind for f in failures] == [
        Kind.TIMEOUT, Kind.GEOMETRY_MAXITER, Kind.NONZERO_EXIT, Kind.INPUT_INVALID,
        Kind.SCF_NOT_CONVERGED, Kind.INCOMPLETE_OUTPUT]
    ok = " AUTOZ failed to generate good internal coordinates.\n" + nw.NORMAL_END
    assert nw.classify_failure(ok, returncode=0, timed_out=False) is None  # NWChem fell back


def test_final_xyz_numbering_atom_order_and_missing_d3(golden, tmp_path):
    frame = "3\n geometry\n{}\n{}\n{}\n"
    atoms = ("C 0 0 0", "N 0 0 1.16", "H 0 0 -1.07")
    (tmp_path / "final-1000.xyz").write_text(frame.format(*atoms))
    for n in (9, 999):  # written later, but numbered lower
        (tmp_path / f"final-{n:03d}.xyz").write_text(frame.format(*atoms))
    assert nw.final_xyz(tmp_path, ["C", "N", "H"]).name == "final-1000.xyz"
    assert nw.final_xyz(tmp_path, ["N", "C", "H"]) is None  # ported K case: atom order
    text = golden.text("nwchem/G03/hono_trans.out").replace("DFT-D3 Model", "")
    mismatches = level_mismatches(XFINE, nw.observe_level(text), version_pin="7.2.3")
    assert [m.split(":")[0] for m in mismatches] == ["dispersion"]  # ported K case: no D3

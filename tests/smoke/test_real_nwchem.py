"""Real NWChem 7.2.3 (WSL, HFAUTO_REAL=1, ranks <= 2; design §10.3).

The shaped-Hessian case uses planar hydroxylamine (N, O and all three H in one plane): at
PBE0-D3BJ/def2-SVP its Hessian has exactly two negative eigenvalues (NH2 inversion near
-1019 cm-1, OH torsion near -666 cm-1). Passing the torsion as ``mode`` makes it the only
negative mode of the shaped initial Hessian, which the saddle follows in autoz with moddir 1.
"""

import re
import shutil
from pathlib import Path

import numpy as np
import pytest

from hfauto.backends.protocols import Capability
from hfauto.chemistry.gates import is_first_order_saddle, is_minimum
from hfauto.chemistry.profile import classify
from hfauto.chemistry.xyz import XYZ, Molecule, read_xyz, read_xyz_trajectory
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import Evidence, FileRef, PathProfile
from hfauto.core.hashing import sha256_file
from hfauto.core.method import MethodSpec

pytestmark = pytest.mark.real
DATA = Path(__file__).parent / "data"
PBE0 = MethodSpec(id="pbe0-d3bj_def2-svp", kind="dft", functional="pbe0", basis="def2-svp",
                  dispersion="d3bj", grid="fine", scf_energy_tol=1e-7)


def _mol(*rows: tuple) -> Molecule:
    return Molecule(XYZ([r[0] for r in rows], np.array([r[1:] for r in rows], dtype=float)), 0, 1)


WATER = _mol(("O", 0, 0, 0.1273), ("H", 0, 0.7872, -0.4692), ("H", 0, -0.7472, -0.4792))
PLANAR_NH3 = _mol(("N", 0, 0, 0), ("H", 1.01, 0, 0), ("H", -0.505, 0.874686, 0),
                  ("H", -0.505, -0.874686, 0))
PLANAR_NH2OH = _mol(("N", 0, 0, 0), ("O", 1.42, 0, 0), ("H", -0.5, 0.866, 0),
                    ("H", -0.5, -0.866, 0), ("H", 1.74, 0.91, 0))
HCN = _mol(("H", 0, 0, -1.066), ("C", 0, 0, 0), ("N", 0, 0, 1.156))
HI = _mol(("H", 0, 0, 0), ("I", 0, 0, 1.61))
OH = Molecule(XYZ(["O", "H"], np.array([[0, 0, 0.00777717], [0, 0, 0.98222283]])), 0, 2)


def _final(tmp_path, ev: Evidence) -> Molecule:
    return Molecule(read_xyz(tmp_path / "run" / ev.final.file.path), 0, 1)


def test_water_opt_then_a_separate_freq_job(real_engine, tmp_path):
    qm = real_engine(Capability.QM, "nwchem")
    opt = qm.optimize(WATER, PBE0)
    assert isinstance(opt, Evidence), opt
    freq = qm.frequencies(_final(tmp_path, opt), PBE0)  # one block, Level as requested
    assert isinstance(freq, Evidence), freq
    assert freq.level == opt.level and is_minimum(freq, opt=opt)


def test_planar_ammonia_saddle_from_the_freq_hessian(real_engine, tmp_path):
    qm = real_engine(Capability.QM, "nwchem")
    saddle = real_engine(Capability.SADDLE, "nwchem_saddle")
    hessian = qm.frequencies(PLANAR_NH3, PBE0)
    ts = saddle.refine(PLANAR_NH3, PBE0, hessian=hessian, mode=hessian.imaginary_modes[0])
    assert isinstance(ts, Evidence) and ts.task == "saddle", ts
    freq = qm.frequencies(_final(tmp_path, ts), PBE0)
    assert is_first_order_saddle(freq, saddle=ts)
    assert sum(f < -50 for f in freq.frequencies_cm1) == 1


def test_shaped_hessian_follows_the_second_negative_mode(real_engine):
    qm = real_engine(Capability.QM, "nwchem")
    saddle = real_engine(Capability.SADDLE, "nwchem_saddle")
    hessian = qm.frequencies(PLANAR_NH2OH, PBE0)
    assert sum(f < -50 for f in hessian.frequencies_cm1) == 2
    ts = saddle.refine(PLANAR_NH2OH, PBE0, hessian=hessian, mode=hessian.imaginary_modes[1])
    assert isinstance(ts, Evidence) and ts.task == "saddle", ts


def test_iodine_def2_ecp_and_a_chloride_atom(real_engine, tmp_path):
    """HI: the def2-ECP of I leaves 26 electrons (alpha 13); Cl-: opt and freq of one atom."""
    qm = real_engine(Capability.QM, "nwchem")
    sp = qm.energy(HI, PBE0)
    assert isinstance(sp, Evidence), sp
    assert re.search(r"Alpha electrons :\s+13\n", (tmp_path / "run" / sp.output.path).read_text())
    chloride = Molecule(XYZ(["Cl"], np.zeros((1, 3))), -1, 1)
    opt = qm.optimize(chloride, PBE0)
    assert isinstance(opt, Evidence), opt
    freq = qm.frequencies(Molecule(_final(tmp_path, opt).xyz, -1, 1), PBE0)
    assert isinstance(freq, Evidence), freq
    assert (freq.frequencies_cm1, freq.n_external) == ((), 3) and is_minimum(freq, opt=opt)


def test_closed_shell_ccsd_t_and_the_wb97x_d3_level(real_engine):
    qm = real_engine(Capability.QM, "nwchem")
    ccsd_t = MethodSpec(id="ccsd-t", kind="wft", wft_method="ccsd(t)", basis="def2-svp")
    ev = qm.energy(HCN, ccsd_t)  # the RHF ccsd module
    assert isinstance(ev, Evidence), ev
    assert (ev.level.method, ev.level.multiplicity, ev.s2) == ("ccsd(t)", 1, None)
    wb97 = MethodSpec(id="wb97x-d3", kind="dft", functional="wb97x-d3", basis="def2-svp",
                      grid="fine")
    ev = qm.energy(WATER, wb97)
    assert isinstance(ev, Evidence), ev
    assert (ev.level.method, ev.level.dispersion, ev.level.grid) == ("wb97x-d3", None, "fine")


def test_open_shell_ccsd_t_is_uhf_through_the_tce(real_engine):
    """OH.: UHF-CCSD(T)/def2-TZVPD through the TCE (P0d)."""
    qm = real_engine(Capability.QM, "nwchem")
    ccsd_t = MethodSpec(id="ccsd-t_def2-tzvpd", kind="wft", wft_method="ccsd(t)",
                        basis="def2-tzvpd")
    ev = qm.energy(OH, ccsd_t)
    assert isinstance(ev, Evidence), ev
    assert (ev.level.method, ev.level.multiplicity, ev.s2) == ("ccsd(t)", 2, None)
    assert ev.energy_hartree == pytest.approx(-75.640465704893, abs=1e-6)


def test_one_string_chunk_is_classified_as_it_stands_hcn_sto3g(real_engine, tmp_path):
    """HCN -> HNC: one ZTS chunk from a low-level path (9 beads, 20 iterations, about 80 s):
    gmax stays near 0.05 Eh/bohr, yet its beads show a single maximum (U5-P4)."""
    run = tmp_path / "run"
    run.mkdir(parents=True, exist_ok=True)
    path = Path(shutil.copy(DATA / "hcn_zts_initial.xyz", run / "hcn_zts_initial.xyz"))
    frames = read_xyz_trajectory(path)
    start, end = (Molecule(frames[i], 0, 1) for i in (0, -1))
    sto3g = MethodSpec(id="pbe0_sto-3g", kind="dft", functional="pbe0", basis="sto-3g",
                       scf_energy_tol=1e-7)
    initial = FileRef(path=path.name, sha256=sha256_file(path))
    string = real_engine(Capability.PATH, "nwchem_string")
    profile = string.find_path(start, end, sto3g, images=9, initial_path=initial)
    assert isinstance(profile, PathProfile), profile
    assert classify(profile.energies_hartree, 1 / HARTREE_TO_KCAL_MOL) == "single"

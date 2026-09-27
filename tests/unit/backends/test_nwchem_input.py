"""NWChem renderer (design §6.3): job split, Hessian input, saddle / string keywords, states."""

import numpy as np
import pytest

from hfauto.backends.nwchem import input as nw
from hfauto.backends.nwchem.output import read_hess
from hfauto.chemistry.xyz import XYZ, Molecule
from hfauto.core.method import MethodSpec

WATER = Molecule(XYZ(["O", "H", "H"], np.array([[0, 0, 0.1173], [0, 0.7572, -0.4692],
                                                [0, -0.7572, -0.4692]])), 0, 1)
PBE0 = MethodSpec(id="pbe0", kind="dft", functional="PBE0", basis="def2-SVPD",
                  dispersion="d3bj", grid="fine", scf_energy_tol=1e-7)
MP2 = MethodSpec(id="mp2", kind="wft", wft_method="mp2", basis="def2-svp")


def test_opt_and_saddle_never_compute_hessians():
    opt, saddle = nw.render_optimize(WATER, PBE0, init_hessian=True), nw.render_saddle(WATER, PBE0)
    for deck in (opt, saddle):
        assert "frequencies" not in deck and "hessian" not in deck
        assert "inhess 2" in deck and "xyz final" in deck and deck.count("task ") == 1
    assert "maxiter 100" in opt and "task dft optimize" in opt
    assert "trust 0.3" in opt  # a QRC side from its TS Hessian (S21)
    plain = nw.render_optimize(WATER, PBE0)  # mode-follow starts must not overshoot the TS
    assert "inhess" not in plain and "trust 0.1" in plain
    freq = nw.render_frequencies(WATER, PBE0)
    assert "task dft frequencies" in freq and freq.count("task ") == 1 and "driver" not in freq


def test_saddle_keywords_and_mode_following():
    free = nw.render_saddle(WATER, PBE0, moddir=0)
    for key in ("trust 0.1", "sadstp 0.1", "maxiter 50", "inhess 2", "task dft saddle"):
        assert key in free
    assert "moddir" not in free and "noautoz" not in free
    lowest = nw.render_saddle(WATER, PBE0, moddir=1)
    assert "  moddir 1" in lowest and "noautoz" not in lowest
    second = nw.render_saddle(WATER, PBE0, moddir=2, cartesian=True)
    assert "  moddir 2" in second and "noautosym noautoz" in second


def test_string_keywords():
    deck = nw.render_string(WATER, WATER, PBE0, nbeads=9, initial_path=True)
    for key in ("nbeads 9", "maxiter 20", "stepsize 0.05", "interpol 3", "tol 1e-5", "impose",
                "freeze1 .true.", "freezeN .true.", f"xyz_path {nw.INITIAL_PATH}",
                "geometry endgeom units angstrom", "task dft string"):
        assert key in deck
    assert "xyz_path" not in nw.render_string(WATER, WATER, PBE0, nbeads=9)


def test_charge_multiplicity_dispersion_solvation_and_frame():
    closed = nw.render_energy(WATER, PBE0)
    for key in ("charge 0", "mult 1", "disp vdw 4", "basis spherical", "* library def2-SVPD",
                "xc pbe0", "grid fine", "convergence energy 1.0e-07", "permanent_dir ."):
        assert key in closed
    assert "odft" not in closed and "noautoz" not in closed and "cosmo" not in closed
    default_tol = nw.render_frequencies(WATER, PBE0.model_copy(update={"scf_energy_tol": None}))
    assert "convergence energy 1.0e-07" in default_tol  # opt and freq share numerics
    other = PBE0.model_copy(update={"dispersion": "d3zero", "solvation": "cosmo:78.4"})
    setup = nw.Setup(cartesian=True, scratch_dir="/scr/job", restart_vectors=True)
    cation = nw.render_energy(Molecule(WATER.xyz, 1, 2), other, setup)
    for key in ("charge 1", "mult 2", "odft", "disp vdw 3", "dielec 78.4", "noautoz",
                "scratch_dir /scr/job", "vectors input job.movecs"):
        assert key in cation
    for deck in (closed, cation, nw.render_wft(WATER, MP2), nw.render_saddle(WATER, PBE0)):
        assert "units angstrom nocenter noautosym" in deck
    with pytest.raises(ValueError):
        nw.render_energy(WATER, PBE0.model_copy(update={"dispersion": "d4"}))


def test_def2_writes_one_ecp_line_per_element_beyond_kr():
    ch3i = Molecule(XYZ(["C", "I", "H", "H", "H"], np.zeros((5, 3))), 0, 1)
    ccsd_t = MP2.model_copy(update={"wft_method": "ccsd(t)", "basis": "def2-tzvpd"})
    for deck in (nw.render_energy(ch3i, PBE0), nw.render_wft(ch3i, ccsd_t),
                 nw.render_string(ch3i, ch3i, PBE0, nbeads=5)):
        assert "end\necp\n  I library def2-ecp\nend" in deck and "* library def2-ecp" not in deck
    assert "ecp" not in nw.render_energy(WATER, PBE0)
    assert nw.ecp(["Te", "I", "Xe", "I", "H"], "def2-svp") == [
        "ecp", "  Te library def2-ecp", "  I library def2-ecp", "  Xe library def2-ecp", "end"]
    with pytest.raises(ValueError, match="no_ecp_for_basis:cc-pVDZ:I"):  # no all-electron I
        nw.render_wft(ch3i, ccsd_t.model_copy(update={"basis": "cc-pVDZ"}))


def test_wft_single_points_and_hess_round_trip(tmp_path):
    mp2 = nw.render_wft(WATER, MP2)
    assert "task mp2 energy" in mp2 and "freeze atomic" in mp2 and "dft" not in mp2
    assert "scf" not in mp2 and "maxiter" not in mp2
    ccsd_t = MP2.model_copy(update={"wft_method": "ccsd(t)"})
    ccsd = nw.render_wft(WATER, ccsd_t, nw.Setup(restart_vectors=True))
    assert "task ccsd(t) energy" in ccsd and "freeze atomic\n  maxiter 50\nend" in ccsd
    assert "scf\n  vectors input job.movecs\nend" in ccsd
    assert "nopen" not in ccsd and "uhf" not in ccsd
    h = np.arange(81.0).reshape(9, 9) * 1e-3
    (tmp_path / "job.hess").write_text(nw.hess_text(h + h.T))
    assert np.allclose(read_hess(tmp_path / "job.hess", 3), h + h.T, rtol=1e-10)

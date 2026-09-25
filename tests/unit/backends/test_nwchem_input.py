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
    assert "maxiter 100" in opt and "task dft optimize" in opt and "tight" not in opt
    assert "trust 0.1" in opt  # QRC / mode-follow starts must not overshoot above the TS
    assert "  tight" in nw.render_optimize(WATER, PBE0, tight=True)
    assert "inhess" not in nw.render_optimize(WATER, PBE0)
    freq = nw.render_frequencies(WATER, PBE0)
    assert "task dft frequencies" in freq and freq.count("task ") == 1 and "driver" not in freq


def test_saddle_keywords_and_mode_following():
    lowest = nw.render_saddle(WATER, PBE0, mode_index=0)
    for key in ("trust 0.1", "sadstp 0.1", "maxiter 50", "inhess 2", "task dft saddle"):
        assert key in lowest
    assert "moddir" not in lowest and "noautoz" not in lowest
    second = nw.render_saddle(WATER, PBE0, mode_index=1)
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


def test_wft_single_points_and_hess_round_trip(tmp_path):
    mp2 = nw.render_wft(WATER, MP2)
    assert "task mp2 energy" in mp2 and "freeze atomic" in mp2 and "dft" not in mp2
    ccsd_t = MP2.model_copy(update={"wft_method": "ccsd(t)"})
    ccsd = nw.render_wft(Molecule(WATER.xyz, 1, 2), ccsd_t)
    assert "task ccsd(t) energy" in ccsd and "nopen 1" in ccsd and "freeze atomic" in ccsd
    h = np.arange(81.0).reshape(9, 9) * 1e-3
    (tmp_path / "job.hess").write_text(nw.hess_text(h + h.T))
    assert np.allclose(read_hess(tmp_path / "job.hess", 3), h + h.T, rtol=1e-10)

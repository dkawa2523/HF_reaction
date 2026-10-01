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
CCSD_T = MethodSpec(id="ccsd-t", kind="wft", wft_method="ccsd(t)", basis="def2-svp")


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


def test_a_fixed_bond_is_a_frozen_zcoord_bond_or_a_cartesian_spring():
    """G2-P4: 1-based atoms inside the geometry; in Cartesian coordinates (the autoz failure's
    continuation) a restraint k (r - r0)² with r0 in bohr (NWChem 7.2.3, M8 probe)."""
    deck = nw.render_optimize(WATER, PBE0, fixed_bond=(0, 2, 1.2345))
    assert "-0.4692000000\n  zcoord\n    bond 1 3 1.2345 rc constant\n  end\nend\ncharge 0" in deck
    assert "constraints" not in deck and "noautoz" not in deck
    spring = nw.render_optimize(WATER, PBE0, nw.Setup(cartesian=True), fixed_bond=(0, 2, 1.2345))
    assert "\nend\n\nconstraints\n  spring bond 1 3 20.0 2.332867\nend\n\ndft\n" in spring
    assert "noautosym noautoz" in spring and "zcoord" not in spring
    for plain in (nw.render_optimize(WATER, PBE0), nw.render_saddle(WATER, PBE0)):
        assert "zcoord" not in plain and "constraints" not in plain


def test_saddle_keywords_always_follow_mode_1():
    deck = nw.render_saddle(WATER, PBE0)  # the shaped Hessian's only negative mode
    for key in ("trust 0.1", "sadstp 0.1", "maxiter 50", "inhess 2\n  moddir 1",
                "task dft saddle"):
        assert key in deck
    assert "noautoz" not in deck
    resumed = nw.render_saddle(WATER, PBE0, nw.Setup(cartesian=True), init_hessian=False)
    assert "  moddir 1" in resumed and "noautosym noautoz" in resumed and "inhess" not in resumed


def test_string_keywords():
    deck = nw.render_string(WATER, WATER, PBE0, nbeads=9, initial_path=True)
    for key in ("nbeads 9", "maxiter 20", "stepsize 0.05", "interpol 3", "tol 1e-5", "impose",
                "freeze1 .true.", "freezeN .true.", f"xyz_path {nw.INITIAL_PATH}",
                "geometry endgeom units angstrom", "task dft string"):
        assert key in deck
    assert "xyz_path" not in nw.render_string(WATER, WATER, PBE0, nbeads=9)


def test_charge_multiplicity_dispersion_and_frame():
    closed = nw.render_energy(WATER, PBE0)
    for key in ("charge 0", "mult 1", "disp vdw 4", "basis spherical", "* library def2-SVPD",
                "xc pbe0", "grid fine", "convergence energy 1.0e-07", "iterations 100",
                "permanent_dir ."):
        assert key in closed
    assert "odft" not in closed and "noautoz" not in closed and "cosmo" not in closed
    assert "cgmin" not in closed
    default_tol = nw.render_frequencies(WATER, PBE0.model_copy(update={"scf_energy_tol": None}))
    assert "convergence energy 1.0e-07" in default_tol  # opt and freq share numerics
    other = PBE0.model_copy(update={"dispersion": "d3zero"})
    setup = nw.Setup(cartesian=True, scratch_dir="/scr/job", restart_vectors=True,
                     scf_rescue=True)
    cation = nw.render_energy(Molecule(WATER.xyz, 1, 2), other, setup)
    for key in ("charge 1", "mult 2", "odft", "disp vdw 3", "noautoz",
                "scratch_dir /scr/job", "vectors input job.movecs\n  cgmin\nend"):
        assert key in cation
    for deck in (nw.render_optimize(WATER, PBE0), nw.render_saddle(WATER, PBE0),
                 nw.render_frequencies(WATER, PBE0),
                 nw.render_string(WATER, WATER, PBE0, nbeads=5)):
        assert "\n  iterations 100\n" in deck
    for deck in (closed, cation, nw.render_wft(WATER, CCSD_T), nw.render_saddle(WATER, PBE0)):
        assert "units angstrom nocenter noautosym" in deck


def test_def2_writes_one_ecp_line_per_element_beyond_kr():
    ch3i = Molecule(XYZ(["C", "I", "H", "H", "H"], np.zeros((5, 3))), 0, 1)
    ccsd_t = CCSD_T.model_copy(update={"basis": "def2-tzvpd"})
    for deck in (nw.render_energy(ch3i, PBE0), nw.render_wft(ch3i, ccsd_t),
                 nw.render_string(ch3i, ch3i, PBE0, nbeads=5)):
        assert "end\necp\n  I library def2-ecp\nend" in deck and "* library def2-ecp" not in deck
    assert "ecp" not in nw.render_energy(WATER, PBE0)
    assert nw.ecp(["Te", "I", "Xe", "I", "H"], "def2-svp") == [
        "ecp", "  Te library def2-ecp", "  I library def2-ecp", "  Xe library def2-ecp", "end"]
    with pytest.raises(ValueError, match="no_ecp_for_basis:cc-pVDZ:I"):  # no all-electron I
        nw.render_wft(ch3i, ccsd_t.model_copy(update={"basis": "cc-pVDZ"}))


def test_ccsd_t_single_point_and_hess_round_trip(tmp_path):
    assert "scf\n  maxiter 100\nend" in nw.render_wft(WATER, CCSD_T)
    ccsd = nw.render_wft(WATER, CCSD_T, nw.Setup(restart_vectors=True))
    assert "task ccsd(t) energy" in ccsd and "dft" not in ccsd and "ccsd\n  freeze atomic\n  maxiter 50\nend" in ccsd
    assert "scf\n  maxiter 100\n  vectors input job.movecs\nend" in ccsd
    assert "nopen" not in ccsd and "rohf" not in ccsd and "tce" not in ccsd
    h = np.arange(81.0).reshape(9, 9) * 1e-3
    (tmp_path / "job.hess").write_text(nw.hess_text(h + h.T))
    assert np.allclose(read_hess(tmp_path / "job.hess", 3), h + h.T, rtol=1e-10)


@pytest.mark.parametrize(("multiplicity", "nopen"), [(2, 1), (3, 2)])
def test_open_shell_ccsd_t_is_rohf_through_the_tce(multiplicity, nopen):
    radical = Molecule(WATER.xyz, 1, multiplicity)
    deck = nw.render_wft(radical, CCSD_T)
    assert f"scf\n  rohf\n  nopen {nopen}\n  maxiter 100\nend" in deck
    assert "tce\n  2eorb\n  2emet 13\n  ccsd(t)\n  freeze atomic\nend" in deck
    assert deck.rstrip().endswith("task tce energy") and "\nccsd\n" not in deck

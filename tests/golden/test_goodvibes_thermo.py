"""species_thermo (GoodVibes 4.3.0 in this process; Linux production extra) against the G06
CSVs of the CLI (`-q --qs grimme`: quasi-harmonic H and S), GoodVibes' own output parser and
textbook values. The point group is chemistry.symmetry's on the excerpt's own Hessian where it
has one; else it is given."""

import csv
import math
import re

import numpy as np
import pytest

from hfauto.backends.nwchem.output import geometry_block, read_hess
from hfauto.chemistry import thermo as th
from hfauto.chemistry.symmetry import Symmetry, analyze
from hfauto.chemistry.xyz import XYZ
from hfauto.core.constants import HARTREE_TO_KCAL_MOL, R_KCAL_MOL_K
from hfauto.core.method import ThermoSettings

pytestmark = pytest.mark.golden
pytest.importorskip("goodvibes")
T, RT = 298.15, R_KCAL_MOL_K * 298.15 / HARTREE_TO_KCAL_MOL
MAIN, G06 = ThermoSettings(), ThermoSettings(vib_scale=0.985)  # G06 ran at vib_scale 0.985
# The CLI read NWChem's rotational constants (7 digits), species_thermo takes them from the
# geometry: -T S_rot differs by 1.3e-8 (HCN, HNC) and 2.0e-8 Eh (TS).
ROTATION_TOL = 3e-8


def _freq(golden, name: str) -> tuple[XYZ, list[float]]:
    """The last geometry and every P.Frequency except the external modes of an NWChem freq."""
    text = golden.text(name)
    lines = re.findall(r"^ P\.Frequency(.*)$", text, re.MULTILINE)
    symbols, coords = geometry_block(text, -1)
    return XYZ(list(symbols), coords), [
        nu for line in lines for nu in map(float, line.split()) if abs(nu) > 1.0]


def _given(xyz: XYZ, group: str = "C1", sigma: int = 1, linear: bool = False) -> Symmetry:
    return Symmetry(group, sigma, linear, 1, xyz.coords)


def _thermo(xyz, modes, sym=None, *, saddle=False, settings=MAIN) -> th.Thermal:
    return th.species_thermo(xyz.symbols, sym or _given(xyz), th.thermo_frequencies(
        modes, saddle=saddle), settings=settings, T=T)


def _csv(golden, name: str) -> dict[str, str]:
    return next(csv.DictReader(golden.text(f"goodvibes/{name}.csv").splitlines()))


@pytest.mark.parametrize(("name", "out", "saddle"), [
    ("hcn", "G05/hcn.out", False), ("hnc", "G04/hnc.out", False), ("ts", "G07/nwchem.out", True)])
def test_g06_linear_minima_and_ts_match_the_cli(golden, name, out, saddle):
    """The TS (G07) is Cs on its own Hessian; HCN and HNC come without one and are C∞v, as
    NWChem reports (symmetry # 1)."""
    row = {k: float(v) for k, v in _csv(golden, f"G06/{name}").items()
           if k in ("scf_energy", "zpe", "qh_enthalpy", "entropy", "qh_gibbs_free_energy")}
    xyz, modes = _freq(golden, f"nwchem/{out}")
    sym = _given(xyz, "Cinfv", 1, linear=True)
    if saddle:
        hessian = read_hess(golden.path("nwchem/G07/hfauto_job.hess"), 3)
        sym = analyze(xyz.symbols, xyz.coords, hessian)
        assert (sym.point_group, sym.sigma, sym.linear, sym.m) == ("Cs", 1, False, 1)
    grimme = _thermo(xyz, modes, sym, saddle=saddle, settings=G06)
    truhlar = _thermo(xyz, modes, sym, saddle=saddle,
                      settings=G06.model_copy(update={"qs": "truhlar"}))
    H = row["qh_enthalpy"] - row["scf_energy"]
    assert (grimme.zpe, grimme.H) == pytest.approx((row["zpe"], H), abs=1e-10)
    assert grimme.G == pytest.approx(row["qh_gibbs_free_energy"] - row["scf_energy"],
                                     abs=ROTATION_TOL)
    assert truhlar.G == pytest.approx(H - T * row["entropy"], abs=ROTATION_TOL)  # modes > 100


def test_g02_came_from_both_frequency_blocks_of_the_ts(golden):
    """G02's 24.11 kcal/mol had the ZPE of both freq blocks of G01; species_thermo takes the
    modes it is given."""
    xyz, modes = _freq(golden, "nwchem/G01/nwchem.out")  # two blocks of six modes
    half = len(modes) // 2
    zpe = [_thermo(xyz, m, saddle=True, settings=G06).zpe for m in (modes[:half], modes[half:])]
    assert float(_csv(golden, "G02/ts")["zpe"]) == pytest.approx(sum(zpe), abs=1e-10)


def test_doublet_and_linear_oh(golden):
    """G25: OH radical, one stiff mode, C∞v (an open shell: the noise criterion holds for a
    diatomic). GoodVibes sees a singlet; the electronic term is thermo.electronic's."""
    xyz, modes = _freq(golden, "nwchem/G25/oh_freq.out")
    hessian = read_hess(golden.path("nwchem/G25/oh_freq.hess"), 2)
    sym = analyze(xyz.symbols, xyz.coords, hessian, open_shell=True)
    assert (sym.point_group, sym.sigma, sym.linear, sym.m) == ("Cinfv", 1, True, 1)
    oh = _thermo(xyz, modes, sym)
    # trans + linear rot + pV; QH leaves (100/3700)^4 of the stretch's ZPE out: 5e-9 Eh
    assert oh.H - oh.zpe == pytest.approx(3.5 * RT, abs=1e-8)
    assert th.electronic(((2, 0.0),), T) == pytest.approx(-RT * math.log(2), rel=1e-6)


def test_argon_standard_entropy():
    """An atom: translation only. JANAF S(298.15 K, 1 bar) = 154.846 J/mol/K; GoodVibes refers
    to 1 atm, and S(1 bar) = S(1 atm) + R ln 1.01325."""
    atom = XYZ(["Ar"], np.zeros((1, 3)))
    ar = _thermo(atom, (), analyze(atom.symbols, atom.coords, np.zeros((3, 3))))
    s_1atm = (ar.H - ar.G) / T * HARTREE_TO_KCAL_MOL * 4184.0  # J/mol/K
    assert s_1atm + 8.314462618 * math.log(1.01325) == pytest.approx(154.846, abs=0.01)
    assert ar.zpe == 0.0 and ar.H == pytest.approx(2.5 * RT, rel=1e-6)


def test_cutoff_and_truhlar_variants_match_goodvibes_parsing_the_output(golden):
    """G23 (TMA-(HF)2 saddle, real modes from 84.6 cm-1): each variant minus the main settings,
    in process and from GoodVibes reading the output itself (S_rot and S_trans cancel)."""
    from goodvibes.api import compute_thermo

    path = str(golden.path("nwchem/G23/nwchem.out"))
    xyz, modes = _freq(golden, "nwchem/G23/nwchem.out")
    real = [nu for nu in modes if nu > 0.0]  # GoodVibes leaves out imaginary modes

    def cli(s: ThermoSettings) -> float:
        return compute_thermo(path, QS=s.qs, QH=True, s_freq_cutoff=s.cutoff_cm1,
                              h_freq_cutoff=s.cutoff_cm1, freq_scale_factor=1.0,
                              zpe_scale_factor=1.0, symm=True).qh_gibbs_free_energy

    for qs, cutoff in (("grimme", 50.0), ("grimme", 150.0), ("truhlar", 50.0), ("truhlar", 150.0)):
        variant = MAIN.model_copy(update={"qs": qs, "cutoff_cm1": cutoff})
        shift = _thermo(xyz, real, settings=variant).G - _thermo(xyz, real).G
        assert abs(shift) > 1e-5 and shift == pytest.approx(cli(variant) - cli(MAIN), abs=1e-10)

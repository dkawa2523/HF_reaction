"""Real NWChem freq -> species_thermo in process against the GoodVibes CLI reading the same
output (WSL, -m real)."""

import json
import sys

import numpy as np
import pytest

from hfauto.backends.protocols import Capability
from hfauto.chemistry.symmetry import analyze
from hfauto.chemistry.thermo import species_thermo, thermal_modes
from hfauto.chemistry.xyz import XYZ, Molecule, read_xyz
from hfauto.core.method import MethodSpec, ThermoSettings
from hfauto.execution.process import Command, run_command

pytestmark = pytest.mark.real
PBE0 = MethodSpec(id="pbe0-d3bj_def2-svp", kind="dft", functional="pbe0", basis="def2-svp",
                  dispersion="d3bj", grid="fine", scf_energy_tol=1e-7)
H2O = XYZ(["O", "H", "H"], np.array([[0, 0, 0.12], [0, 0.76, -0.47], [0, -0.76, -0.47]]))


def test_in_process_thermo_matches_the_cli(real_engine, tmp_path):
    qm, run = real_engine(Capability.QM, "nwchem"), tmp_path / "run"
    opt = qm.optimize(Molecule(H2O, 0, 1), PBE0)
    freq = qm.frequencies(Molecule(read_xyz(run / opt.final.file.path), 0, 1), PBE0)
    xyz, hessian = read_xyz(run / freq.final.file.path), np.load(run / freq.hessian.path)
    sym = analyze(xyz.symbols, xyz.coords, hessian)
    assert (sym.point_group, sym.sigma) == ("C2v", 2)
    modes = thermal_modes(xyz.symbols, xyz.coords, hessian, linear=sym.linear, saddle=False)
    thermal = species_thermo(xyz.symbols, sym, modes, settings=ThermoSettings(), T=298.15)
    (tmp_path / "freq.out").write_bytes((run / freq.output.path).read_bytes())  # CLI: *.out
    argv = (sys.executable, "-m", "goodvibes", "freq.out", "-v", "1", "--zpe-vscal", "1",
            "--symm", "-q", "--json=gv")  # species_thermo is QH
    assert run_command(Command(argv=argv, cwd=tmp_path), timeout_s=600).returncode == 0
    cli = json.loads((tmp_path / "gv").read_text())["results"][0]["thermo"]
    # the CLI reads NWChem's printed frequencies and rotational constants
    assert abs(freq.energy_hartree + thermal.G - cli["qh_gibbs_free_energy"]) < 1e-5

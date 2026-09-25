"""Real GoodVibes on NWChem freq jobs (WSL, -m real): API vs CLI; S_rot > 0 for linear HNC."""

import json
import sys

import numpy as np
import pytest

from hfauto.backends.protocols import Capability
from hfauto.chemistry.xyz import XYZ, Molecule, read_xyz
from hfauto.core.method import MethodSpec, ThermoSettings
from hfauto.execution.process import Command, run_command

pytestmark = pytest.mark.real
PBE0 = MethodSpec(id="pbe0-d3bj_def2-svp", kind="dft", functional="pbe0", basis="def2-svp",
                  dispersion="d3bj", grid="fine", scf_energy_tol=1e-7)
UNSCALED = ThermoSettings(vib_scale=1.0, zpe_scale=1.0)
MOLECULES = {"h2o": (["O", "H", "H"], [[0, 0, 0.12], [0, 0.76, -0.47], [0, -0.76, -0.47]]),
             "hnc": (["H", "N", "C"], [[0, 0, -1.0], [0, 0, 0], [0, 0, 1.17]])}


@pytest.mark.parametrize(("name", "symbols", "coords"), [(k, *v) for k, v in MOLECULES.items()])
def test_api_thermo_matches_the_cli(real_engine, real_site, tmp_path, name, symbols, coords):
    qm, run = real_engine(Capability.QM, "nwchem"), tmp_path / "run"
    opt = qm.optimize(Molecule(XYZ(symbols, np.array(coords, dtype=float)), 0, 1), PBE0)
    freq = qm.frequencies(Molecule(read_xyz(run / opt.final.file.path), 0, 1), PBE0)
    [api] = real_engine(Capability.THERMO, "goodvibes").thermo(freq, [UNSCALED])
    assert api.S_rot > 0 and freq.n_external == (5 if name == "hnc" else 6)
    if name == "h2o":  # NWChem prints no finite A for a linear molecule, so the CLI reads H2O
        (tmp_path / "freq.out").write_bytes((run / freq.output.path).read_bytes())  # CLI: *.out
        argv = (real_site.engines["goodvibes"].python or sys.executable, "-m", "goodvibes",
                "freq.out", "-v", "1", "--zpe-vscal", "1", "--symm", "--json=gv")
        assert run_command(Command(argv=argv, cwd=tmp_path), timeout_s=600).returncode == 0
        cli = json.loads((tmp_path / "gv").read_text())["results"][0]["thermo"]
        assert abs(api.G_hartree - cli["qh_gibbs_free_energy"]) < 1e-5

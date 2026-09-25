"""Real xTB and pysisyphus (§10.3): run in WSL with HFAUTO_REAL=1 and HFAUTO_SITE."""

import numpy as np
import pytest

from hfauto.backends.protocols import Capability
from hfauto.backends.xtb import XTBEngine
from hfauto.chemistry.xyz import XYZ, Molecule, read_xyz
from hfauto.core.method import ExecutionSpec, MethodSpec

pytestmark = pytest.mark.real
GFN2 = MethodSpec(id="gfn2", kind="xtb", gfn=2)


def mol(symbols, *rows):
    return Molecule(XYZ(list(symbols), np.array(rows, dtype=float)), 0, 1)


def test_xtb_water_opt_hess_units_and_cycle_limit(real_engine, tmp_path):
    xtb, run = real_engine(Capability.QM, "xtb"), tmp_path / "run"
    opt = xtb.optimize(mol("OHH", [0, 0, 0], [1.05, 0, 0], [-0.3, 0.95, 0.05]), GFN2)
    freq = xtb.frequencies(Molecule(read_xyz(run / opt.final.file.path), 0, 1), GFN2)
    text = (run / freq.hessian.path).with_name("vibspectrum").read_text()
    rows = [r.split() for r in text.splitlines()]
    own = sorted(float(r[2]) for r in rows if len(r) == 5 and r[0].isdigit())  # xtb's own
    assert freq.n_external == 6 and np.allclose(freq.frequencies_cm1, own, rtol=2e-3)  # Eh/bohr²
    capped = XTBEngine(jobs=xtb.jobs, site=xtb.site.model_copy(
        update={"execution": ExecutionSpec(maxiter=2)}))  # --cycles 2
    far = mol("OHH", [0, 0, 0], [1.4, 0, 0], [1.0, 1.3, 0.4])
    assert capped.optimize(far, GFN2).kind == "geometry_maxiter"


def test_pysis_gs_hcn_to_hnc_ts_has_one_imaginary_mode(real_engine, tmp_path):
    hcn = mol("CNH", [0, 0, 0], [0, 0, 1.156], [0, 0, -1.066])  # exactly linear ends
    hnc = mol("CNH", [0, 0, 0], [0, 0, 1.17], [0, 0, 2.17])
    path = real_engine(Capability.PATH, "pysis_gs").find_path(hcn, hnc, GFN2, images=11,
                                                                refine_ts=True)
    assert len(path.energies_hartree) == 11 and path.ts is not None
    ts = Molecule(read_xyz(tmp_path / "run" / path.ts.file.path), 0, 1)
    freq = real_engine(Capability.QM, "xtb").frequencies(ts, GFN2)
    assert sum(f < -50 for f in freq.frequencies_cm1) == 1

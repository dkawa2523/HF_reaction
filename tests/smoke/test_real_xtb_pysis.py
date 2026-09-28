"""Real xTB and pysisyphus (§10.3): run in WSL with HFAUTO_REAL=1 and HFAUTO_SITE."""

import numpy as np
import pytest

from hfauto.backends.protocols import Capability
from hfauto.chemistry.identity import mapped_rmsd
from hfauto.chemistry.interpolation import idpp
from hfauto.chemistry.xyz import XYZ, Molecule, read_xyz, read_xyz_trajectory, write_xyz_trajectory
from hfauto.core.evidence import FileRef
from hfauto.core.hashing import sha256_file
from hfauto.core.method import MethodSpec

pytestmark = pytest.mark.real
GFN2 = MethodSpec(id="gfn2", kind="xtb", gfn=2)


def mol(symbols, *rows):
    return Molecule(XYZ(list(symbols), np.array(rows, dtype=float)), 0, 1)


def test_xtb_water_opt_hess_units(real_engine, tmp_path):
    xtb, run = real_engine(Capability.QM, "xtb"), tmp_path / "run"
    opt = xtb.optimize(mol("OHH", [0, 0, 0], [1.05, 0, 0], [-0.3, 0.95, 0.05]), GFN2)
    freq = xtb.frequencies(Molecule(read_xyz(run / opt.final.file.path), 0, 1), GFN2)
    text = (run / freq.hessian.path).with_name("vibspectrum").read_text()
    rows = [r.split() for r in text.splitlines()]
    own = sorted(float(r[2]) for r in rows if len(r) == 5 and r[0].isdigit())  # xtb's own
    assert freq.n_external == 6 and np.allclose(freq.frequencies_cm1, own, rtol=2e-3)  # Eh/bohr²


def test_pysis_neb_hcn_to_hnc_keeps_its_ends_and_its_ts_has_one_imaginary_mode(real_engine,
                                                                               tmp_path):
    """Fixed-end CI-NEB from hfauto's IDPP between exactly linear ends; TSOpt from the CI."""
    hcn = mol("CNH", [0, 0, 0], [0, 0, 1.156], [0, 0, -1.066])
    hnc = mol("CNH", [0, 0, 0], [0, 0, 1.17], [0, 0, 2.17])
    run, symbols = tmp_path / "run", list(hcn.xyz.symbols)
    frames = idpp(symbols, hcn.xyz.coords, hnc.xyz.coords, 11)
    initial = write_xyz_trajectory([XYZ(symbols, f) for f in frames], run / "idpp.xyz")
    ref = FileRef(path="idpp.xyz", sha256=sha256_file(initial))
    path = real_engine(Capability.PATH, "pysis_neb").find_path(hcn, hnc, GFN2, images=11,
                                                                initial_path=ref)
    images = read_xyz_trajectory(run / path.images.path)
    assert len(images) == 11 and path.ts is not None
    assert mapped_rmsd(images[0].coords, hcn.xyz.coords) < 1e-4
    assert mapped_rmsd(images[-1].coords, hnc.xyz.coords) < 1e-4
    ts = Molecule(read_xyz(run / path.ts.file.path), 0, 1)
    freq = real_engine(Capability.QM, "xtb").frequencies(ts, GFN2)
    assert sum(f < -50 for f in freq.frequencies_cm1) == 1

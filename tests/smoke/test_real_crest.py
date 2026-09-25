"""Real CREST (WSL, HFAUTO_REAL=1): --chrg / --uhf reach the program for NH3·HF."""

import numpy as np
import pytest

from hfauto.backends.protocols import Capability, ConformerEnsemble, ConformerSettings
from hfauto.chemistry.placement import seeds
from hfauto.chemistry.xyz import XYZ, Molecule
from hfauto.core.method import MethodSpec


@pytest.mark.real
def test_crest_receives_charge_and_spin(real_engine, tmp_path):
    nh3 = XYZ(["N", "H", "H", "H"], np.array([[0, 0, .1], [.94, 0, -.25], [-.47, .81, -.25],
                                              [-.47, -.81, -.25]]))
    seed = seeds(nh3, [XYZ(["H", "F"], np.array([[0, 0, 0], [0, 0, .92]]))], n_seeds=1)[0]
    crest, gfn2 = real_engine(Capability.CONFORMERS, "crest"), MethodSpec(id="g", kind="xtb", gfn=2)
    settings = ConformerSettings(nci=True, threads=2, topology="noref")
    neu, cat = (crest.search(Molecule(seed, q, m), gfn2, settings) for q, m in [(0, 1), (1, 2)])
    assert isinstance(neu, ConformerEnsemble) and neu.members, neu
    key = cat.job_key  # the radical cation may fail to optimize; its input must be right
    out = (tmp_path / "run" / "jobs" / key[:2] / key / "attempt_00" / "stdout.txt").read_text()
    flat = " ".join(out.split())
    assert "Molecular charge : 1" in flat and "UHF parameter : 1" in flat

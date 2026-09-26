"""Real CREST (WSL, HFAUTO_REAL=1): --chrg / --uhf reach the program for NH3·HF, and an anion
composition passes the initial topology check with --noopt (M4)."""

import numpy as np
import pytest

from hfauto.backends.protocols import Capability, ConformerEnsemble, ConformerSettings
from hfauto.chemistry.placement import seeds
from hfauto.chemistry.xyz import XYZ, Molecule
from hfauto.core.method import MethodSpec

HF = XYZ(["H", "F"], np.array([[0, 0, 0], [0, 0, .92]]))
GFN2 = MethodSpec(id="g", kind="xtb", gfn=2)


def _stdout(tmp_path, key):
    out = (tmp_path / "run" / "jobs" / key[:2] / key / "attempt_00" / "stdout.txt").read_text()
    return " ".join(out.split())


@pytest.mark.real
def test_crest_receives_charge_and_spin(real_engine, tmp_path):
    nh3 = XYZ(["N", "H", "H", "H"], np.array([[0, 0, .1], [.94, 0, -.25], [-.47, .81, -.25],
                                              [-.47, -.81, -.25]]))
    seed = seeds(nh3, [HF], n_seeds=1)[0]
    crest = real_engine(Capability.CONFORMERS, "crest")
    neu, cat = (crest.search(Molecule(seed, q, m), GFN2, ConformerSettings(nci=True))
                for q, m in [(0, 1), (1, 2)])
    assert isinstance(neu, ConformerEnsemble) and neu.members, neu
    flat = _stdout(tmp_path, cat.job_key)  # the radical cation may fail; its input must be right
    assert "Molecular charge : 1" in flat and "UHF parameter : 1" in flat


@pytest.mark.real
def test_anion_composition_runs_with_noopt(real_engine, tmp_path):
    """F-·(HF)2 stopped on the initial topology check without --noopt (about 25 s)."""
    seed = seeds(HF, [HF, XYZ(["F"], np.zeros((1, 3)))], n_seeds=1)[0]
    settings = ConformerSettings(nci=True, notopo_atoms=tuple(range(5)))  # every H and F
    ens = real_engine(Capability.CONFORMERS, "crest").search(Molecule(seed, -1, 1), GFN2, settings)
    assert isinstance(ens, ConformerEnsemble) and len(ens.members) >= 1, ens
    assert "--noopt" in _stdout(tmp_path, ens.job_key)

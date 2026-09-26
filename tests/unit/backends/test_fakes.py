"""Fakes on analytic surfaces: frequencies, saddle following, paths and scripts (§6.4)."""

import fakes
import numpy as np
import pytest

from hfauto.chemistry.profile import energies_settled, shape
from hfauto.chemistry.xyz import Molecule
from hfauto.core.evidence import FailureKind
from hfauto.core.method import MethodSpec

M = MethodSpec(id="pbe0", kind="dft", functional="pbe0", basis="def2-svp")


def test_harmonic_minimum_has_only_real_frequencies(tmp_run) -> None:
    p = fakes.harmonic()
    qm = fakes.FakeQM(tmp_run, p)
    opt = qm.optimize(p.molecule("start"), M)
    assert opt.trajectory_energies_hartree[0] == p.energy(p.points["start"]) > opt.energy_hartree
    freq = qm.frequencies(Molecule(fakes.xyz_loader(tmp_run)(opt.final), 0, 1), M)
    assert freq.start.fingerprint == opt.final.fingerprint and freq.level.program == "fake"
    assert freq.n_external == 6 and len(freq.frequencies_cm1) == 3 and min(freq.frequencies_cm1) > 0
    assert np.load(tmp_run / freq.hessian.path).shape == (9, 9)


def test_double_well_ts_has_one_imaginary_mode_and_saddle_finds_it(tmp_run) -> None:
    pes = fakes.double_well()
    qm, saddle = fakes.FakeQM(tmp_run, pes), fakes.FakeSaddle(tmp_run, pes)
    ts = qm.frequencies(pes.molecule("ts"), M)
    assert sum(nu < -50 for nu in ts.frequencies_cm1) == 1 and len(ts.imaginary_modes) == 1
    seed = pes.molecule("ts")
    seed.xyz.coords[1] += (0.05, 0.03, 0.0)
    found = saddle.refine(seed, M, hessian=qm.frequencies(seed, M))
    assert found.task == "saddle" and found.energy_hartree == pytest.approx(ts.energy_hartree)


@pytest.mark.parametrize("name,want", [("double_well", "single_max"), ("triple_well", "multi_max"),
                                       ("flat_uphill", "monotonic")])
def test_paths_follow_the_pes_or_the_script(tmp_run, name, want) -> None:
    pes = getattr(fakes, name)()
    path = fakes.FakePath(tmp_run, pes, script=["pes", "monotonic", "failed"])
    ends = pes.molecule("reactant"), pes.molecule("product")
    profile = path.find_path(*ends, M, images=11, refine_ts=want == "single_max")
    assert shape(profile.energies_hartree, 1e-4) == want
    # Like a real ZTS: gmax stays flat, the bead energies settle (not for "unconverged").
    assert energies_settled(profile.energy_history, 1e-4) and min(profile.gmax_history) > 1e-3
    if profile.ts is not None:
        assert profile.ts_energy_hartree == pytest.approx(pes.energy(pes.points["ts"]))
    assert shape(path.find_path(*ends, M, images=11).energies_hartree, 1e-4) == "monotonic"
    assert path.find_path(*ends, M, images=11).kind is FailureKind.NONZERO_EXIT
    stuck = fakes.FakePath(tmp_run, pes, script=["unconverged"]).find_path(*ends, M, images=11)
    assert not energies_settled(stuck.energy_history, 1e-3)

"""Fakes on analytic surfaces: frequencies, saddle following, paths and scripts (§6.4)."""

import fakes
import numpy as np
import pytest

from hfauto.chemistry.profile import classify
from hfauto.chemistry.xyz import XYZ, Molecule, read_xyz, write_xyz_trajectory
from hfauto.core.evidence import FailureKind, FileRef
from hfauto.core.hashing import sha256_file
from hfauto.core.method import MethodSpec

M = MethodSpec(id="pbe0", kind="dft", functional="pbe0", basis="def2-svp")


def test_harmonic_minimum_has_only_real_frequencies(tmp_run) -> None:
    p = fakes.harmonic()
    qm = fakes.FakeQM(tmp_run, p)
    opt = qm.optimize(p.molecule("start"), M)
    assert p.energy(p.points["start"]) > opt.energy_hartree
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
    hessian = qm.frequencies(seed, M)
    found = saddle.refine(seed, M, hessian=hessian, mode=hessian.imaginary_modes[0])
    assert found.task == "saddle" and found.energy_hartree == pytest.approx(ts.energy_hartree)


@pytest.mark.parametrize("name,want", [("double_well", "single"), ("triple_well", "intermediate"),
                                       ("flat_uphill", "barrierless")])
def test_paths_follow_the_pes_or_the_script(tmp_run, name, want) -> None:
    pes = getattr(fakes, name)()
    ends = pes.molecule("reactant"), pes.molecule("product")
    a, b = (m.xyz.coords for m in ends)
    frames = [XYZ(list(pes.symbols), a + t * (b - a)) for t in np.linspace(0.0, 1.0, 11)]
    initial = write_xyz_trajectory(frames, tmp_run / "initial.xyz")
    ref = FileRef(path="initial.xyz", sha256=sha256_file(initial))
    path = fakes.FakePath(tmp_run, pes, script=["pes", "barrierless", "failed"],
                          tsopt=want == "single")
    profile = path.find_path(*ends, M, images=11, initial_path=ref)  # the initial path's PES
    assert classify(profile.energies_hartree, 1e-4) == want
    if want == "single":  # like pysis_neb: a TS optimized from the highest image
        ts = read_xyz(tmp_run / profile.ts.file.path).coords
        assert pes.energy(ts) == pytest.approx(pes.energy(pes.points["ts"]))
    else:
        assert profile.ts is None
    again = path.find_path(*ends, M, images=11, initial_path=ref)
    assert classify(again.energies_hartree, 1e-4) == "barrierless"
    assert path.find_path(*ends, M, images=11, initial_path=ref).kind is FailureKind.NONZERO_EXIT

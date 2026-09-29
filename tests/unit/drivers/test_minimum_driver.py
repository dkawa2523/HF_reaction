"""MinimumDriver and Registry on analytic fake surfaces (§7.2)."""

from dataclasses import replace

import fakes
import numpy as np
import pytest
from scipy.spatial.distance import pdist

from hfauto.chemistry.xyz import composition_key
from hfauto.core.method import MethodSpec
from hfauto.core.records import SpeciesRecord
from hfauto.drivers.minimum import Registry, relax_to_minimum

M = MethodSpec(id="pbe0", kind="dft", functional="pbe0", basis="def2-svp")
CHFCLBR = np.array(
    [[0.0, 0.0, 0.0], [0.63, 0.63, 0.63], [-0.8, -0.8, 0.8], [-1.0, 1.0, -1.0], [1.1, -1.1, -1.1]]
)


def relax(root, pes, point, qm=None, **kw):
    kw.setdefault("load_xyz", fakes.xyz_loader(root))
    return relax_to_minimum(pes.molecule(point), M, qm or fakes.FakeQM(root, pes), **kw)


def add(registry, outcome, sid):
    g = outcome.opt.final
    return registry.add(outcome, SpeciesRecord(species_id=sid, state_label="x",
                        composition_id=composition_key(g.symbols, 0, 1), charge=0,
                        multiplicity=1, geometry=g, source="input"), tier="dft")


class SoftQM(fakes.FakeQM):
    def frequencies(self, mol, method, *, deadline=None):  # a spurious −30 cm⁻¹ mode at minima
        ev = super().frequencies(mol, method, deadline=deadline)
        if ev.imaginary_modes:
            return ev
        return ev.model_copy(update={"frequencies_cm1": (-30.0, *ev.frequencies_cm1[1:]),
                                     "imaginary_modes": (tuple(np.eye(9)[3]),)})


def test_harmonic_start_gives_a_minimum_or_with_a_stuck_soft_mode_a_soft_minimum(tmp_run):
    out = relax(tmp_run, fakes.harmonic(), "start")
    assert out.status == "minimum" and out.history == ("opt", "freq:none") and not out.ts_candidate
    assert out.freq.start.fingerprint == out.opt.final.fingerprint
    soft = relax(tmp_run, fakes.harmonic(), "start", SoftQM(tmp_run, fakes.harmonic()))
    assert soft.status == "soft_minimum" and soft.notes == ("soft_imaginary_mode",)
    assert soft.history[-1] == "soft:persisted"


def test_a_converged_opt_is_taken_as_it_is(tmp_run):
    pes = fakes.harmonic()
    qm = fakes.FakeQM(tmp_run, pes)
    opt = qm.optimize(pes.molecule("start"), M)
    out = relax(tmp_run, pes, "start", qm, opt=opt)
    assert out.status == "minimum" and out.history == ("opt:reused", "freq:none")
    assert out.opt is opt and qm.calls == ["optimize", "frequencies"]  # no second optimize


def springs(symbols, ref, points):
    """Pairwise springs at rest in ``ref``: every distance of ``ref`` is a minimum."""
    d0 = pdist(ref)
    return fakes.PES(symbols, lambda x: float(np.sum((pdist(np.reshape(x, (-1, 3))) - d0) ** 2)),
                     points)


def two_barrier_tops() -> fakes.PES:
    """Two F-H-F units, each H on top of its own double well (like the two methyl torsions of
    eclipsed C2v dimethyl ether): a second-order saddle. Springs hold the four F atoms."""
    first, second = fakes.symmetric_double_well(), fakes.symmetric_double_well(height=0.02)
    top = first.points["ts"]
    x0, heavy = np.vstack([top, top[:, ::-1] + [0.0, 3.0, 0.0]]), [0, 2, 3, 5]
    frame = springs(("F",) * 4, x0[heavy], {}).energy

    def energy(x):
        x = np.reshape(x, (-1, 3))
        return first.energy(x[:3]) + second.energy(x[3:]) + frame(x[heavy])

    return fakes.PES(("F", "H", "F") * 2, energy, {"top": x0})


@pytest.mark.parametrize("pes", [fakes.double_well(), fakes.symmetric_double_well()])
def test_a_first_order_saddle_goes_both_ways_from_its_own_hessian(tmp_run, pes) -> None:
    qm = fakes.FakeQM(tmp_run, pes)
    out = relax(tmp_run, pes, "ts", qm)
    assert out.status == "saddle" and out.history[-1] == "follow1:ts_candidate"
    ends = sorted(pes.energy(pes.points[p]) for p in ("reactant", "product"))
    assert sorted(side.opt.energy_hartree for side in out.ts_candidate) == pytest.approx(ends)
    for side in out.ts_candidate:  # each side's outcome, from its one opt and freq
        assert side.status == "minimum" and side.history == ("opt:follow1", "freq:none")
        assert side.freq.start.fingerprint == side.opt.final.fingerprint
    assert qm.calls == ["optimize", "frequencies", *["optimize+init_hessian", "frequencies"] * 2]


def test_a_soft_side_of_a_ts_candidate_gets_its_own_push(tmp_run) -> None:
    pes = fakes.double_well()
    out = relax(tmp_run, pes, "ts", SoftQM(tmp_run, pes))
    assert out.status == "saddle" and out.history[-1] == "follow1:ts_candidate"
    for side in out.ts_candidate:
        assert side.status == "soft_minimum" and side.notes == ("soft_imaginary_mode",)
        assert side.history == ("opt:follow1", "freq:soft", "soft:persisted")


@pytest.mark.parametrize("symbols", [("N", "H", "H", "H"), ("N", "H", "F", "Cl")])
def test_a_planar_amine_gives_an_inversion_ts_candidate(tmp_run, symbols) -> None:
    """The two pyramids are one basin, relabelled by a proper rotation (NH3) or mirror images
    (NHFCl): a degenerate TS between two structures."""
    pyramid = np.array([[0.0, 0.0, 0.38], [0.94, 0.0, 0.0], [-0.47, 0.814, 0.0],
                        [-0.47, -0.814, 0.0]])
    pes = springs(symbols, pyramid, {"planar": pyramid * [1.0, 1.0, 0.0]})
    out = relax(tmp_run, pes, "planar")
    assert out.status == "saddle" and out.history[-1] == "follow1:ts_candidate"
    plus, minus = (fakes.xyz_loader(tmp_run)(side.opt.final).coords for side in out.ts_candidate)
    assert plus[0, 2] * minus[0, 2] < 0  # N above and below the H3 plane


def test_a_second_order_saddle_descends_one_side_without_a_ts_candidate(tmp_run) -> None:
    pes = two_barrier_tops()
    qm = fakes.FakeQM(tmp_run, pes)
    out = relax(tmp_run, pes, "top", qm)
    assert out.history == ("opt", "freq:saddle", "follow1:one_side")
    assert out.status == "minimum" and out.ts_candidate is None
    assert qm.calls == ["optimize", "frequencies", "optimize+init_hessian", "frequencies"]
    x = fakes.xyz_loader(tmp_run)(out.opt.final).coords.reshape(2, 3, 3)
    off_centre = np.linalg.norm(x[:, 1] - 0.5 * (x[:, 0] + x[:, 2]), axis=1)
    assert (off_centre > 0.35).all()  # each H in a well of its own unit (0.4 Å off centre)


def test_registry_known_new_joined_and_init_hessian(tmp_run) -> None:
    pes, load = fakes.double_well(), fakes.xyz_loader(tmp_run)
    qm, registry = fakes.FakeQM(tmp_run, pes), Registry([], load)
    a, b = relax(tmp_run, pes, "reactant", qm), relax(tmp_run, pes, "product", qm)
    ra, rb = add(registry, a, "a"), add(registry, b, "b")
    assert ra.basin_id != rb.basin_id and ra.members == ("a",)
    hess = qm.frequencies(pes.molecule("reactant"), M)
    known = relax_to_minimum(pes.molecule("reactant"), M, qm, known=registry, init_hessian=hess)
    assert known.status == "known" and known.known_basin == ra.basin_id
    assert qm.calls[-1] == "optimize+init_hessian"  # used, and no freq job after it
    assert add(registry, known, "a2").basin_id == ra.basin_id
    calls = len(qm.calls)
    reused = relax_to_minimum(pes.molecule("reactant"), M, qm, known=registry, opt=b.opt)
    assert reused.known_basin == rb.basin_id and len(qm.calls) == calls  # no job at all
    image = fakes.xyz_loader(tmp_run)(b.opt.final).coords * [-1.0, 1.0, 1.0]
    assert registry.find(a.opt, coords=image) is None  # the coordinates, not the opt's final
    assert registry.find(b.opt, coords=image) == rb.basin_id  # a mirror image is one basin
    assert relax(tmp_run, pes, "product", qm, init_hessian=hess).failure.kind == "input_invalid"
    energy = ra.energy_hartree
    shifted = replace(a, opt=a.opt.model_copy(update={"energy_hartree": energy + 3e-5}))
    joined = add(registry, shifted, "c")  # one criterion: assign
    assert joined.basin_id == ra.basin_id and joined.members == ("a", "a2", "c")
    far = replace(a, opt=a.opt.model_copy(update={"energy_hartree": energy + 6e-5}))
    assert add(registry, far, "d").basin_id not in (ra.basin_id, rb.basin_id)
    rebuilt = Registry([(rb, b.opt.final)], load)  # (record, geometry)
    assert rebuilt.find(b.opt) == rb.basin_id and not rb.chiral


def test_mirror_images_share_one_chiral_basin_of_one_spin_state(tmp_run) -> None:
    pes = springs(("C", "H", "F", "Cl", "Br"), CHFCLBR,  # both enantiomers are exact minima
                  {"r": CHFCLBR, "s": CHFCLBR * [-1.0, 1.0, 1.0]})
    qm, registry = fakes.FakeQM(tmp_run, pes), Registry([], fakes.xyz_loader(tmp_run))
    r = add(registry, relax(tmp_run, pes, "r", qm), "r")
    s = relax_to_minimum(pes.molecule("s"), M, qm, known=registry)
    assert r.chiral and s.status == "known" and s.known_basin == r.basin_id
    triplet = relax_to_minimum(pes.molecule("s", multiplicity=3), M, qm, known=registry)
    assert triplet.status == "minimum"  # never known in the singlet basin (U3-I4)

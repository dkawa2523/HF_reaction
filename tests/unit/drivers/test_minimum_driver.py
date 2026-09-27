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
    def frequencies(self, mol, method, *, deadline=None):  # a spurious −30 cm⁻¹ mode
        ev = super().frequencies(mol, method, deadline=deadline)
        return ev.model_copy(update={"frequencies_cm1": (-30.0, *ev.frequencies_cm1[1:]),
                                     "imaginary_modes": (tuple(np.eye(9)[3]),)})


def test_harmonic_start_gives_a_minimum_or_with_a_stuck_soft_mode_a_soft_minimum(tmp_run):
    out = relax(tmp_run, fakes.harmonic(), "start")
    assert out.status == "minimum" and out.history == ("opt", "freq:none") and not out.ts_candidate
    assert out.freq.start.fingerprint == out.opt.final.fingerprint
    soft = relax(tmp_run, fakes.harmonic(), "start", SoftQM(tmp_run, fakes.harmonic()))
    assert soft.status == "soft_minimum" and soft.notes == ("soft_imaginary_mode",)
    assert soft.history[-1] == "soft:persisted"


@pytest.mark.parametrize("name,status,last", [("double_well", "saddle", "ts_candidate"),
                                              ("symmetric_double_well", "minimum", "replace")])
def test_mode_follow_from_a_transition_state(tmp_run, name, status, last) -> None:
    pes = getattr(fakes, name)()
    out = relax(tmp_run, pes, "ts")
    assert out.status == status and out.history[-1] == f"follow1:{last}"
    if status == "saddle":  # ± fell into two different minima: a TS candidate for free
        ends = sorted(pes.energy(pes.points[p]) for p in ("reactant", "product"))
        assert sorted(ev.energy_hartree for ev in out.ts_candidate) == pytest.approx(ends)


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
    assert relax(tmp_run, pes, "product", qm, init_hessian=hess).failure.kind == "input_invalid"
    energy = ra.energy_hartree
    shifted = replace(a, opt=a.opt.model_copy(update={"energy_hartree": energy + 3e-5}))
    assert add(registry, shifted, "c").basin_id == ra.basin_id  # one criterion: assign
    far = replace(a, opt=a.opt.model_copy(update={"energy_hartree": energy + 6e-5}))
    assert add(registry, far, "d").basin_id not in (ra.basin_id, rb.basin_id)
    assert registry.members(ra.basin_id) == ("a", "a2", "c")
    rebuilt = Registry([(rb, b.opt.final)], load)  # (record, geometry)
    assert rebuilt.find(b.opt) == rb.basin_id and not rb.chiral


def test_mirror_images_share_one_chiral_basin_of_one_spin_state(tmp_run) -> None:
    d0 = pdist(CHFCLBR)  # pairwise springs: both enantiomers are exact minima
    pes = fakes.PES(("C", "H", "F", "Cl", "Br"),
                    lambda x: float(np.sum((pdist(np.reshape(x, (-1, 3))) - d0) ** 2)),
                    {"r": CHFCLBR, "s": CHFCLBR * [-1.0, 1.0, 1.0]})
    qm, registry = fakes.FakeQM(tmp_run, pes), Registry([], fakes.xyz_loader(tmp_run))
    r = add(registry, relax(tmp_run, pes, "r", qm), "r")
    s = relax_to_minimum(pes.molecule("s"), M, qm, known=registry)
    assert r.chiral and s.status == "known" and s.known_basin == r.basin_id
    triplet = relax_to_minimum(pes.molecule("s", multiplicity=3), M, qm, known=registry)
    assert triplet.status == "minimum"  # never known in the singlet basin (U3-I4)

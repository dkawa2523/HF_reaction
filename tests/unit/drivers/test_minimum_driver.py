"""MinimumDriver and Registry on analytic fake surfaces (§7.2)."""

from dataclasses import replace

import fakes
import numpy as np
import pytest

from hfauto.chemistry.xyz import composition_key
from hfauto.core.method import MethodSpec
from hfauto.core.records import SpeciesRecord
from hfauto.drivers.minimum import Registry, relax_to_minimum

M = MethodSpec(id="pbe0", kind="dft", functional="pbe0", basis="def2-svp")


def relax(root, pes, point, qm=None, **kw):
    kw.setdefault("load_xyz", fakes.xyz_loader(root))
    return relax_to_minimum(pes.molecule(point), M, qm or fakes.FakeQM(root, pes), **kw)


def add(registry, outcome, sid):
    g = outcome.opt.final
    return registry.add(outcome, SpeciesRecord(species_id=sid, formula="x", state_label="x",
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


def test_registry_known_new_same_ambiguous_and_init_hessian(tmp_run) -> None:
    pes, load = fakes.double_well(), fakes.xyz_loader(tmp_run)
    qm, registry = fakes.FakeQM(tmp_run, pes), Registry([], load)
    a, b = relax(tmp_run, pes, "reactant", qm), relax(tmp_run, pes, "product", qm)
    (ra, va), (rb, vb) = add(registry, a, "a"), add(registry, b, "b")
    assert (va, vb) == ("new", "new") and ra.members == ("a",) and ra.n_fragments == 2
    hess = qm.frequencies(pes.molecule("reactant"), M)
    known = relax_to_minimum(pes.molecule("reactant"), M, qm, known=registry, init_hessian=hess)
    assert known.status == "known" and known.known_basin == ra.basin_id
    assert qm.calls[-1] == "optimize+init_hessian"  # used, and no freq job after it
    assert add(registry, known, "a2")[1] == "same"
    assert relax(tmp_run, pes, "product", qm, init_hessian=hess).failure.kind == "input_invalid"
    shifted = replace(a, opt=a.opt.model_copy(update={"energy_hartree": ra.energy_hartree + 3e-5}))
    rc, vc = add(registry, shifted, "c")
    assert vc == "ambiguous" and rc.basin_id not in (ra.basin_id, rb.basin_id)
    assert add(registry, a, "c")[1] == "same"  # re-judged after a (here: no) re-optimization
    assert registry.members(ra.basin_id) == ("a", "a2", "c")
    x, rebuilt = load(b.opt.final), Registry([(rb, b.opt.final)], load)  # (record, geometry)
    assert rebuilt.find(x.symbols, x.coords, rb.energy_hartree) == rb.basin_id

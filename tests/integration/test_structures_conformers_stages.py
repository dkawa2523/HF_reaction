"""structures → conformers on fake_runtime with FakeConformers (§8.2, CH-18, CH-23)."""

from functools import partial

import numpy as np
import pytest
from fakes import FakeConformers, write_geometry
from pydantic import ValidationError

from hfauto.backends.protocols import Capability, ConformerEnsemble
from hfauto.chemistry.xyz import XYZ, write_xyz
from hfauto.core.evidence import Failure, FailureKind
from hfauto.core.manifest import Manifest
from hfauto.core.method import MethodSpec
from hfauto.core.system import CompositionInput, ReactionInput, SpeciesInput, SystemConfig
from hfauto.stages.conformer_search import ConformersConfig, ConformersStage
from hfauto.stages.structures import StructuresConfig, StructuresStage

GEOMS = {"nh3": ("N H H H", [[0, 0, .1], [.94, 0, -.25], [-.47, .81, -.25], [-.47, -.81, -.25]]),
         "hf": ("H F", [[0, 0, 0], [0, 0, .92]]),
         "hono": ("H O N O", [[.95, .3, 0], [0, 0, 0], [-.5, 1.3, 0], [-1.7, 1.3, 0]]),
         "hcn": ("H C N", [[0, 0, -1.07], [0, 0, 0], [0, 0, 1.16]]),
         "cnh": ("C N H", [[0, 0, 0], [0, 0, 1.17], [0, 0, 2.16]])}
ENSEMBLE = partial(ConformerEnsemble, members=(), topology_stops=(), version="0", job_key="")
GFN2 = {"gfn2": MethodSpec(id="gfn2", kind="xtb", gfn=2)}


def _species(tmp_path, keys):
    return [SpeciesInput(id=k, xyz=write_xyz(XYZ(s.split(), np.array(x)), tmp_path / f"{k}.xyz"),
                         role="endpoint" if k in ("hcn", "cnh") else "monomer")
            for k, (s, x) in GEOMS.items() if k in keys]


def _proton_shift(x):  # HONO's H moved from O to N: CREST's topology stop
    return np.vstack([x[2] + [.3, .95, 0], x[1:]])


def _crest(rt, calls, n_stops):
    """Monomers stop ``n_stops`` times on the proton shift; HF·HF fails; else 2 members."""
    def search(mol, method, settings):
        x, sym = mol.xyz.coords, mol.xyz.symbols
        calls.append((mol, settings))
        key = f"f/{len(sym)}_{len(calls)}"
        if sym == ["H", "F", "H", "F"]:
            return Failure(kind=FailureKind.NONZERO_EXIT, reason="scripted")
        if not settings.nci and sum(not s.nci for _, s in calls) <= n_stops:
            stop = write_geometry(rt().run_dir, f"{key}s.xyz", sym, _proton_shift(x))
            return ENSEMBLE(topology_stops=(stop,))
        g = [write_geometry(rt().run_dir, f"{key}{i}.xyz", sym, x + i * 1e-3) for i in (0, 1)]
        return ENSEMBLE(members=((g[0], -10.0), (g[1], -9.99)))  # CREGEN-unique; +6.3 kcal/mol
    return FakeConformers(search)


def test_structures_then_conformers(tmp_path, fake_runtime):
    species = _species(tmp_path, GEOMS)
    system = SystemConfig(system_id="t", species=[*species, SpeciesInput(id="pair", smiles="N.F")],
                          compositions=[CompositionInput(id="nh", components={"nh3": 1, "hf": 1}),
                                        CompositionInput(id="hf2", components={"hf": 2})],
                          reactions=[ReactionInput(id="bad", reactant="hcn", product="cnh")])
    arts = StructuresStage().run(Manifest(run_id="r", stage_id="s", created_at=""),
                                 StructuresConfig(), fake_runtime(system, {}, stage_id="s"))
    failed = {a.artifact_id: a.failure.reason for a in arts if a.failure}
    assert set(failed) == {"species_pair", "species_hcn", "species_cnh"}
    assert "composition" in failed["species_pair"] and "atom order" in failed["species_hcn"]
    nh3 = next(a.payload for a in arts if a.artifact_id == "species_nh3")
    assert nh3.geometry.file.path == "s/xyz/nh3.xyz" and nh3.composition_id == "H3N_q0_m1"

    calls = []
    rt = fake_runtime(system, {(Capability.CONFORMERS, "crest"): _crest(lambda: rt, calls, 1)},
                      methods=GFN2, stage_id="c")
    with pytest.raises(ValidationError):  # no energy window, no engine settings
        ConformersConfig(engine="crest", method="gfn2", window_kcal=4.0)
    out = ConformersStage().run(Manifest(run_id="r", stage_id="s", created_at="", artifacts=arts),
                                ConformersConfig(engine="crest", method="gfn2", ewin_kcal=5.0), rt)
    got = sorted(a.payload.species_id for a in out if a.payload)  # tags: c/t/p = source
    assert [k for k in got if not k.startswith("hf2_p")] == [  # +6.3 kcal/mol kept (c01); the
        "hf_c00", "hono_c00", "hono_c01", "nh3_c00", "nh_c00", "nh_c01"]  # stop has no energy
    assert len(got) > 6 and [a.artifact_id for a in out if a.failure] == ["conformers_hf2"]
    (first, on), (retry, again) = calls[:2]  # HONO retries once from the stop, same settings
    assert np.allclose(retry.xyz.coords, _proton_shift(first.xyz.coords)) and again == on
    assert (retry.charge, retry.multiplicity) == (first.charge, first.multiplicity) and not on.nci
    compositions = {len(m.xyz.symbols): s for m, s in calls[2:]}  # searched concurrently
    assert all(s.nci and s.notopo_atoms == tuple(range(n)) for n, s in compositions.items())
    assert {s.ewin_kcal for _, s in calls} == {5.0} and set(compositions) == {4, 6}


def test_a_repeated_topology_stop_fails_and_the_stop_seeds_compositions(tmp_path, fake_runtime):
    system = SystemConfig(system_id="t", species=_species(tmp_path, ("hono", "hf")),
                          compositions=[CompositionInput(id="ohf", components={"hono": 1, "hf": 1})])
    arts = StructuresStage().run(Manifest(run_id="r", stage_id="s", created_at=""),
                                 StructuresConfig(), fake_runtime(system, {}, stage_id="s"))
    calls = []
    rt = fake_runtime(system, {(Capability.CONFORMERS, "crest"): _crest(lambda: rt, calls, 2)},
                      methods=GFN2, stage_id="c")
    out = ConformersStage().run(Manifest(run_id="r", stage_id="s", created_at="", artifacts=arts),
                                ConformersConfig(engine="crest", method="gfn2"), rt)
    failed = [(a.artifact_id, a.failure.kind, a.failure.reason) for a in out if a.failure]
    assert failed == [("conformers_hono", FailureKind.INCOMPLETE_OUTPUT, "topology_stop_repeated")]
    hono = {a.payload.species_id: a.payload.source for a in out
            if a.payload and a.payload.species_id.startswith("hono")}
    assert hono == {"hono_t00": "crest_topology", "hono_c01": "conformer"}  # the stop; the input
    assert len(calls) == 3 and calls[1][1] == calls[0][1] and calls[2][1].nci  # no second retry
    stop = _proton_shift(calls[0][0].xyz.coords)  # the lowest HONO seeds the composition
    assert np.allclose(calls[2][0].xyz.coords[:4], stop)


def test_a_composition_takes_the_declared_or_the_only_coupled_multiplicity(tmp_path, fake_runtime):
    oh = write_xyz(XYZ(["O", "H"], np.array([[0, 0, 0], [0, 0, .97]])), tmp_path / "oh.xyz")
    comps = [CompositionInput(id=k, components=c, multiplicity=m) for k, c, m in (
        ("free", {"oh": 2}, None), ("triplet", {"oh": 2}, 3), ("singlet", {"oh": 2}, 1),
        ("ohf", {"oh": 1, "hf": 1}, None))]
    species = [SpeciesInput(id="oh", xyz=oh, multiplicity=2), *_species(tmp_path, ("hf",))]
    system = SystemConfig(system_id="t", species=species, compositions=comps)
    arts = StructuresStage().run(Manifest(run_id="r", stage_id="s", created_at=""),
                                 StructuresConfig(), fake_runtime(system, {}, stage_id="s"))
    calls = []
    rt = fake_runtime(system, {(Capability.CONFORMERS, "crest"): _crest(lambda: rt, calls, 0)},
                      methods=GFN2, stage_id="c")
    out = ConformersStage().run(Manifest(run_id="r", stage_id="s", created_at="", artifacts=arts),
                                ConformersConfig(engine="crest", method="gfn2"), rt)
    failed = {a.artifact_id: a.failure.reason for a in out if a.failure}
    assert all(a.failure.kind == FailureKind.INPUT_INVALID for a in out if a.failure)
    assert set(failed) == {"conformers_free", "conformers_singlet"}
    assert failed["conformers_free"] == "declare_multiplicity: candidates (1, 3)"
    assert failed["conformers_singlet"].startswith("open_shell_singlet_unsupported")
    assert sorted(m.multiplicity for m, _ in calls) == [2, 3]  # OH·HF doublet, (OH)2 triplet

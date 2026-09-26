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
ENSEMBLE = partial(ConformerEnsemble, members=(), topology_removed=2, topology_stops=(),
                   version="0", job_key="")


def test_structures_then_conformers(tmp_path, fake_runtime):
    species = [SpeciesInput(id=k, xyz=write_xyz(XYZ(s.split(), np.array(x)), tmp_path / f"{k}.xyz"),
                            role="endpoint" if k in ("hcn", "cnh") else "monomer")
               for k, (s, x) in GEOMS.items()]
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

    def script(mol, method, settings):  # HONO stops once on a proton shift; HF·HF fails
        x, sym = mol.xyz.coords, mol.xyz.symbols
        key = f"f/{len(sym)}{settings.topology}"
        calls.append((len(sym), settings))
        if sym == ["H", "F", "H", "F"]:
            return Failure(kind=FailureKind.NONZERO_EXIT, reason="scripted")
        if settings.topology == "on" and not settings.nci:  # H moved from O to N
            stop = np.vstack([x[2] + [.3, .95, 0], x[1:]])
            return ENSEMBLE(topology_stops=(write_geometry(rt.run_dir, f"{key}s.xyz", sym, stop),))
        g = [write_geometry(rt.run_dir, f"{key}{i}.xyz", sym, x + i * 1e-3) for i in (0, 1)]
        return ENSEMBLE(members=((g[0], -10.0), (g[1], -10.0 + 1e-6), (g[0], -9.99)))  # dup; +6.3

    calls = []
    rt = fake_runtime(system, {(Capability.CONFORMERS, "crest"): FakeConformers(script)},
                      methods={"gfn2": MethodSpec(id="gfn2", kind="xtb", gfn=2)}, stage_id="c")
    with pytest.raises(ValidationError):  # no energy window, no engine settings
        ConformersConfig(engine="crest", method="gfn2", window_kcal=4.0)
    out = ConformersStage().run(Manifest(run_id="r", stage_id="s", created_at="", artifacts=arts),
                                ConformersConfig(engine="crest", method="gfn2", ewin_kcal=5.0), rt)
    got = sorted(a.payload.species_id for a in out if a.payload)  # tags: c/t/p = source
    assert [k for k in got if not k.startswith("hf2_p")] == [  # +6.3 kcal/mol kept (c01)
        "hf_c00", "hono_c00", "hono_c01", "hono_t02", "nh3_c00", "nh_c00", "nh_c01"]
    assert len(got) > 7 and [a.artifact_id for a in out if a.failure] == ["conformers_hf2"]
    assert [s.topology for _, s in calls[:2]] == ["on", "noref"] and not calls[0][1].nci
    compositions = dict(calls[2:])  # searched concurrently: keyed by atom count
    assert compositions[6].nci and compositions[6].notopo_atoms == (0, 1, 2, 3, 4, 5)
    assert compositions[4].nci and {s.ewin_kcal for _, s in calls} == {5.0}

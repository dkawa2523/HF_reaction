"""structures: declared-reaction endpoint checks (review U1-I2 b)."""

import numpy as np

from hfauto.chemistry.xyz import XYZ, write_xyz
from hfauto.core.evidence import FailureKind
from hfauto.core.manifest import Manifest
from hfauto.core.records import CoordinateTerm
from hfauto.core.system import ReactionInput, SpeciesInput, SystemConfig
from hfauto.stages.structures import StructuresConfig, StructuresStage


def test_endpoint_mismatch_and_coordinate_range_fail_both_ends(tmp_path, fake_runtime):
    hf = write_xyz(XYZ(["H", "F"], np.array([[0, 0, 0], [0, 0, 0.92]])), tmp_path / "hf.xyz")
    cation = {"charge": 1, "multiplicity": 2}
    species = [SpeciesInput(id=k, xyz=hf, role="endpoint", **(cation if k == "q" else {}))
               for k in ("p", "q", "r", "s", "u", "v")]
    bond = [CoordinateTerm(kind="distance", atoms=(0, 1))]
    reactions = [ReactionInput(id="charge", reactant="p", product="q", coordinate=bond),
                 ReactionInput(id="index", reactant="r", product="s",
                               coordinate=[CoordinateTerm(kind="distance", atoms=(0, 2))]),
                 ReactionInput(id="ok", reactant="u", product="v", coordinate=bond)]
    system = SystemConfig(system_id="t", species=species, reactions=reactions)
    arts = StructuresStage().run(Manifest(run_id="r", stage_id="s", created_at=""),
                                 StructuresConfig(), fake_runtime(system, {}, stage_id="s"))
    failed = {a.artifact_id: a.failure for a in arts if a.failure}
    assert set(failed) == {"species_p", "species_q", "species_r", "species_s"}
    assert all(f.kind == FailureKind.INPUT_INVALID for f in failed.values())
    assert "reaction charge: endpoints differ" in failed["species_q"].reason
    assert "atom order" in failed["species_p"].reason
    assert failed["species_r"].reason == "reaction index: coordinate atom index out of range"

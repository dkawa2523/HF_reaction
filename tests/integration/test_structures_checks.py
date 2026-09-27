"""structures: declared-reaction endpoint checks (review U1-I2 b)."""

import numpy as np
import pytest

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


def test_elements_and_spin_are_checked_once_at_the_entrance(tmp_path, fake_runtime):
    pytest.importorskip("rdkit")

    def xyz(name, symbols):
        return write_xyz(XYZ(symbols, np.eye(len(symbols), 3) * 2.2), tmp_path / f"{name}.xyz")

    fecl3 = xyz("fecl3", ["Fe", "Cl", "Cl", "Cl"])
    species = [SpeciesInput(id="ce", xyz=xyz("ce", ["Ce", "Cl", "Cl", "Cl"])),
               SpeciesInput(id="fr", xyz=xyz("fr", ["Fr"]), charge=1),
               SpeciesInput(id="fe", xyz=fecl3), SpeciesInput(id="fe6", xyz=fecl3, multiplicity=6),
               SpeciesInput(id="fe2", smiles="[Fe+2]", charge=2),
               SpeciesInput(id="o2", smiles="[O][O]")]
    system = SystemConfig(system_id="t", species=species)
    arts = StructuresStage().run(Manifest(run_id="r", stage_id="s", created_at=""),
                                 StructuresConfig(), fake_runtime(system, {}, stage_id="s"))
    failed = {a.artifact_id[8:]: a.failure.reason for a in arts if a.failure}
    assert all(a.failure.kind == FailureKind.INPUT_INVALID for a in arts if a.failure)
    assert failed == {"ce": "unsupported_element:Ce", "fr": "unsupported_element:Fr",
                      "fe": "declare_multiplicity: d-block element(s) Fe",
                      "fe2": "declare_multiplicity: d-block element(s) Fe"}
    ok = {a.payload.species_id: a.payload.composition_id for a in arts if a.payload}
    assert ok == {"fe6": "Cl3Fe_q0_m6", "o2": "O2_q0_m3"}  # [O][O]: radicals + 1

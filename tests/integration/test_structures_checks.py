"""structures: declared-reaction endpoint checks (review U1-I2 b)."""

import numpy as np
import pytest
from pydantic import ValidationError

from hfauto.chemistry.xyz import XYZ, write_xyz
from hfauto.core.evidence import FailureKind
from hfauto.core.manifest import Manifest
from hfauto.core.records import CoordinateTerm
from hfauto.core.system import ReactionInput, SpeciesInput, SystemConfig
from hfauto.stages.structures import StructuresConfig, StructuresStage


def test_endpoint_mismatch_and_coordinate_range_fail_both_ends(tmp_path, fake_runtime):
    hf = write_xyz(XYZ(["H", "F"], np.array([[0, 0, 0], [0, 0, 0.92]])), tmp_path / "hf.xyz")
    spin = {"q": {"charge": 1, "multiplicity": 2}, "x": {"multiplicity": 3}}
    species = [SpeciesInput(id=k, xyz=hf, role="endpoint",
                            **({"multiplicity": 1} | spin.get(k, {})))
               for k in ("p", "q", "r", "s", "u", "v", "w", "x")]
    bond = [CoordinateTerm(kind="distance", atoms=(0, 1))]
    reactions = [ReactionInput(id="charge", reactant="p", product="q", coordinate=bond),
                 ReactionInput(id="index", reactant="r", product="s",
                               coordinate=[CoordinateTerm(kind="distance", atoms=(0, 2))]),
                 ReactionInput(id="ok", reactant="u", product="v", coordinate=bond),
                 ReactionInput(id="spin", reactant="w", product="x", coordinate=bond)]
    system = SystemConfig(system_id="t", species=species, reactions=reactions)
    arts = StructuresStage().run(Manifest(run_id="r", stage_id="s", created_at=""),
                                 StructuresConfig(), fake_runtime(system, {}, stage_id="s"))
    failed = {a.artifact_id: a.failure for a in arts if a.failure}
    assert set(failed) == {f"species_{k}" for k in ("p", "q", "r", "s", "w", "x")}
    assert all(f.kind == FailureKind.INPUT_INVALID for f in failed.values())
    assert failed["species_q"].reason == "reaction charge: endpoints differ in charge or atom order"
    assert failed["species_r"].reason == "reaction index: coordinate atom index out of range"
    assert failed["species_w"].reason == failed["species_x"].reason == (  # only the spin differs
        "spin_crossing_reaction_unsupported: reaction spin joins multiplicities 1 and 3")


def test_the_declared_multiplicity_is_the_only_spin_source(tmp_path, fake_runtime):
    pytest.importorskip("rdkit")

    def xyz(name, symbols):
        return write_xyz(XYZ(symbols, np.eye(len(symbols), 3) * 2.2), tmp_path / f"{name}.xyz")

    fecl3 = xyz("fecl3", ["Fe", "Cl", "Cl", "Cl"])
    with pytest.raises(ValidationError, match="multiplicity"):  # no default, xyz or SMILES
        SpeciesInput(id="fe", xyz=fecl3)
    with pytest.raises(ValidationError, match="multiplicity"):
        SpeciesInput(id="ch3", smiles="[CH3]")
    species = [SpeciesInput(id="ce", xyz=xyz("ce", ["Ce", "Cl", "Cl", "Cl"]), multiplicity=2),
               SpeciesInput(id="fr", xyz=xyz("fr", ["Fr"]), charge=1, multiplicity=1),
               SpeciesInput(id="fe6", xyz=fecl3, multiplicity=6),
               SpeciesInput(id="cf2", smiles="F[C]F", multiplicity=1),  # singlet ground state
               SpeciesInput(id="ch3", smiles="[CH3]", multiplicity=1),
               SpeciesInput(id="o2", smiles="[O][O]", multiplicity=3)]
    system = SystemConfig(system_id="t", species=species)
    arts = StructuresStage().run(Manifest(run_id="r", stage_id="s", created_at=""),
                                 StructuresConfig(), fake_runtime(system, {}, stage_id="s"))
    failed = {a.artifact_id[8:]: a.failure.reason for a in arts if a.failure}
    assert all(a.failure.kind == FailureKind.INPUT_INVALID for a in arts if a.failure)
    assert failed == {"ce": "unsupported_element:Ce", "fr": "unsupported_element:Fr",
                      "ch3": "Electron-count parity is incompatible with multiplicity: "
                             "electrons=9, multiplicity=1"}
    ok = {a.payload.species_id: a.payload.composition_id for a in arts if a.payload}
    assert ok == {"fe6": "Cl3Fe_q0_m6", "cf2": "CF2_q0_m1", "o2": "O2_q0_m3"}

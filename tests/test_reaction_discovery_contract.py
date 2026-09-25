from __future__ import annotations

from pathlib import Path

import numpy as np

from hfauto.backends.reaction_discovery.base import DiscoveryResult
from hfauto.backends.reaction_discovery.readuct import _afir_bias_target
from hfauto.chemistry.xyz import XYZ, write_xyz
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.chemistry import (
    DrivingAtomPairRecord,
    ReactionTrialRecord,
)
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import StageContext
from hfauto.stages.explore_reactions import (
    ExploreReactionsStage,
    _select_product_endpoint,
)
from hfauto.stages.minimum_registry import MinimumRegistryStage
from hfauto.stages.reaction_plan import ReactionPlanStage


class _RecordingDiscoveryBackend:
    name = "recording"

    def __init__(self, product_xyz: str) -> None:
        self.product_xyz = product_xyz

    def explore(self, trial, source, settings, workdir, *, driver):
        return DiscoveryResult(
            driver=driver,
            success=True,
            product_xyz_path=self.product_xyz,
            low_level_ts_validated=True,
            low_level_irc_connected=True,
            imaginary_mode_count=1,
            biased_energy_used_as_barrier=False,
            paths={"product": self.product_xyz},
        )


class _CheckpointDiscoveryBackend:
    name = "checkpoint"

    def __init__(self) -> None:
        self.calls = 0

    def explore(self, trial, source, settings, workdir, *, driver):
        self.calls += 1
        raw = Path(workdir) / "raw_attempt.txt"
        raw.parent.mkdir(parents=True, exist_ok=True)
        raw.write_text("no product\n", encoding="utf-8")
        return DiscoveryResult(
            driver=driver,
            success=False,
            failure_reason="bounded_search_no_product",
            paths={"raw": str(raw)},
        )


def test_afir_mixed_trial_biases_association_and_retains_dissociation_evidence() -> None:
    trial = ReactionTrialRecord(
        trial_id="mixed",
        source_species_id="reactant",
        associations=[DrivingAtomPairRecord(atoms=(1, 2))],
        dissociations=[DrivingAtomPairRecord(atoms=(2, 3))],
        driver_order=["afir"],
    )

    pair, attractive, evidence = _afir_bias_target(trial)

    assert pair == trial.associations[0]
    assert attractive is True
    assert evidence == {
        "afir_bias_target": "association",
        "directly_biased_pair": [1, 2],
        "unbiased_dissociation_pair_count": 1,
    }


def _write_hcn_pair(tmp_path: Path) -> tuple[Path, Path]:
    reactant = write_xyz(
        XYZ(
            ["C", "N", "H"],
            np.asarray([[0.0, 0.0, 0.0], [1.16, 0.0, 0.0], [-1.06, 0.0, 0.0]]),
        ),
        tmp_path / "hcn.xyz",
    )
    product = write_xyz(
        XYZ(
            ["C", "N", "H"],
            np.asarray([[0.0, 0.0, 0.0], [1.18, 0.0, 0.0], [2.18, 0.0, 0.0]]),
        ),
        tmp_path / "hnc.xyz",
    )
    return reactant, product


def test_explore_stage_emits_candidate_but_not_reaction(
    tmp_path: Path, monkeypatch
) -> None:
    reactant_path, product_path = _write_hcn_pair(tmp_path)
    source = Artifact(
        artifact_id="spc_hcn",
        artifact_type="species_preopt",
        paths={"xyz": str(reactant_path)},
        data={
            "species_id": "spc_hcn",
            "state": "candidate",
            "charge": 0,
            "multiplicity": 1,
            "xyz_path": str(reactant_path),
        },
    )
    trial = ReactionTrialRecord(
        trial_id="trial_hcn_hnc",
        source_species_id="spc_hcn",
        associations=[DrivingAtomPairRecord(atoms=(1, 2))],
        dissociations=[DrivingAtomPairRecord(atoms=(0, 2))],
        driver_order=["nt2"],
        max_attempts=1,
    )
    trial_artifact = Artifact(
        artifact_id=trial.trial_id,
        artifact_type="reaction_trial",
        parents=[source.artifact_id],
        data=trial.model_dump(mode="json"),
    )
    manifest = Manifest.new(run_id="test", stage="generate-reactions")
    manifest.extend([source, trial_artifact])
    monkeypatch.setattr(
        "hfauto.stages.explore_reactions.get_reaction_discovery_engine",
        lambda *_args, **_kwargs: _RecordingDiscoveryBackend(str(product_path)),
    )

    result = ExploreReactionsStage().run(
        manifest,
        {"engine": "recording"},
        StageContext(tmp_path / "explore", "test", {}),
    )

    assert len(result.latest_artifacts("reaction_candidate")) == 1
    assert result.latest_artifacts("reaction") == []
    attempt = result.latest_artifacts("reaction_discovery_attempt")[0]
    assert attempt.qc["biased_energy_used_as_barrier"] is False
    assert attempt.data["low_level_ts_validated"] is True
    summary = result.find("reaction_discovery_summary")
    assert summary is not None
    assert summary.data["outcome"] == "reaction_candidates_found"
    assert summary.qc["path_search_eligible"] is True


def test_explore_stage_does_not_rerun_completed_unbiased_relaxation(
    tmp_path: Path, monkeypatch
) -> None:
    reactant_path, _product_path = _write_hcn_pair(tmp_path)
    source = Artifact(
        artifact_id="spc_hcn",
        artifact_type="species_preopt",
        paths={"xyz": str(reactant_path)},
        data={
            "species_id": "spc_hcn",
            "state": "candidate",
            "charge": 0,
            "multiplicity": 1,
            "xyz_path": str(reactant_path),
        },
    )
    trial = ReactionTrialRecord(
        trial_id="trial_relaxed",
        source_species_id="spc_hcn",
        associations=[DrivingAtomPairRecord(atoms=(1, 2))],
        driver_order=["unbiased_preopt"],
        max_attempts=1,
    )
    data = trial.model_dump(mode="json")
    data["discovery_complete"] = True
    manifest = Manifest.new(run_id="relaxed", stage="relaxation-discovery")
    manifest.extend(
        [
            source,
            Artifact(
                artifact_id=trial.trial_id,
                artifact_type="reaction_trial",
                parents=[source.artifact_id],
                data=data,
            ),
        ]
    )
    backend = _CheckpointDiscoveryBackend()
    monkeypatch.setattr(
        "hfauto.stages.explore_reactions.get_reaction_discovery_engine",
        lambda *_args, **_kwargs: backend,
    )

    ExploreReactionsStage().run(
        manifest,
        {"engine": "checkpoint", "stop_after_first_candidate": False},
        StageContext(tmp_path / "explore", "relaxed", {}),
    )

    assert backend.calls == 0


def test_explore_stage_resumes_integrity_checked_attempt_checkpoint(
    tmp_path: Path, monkeypatch
) -> None:
    reactant_path, _ = _write_hcn_pair(tmp_path)
    source = Artifact(
        artifact_id="spc_hcn",
        artifact_type="species_preopt",
        paths={"xyz": str(reactant_path)},
        data={
            "species_id": "spc_hcn",
            "state": "candidate",
            "charge": 0,
            "multiplicity": 1,
            "xyz_path": str(reactant_path),
        },
    )
    trial = ReactionTrialRecord(
        trial_id="trial_checkpoint",
        source_species_id="spc_hcn",
        associations=[DrivingAtomPairRecord(atoms=(1, 2))],
        dissociations=[DrivingAtomPairRecord(atoms=(0, 2))],
        driver_order=["nt2"],
        max_attempts=1,
    )
    manifest = Manifest.new(run_id="checkpoint", stage="generate-reactions")
    manifest.extend(
        [
            source,
            Artifact(
                artifact_id=trial.trial_id,
                artifact_type="reaction_trial",
                parents=[source.artifact_id],
                data=trial.model_dump(mode="json"),
            ),
        ]
    )
    backend = _CheckpointDiscoveryBackend()
    monkeypatch.setattr(
        "hfauto.stages.explore_reactions.get_reaction_discovery_engine",
        lambda *_args, **_kwargs: backend,
    )
    context = StageContext(
        tmp_path / "explore",
        "checkpoint",
        {"mode": "production"},
    )

    ExploreReactionsStage().run(manifest, {"engine": "checkpoint"}, context)
    resumed = ExploreReactionsStage().run(
        manifest, {"engine": "checkpoint"}, context
    )

    assert backend.calls == 1
    assert resumed.metadata["reaction_discovery_attempt_count"] == 1
    attempt = resumed.latest_artifacts("reaction_discovery_attempt")[0]
    assert attempt.data["request_fingerprint"]
    assert attempt.data["path_hashes"]


def test_bond_driven_trial_rejects_large_same_graph_relaxation(
    tmp_path: Path,
) -> None:
    source = write_xyz(
        XYZ(
            ["O", "H", "H"],
            np.asarray([[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]]),
        ),
        tmp_path / "water_source.xyz",
    )
    relaxed = write_xyz(
        XYZ(
            ["O", "H", "H"],
            np.asarray([[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.75, 0.60, 0.0]]),
        ),
        tmp_path / "water_relaxed.xyz",
    )
    trial = ReactionTrialRecord(
        trial_id="trial_water",
        source_species_id="water",
        associations=[DrivingAtomPairRecord(atoms=(1, 2))],
        driver_order=["afir"],
    )

    path, evidence = _select_product_endpoint(
        trial,
        str(source),
        DiscoveryResult(
            driver="afir", success=True, product_xyz_path=str(relaxed)
        ),
        same_graph_product_rmsd_A=0.01,
    )

    assert path is None
    assert evidence["changed_graph"] is False
    assert evidence["structural_change_detected"] is False
    assert evidence["reason"] == "requested_bond_changes_not_realized"


def test_coordinate_trial_requires_declared_internal_coordinate_progress(
    tmp_path: Path,
) -> None:
    source = write_xyz(
        XYZ(
            ["H", "O", "N", "O"],
            np.asarray(
                [
                    [-0.20, -0.95, 0.0],
                    [0.00, 0.00, 0.0],
                    [1.43, 0.00, 0.0],
                    [1.83, 1.11, 0.0],
                ]
            ),
        ),
        tmp_path / "trans_hono.xyz",
    )
    product = write_xyz(
        XYZ(
            ["H", "O", "N", "O"],
            np.asarray(
                [
                    [-0.20, 0.95, 0.0],
                    [0.00, 0.00, 0.0],
                    [1.43, 0.00, 0.0],
                    [1.83, 1.11, 0.0],
                ]
            ),
        ),
        tmp_path / "cis_hono.xyz",
    )
    trial = ReactionTrialRecord(
        trial_id="trial_hono_torsion",
        source_species_id="trans_hono",
        driver_order=["dft_minimum_ensemble"],
        reaction_coordinate={
            "min_change": 1.0,
            "terms": [
                {
                    "kind": "dihedral",
                    "atoms": [0, 1, 2, 3],
                    "coefficient": 1.0,
                }
            ],
        },
    )

    path, evidence = _select_product_endpoint(
        trial,
        str(source),
        DiscoveryResult(
            driver="dft_minimum_ensemble",
            success=True,
            product_xyz_path=str(product),
        ),
        same_graph_product_rmsd_A=10.0,
    )

    assert path == str(product)
    assert evidence["structural_change_basis"] == (
        "declared_reaction_coordinate_change"
    )
    assert evidence["reaction_coordinate_evidence"]["accepted"] is True


def test_bond_driven_trial_rejects_equivalent_atom_exchange(
    tmp_path: Path,
) -> None:
    source = write_xyz(
        XYZ(
            ["H", "F", "H", "F"],
            np.asarray(
                [[0.0, 0.0, 0.0], [0.92, 0.0, 0.0], [5.0, 0.0, 0.0], [5.92, 0.0, 0.0]]
            ),
        ),
        tmp_path / "two_hf_source.xyz",
    )
    exchanged = write_xyz(
        XYZ(
            ["H", "F", "H", "F"],
            np.asarray(
                [[0.0, 0.0, 0.0], [5.92, 0.0, 0.0], [5.0, 0.0, 0.0], [0.92, 0.0, 0.0]]
            ),
        ),
        tmp_path / "two_hf_exchanged.xyz",
    )
    trial = ReactionTrialRecord(
        trial_id="trial_hf_exchange",
        source_species_id="two_hf",
        associations=[DrivingAtomPairRecord(atoms=(0, 3))],
        dissociations=[DrivingAtomPairRecord(atoms=(0, 1))],
        driver_order=["nt2"],
    )

    path, evidence = _select_product_endpoint(
        trial,
        str(source),
        DiscoveryResult(
            driver="nt2",
            success=True,
            product_xyz_path=str(exchanged),
        ),
        same_graph_product_rmsd_A=0.01,
    )

    assert path is None
    assert evidence["atom_mapped_graph_changed"] is True
    assert evidence["permutation_invariant_graph_isomorphic"] is True
    assert evidence["changed_graph"] is False
    assert evidence["structural_change_detected"] is False
    assert evidence["reason"] == "equivalent_atom_permutation_only"


def _minimum(
    artifact_id: str,
    species_id: str,
    xyz_path: Path,
    energy: float,
) -> tuple[Artifact, Artifact]:
    calculation_id = f"calc_{species_id}"
    species = Artifact(
        artifact_id=artifact_id,
        artifact_type="species_optimized",
        parents=[calculation_id],
        paths={"xyz": str(xyz_path)},
        data={
            "species_id": species_id,
            "state": "candidate",
            "charge": 0,
            "multiplicity": 1,
            "xyz_path": str(xyz_path),
        },
        method={"stage": "dft-minima", "engine": "nwchem"},
        provenance={
            "validated_qm_method": {
                "engine": "nwchem",
                "functional": "pbe0",
                "basis": "def2-svpd",
                "disp_vdw": 3,
                "program_version": "7.2.3",
            }
        },
        qc={
            "is_minimum": True,
            "minimum_accepted": True,
            "scf_converged": True,
            "geometry_converged": True,
            "geometry_sane": True,
            "fallback_dummy": False,
            "real_qm_executed": True,
        },
    )
    calculation = Artifact(
        artifact_id=calculation_id,
        artifact_type="calculation",
        parents=[artifact_id],
        paths={"final_xyz": str(xyz_path)},
        data={
            "species_id": species_id,
            "task": "opt_freq",
            "n_imag": 0,
            "frequency_count_complete": True,
            "electronic_energy_hartree": energy,
        },
        method={
            "stage": "dft-minima",
            "engine": "nwchem",
            "functional": "pbe0",
            "basis": "def2-svpd",
            "disp_vdw": 3,
            "program_version": "7.2.3",
        },
        qc={
            "is_minimum": True,
            "minimum_accepted": True,
            "frequency_count_complete": True,
            "real_qm_executed": True,
            "fallback_dummy": False,
        },
    )
    return species, calculation


def test_registry_promotes_only_distinct_dft_minima(tmp_path: Path) -> None:
    reactant_path, product_path = _write_hcn_pair(tmp_path)
    reactant, reactant_calc = _minimum(
        "opt_hcn", "spc_hcn", reactant_path, -93.0
    )
    product, product_calc = _minimum(
        "opt_hnc", "spc_hnc", product_path, -92.98
    )
    candidate = Artifact(
        artifact_id="candidate_hcn_hnc",
        artifact_type="reaction_candidate",
        data={
            "candidate_id": "candidate_hcn_hnc",
            "trial_id": "trial_hcn_hnc",
            "reactant_species_id": "spc_hcn",
            "product_species_id": "spc_hnc",
            "driver": "nt2",
            "bond_changes": [
                {"kind": "break", "atoms": [0, 2]},
                {"kind": "form", "atoms": [1, 2]},
            ],
            "composition_preserved": True,
            "charge": 0,
            "multiplicity": 1,
            "low_level_endpoint_distinct": True,
            "evidence_artifact_id": "attempt_hcn_hnc",
            "mechanism_family": "isomerization",
            "reaction_coordinate": {
                "terms": [
                    {"kind": "distance", "atoms": [0, 2], "coefficient": 1.0},
                    {"kind": "distance", "atoms": [1, 2], "coefficient": -1.0},
                ]
            },
        },
    )
    manifest = Manifest.new(run_id="test", stage="dft-minima")
    manifest.extend(
        [candidate, reactant, reactant_calc, product, product_calc]
    )

    registered = MinimumRegistryStage().run(
        manifest,
        {"require_real_qm": True},
        StageContext(tmp_path / "registry", "test", {}),
    )
    planned = ReactionPlanStage().run(
        registered,
        {},
        StageContext(tmp_path / "plan", "test", {}),
    )

    assert len(planned.latest_artifacts("minimum_basin")) == 2
    reactions = planned.latest_artifacts("reaction")
    assert len(reactions) == 1
    assert reactions[0].qc["distinct_registry_basins"] is True
    assessment = reactions[0].data["basin_assessment"]
    assert assessment["reactant_evidence"]["electronic_energy_hartree"] == -93.0
    assert assessment["product_evidence"]["electronic_energy_hartree"] == -92.98
    case = planned.latest_artifacts("reaction_case")[0]
    assert case.data["state"] == "path_search_ready"
    assert case.data["ts_search_allowed"] is True


def test_registry_same_basin_does_not_promote_reaction(tmp_path: Path) -> None:
    reactant_path, _ = _write_hcn_pair(tmp_path)
    reactant, reactant_calc = _minimum(
        "opt_hcn", "spc_hcn", reactant_path, -93.0
    )
    duplicate, duplicate_calc = _minimum(
        "opt_hcn_duplicate", "spc_hcn_duplicate", reactant_path, -93.0
    )
    candidate = Artifact(
        artifact_id="candidate_same",
        artifact_type="reaction_candidate",
        data={
            "candidate_id": "candidate_same",
            "trial_id": "trial_same",
            "reactant_species_id": "spc_hcn",
            "product_species_id": "spc_hcn_duplicate",
            "driver": "nt2",
            "bond_changes": [{"kind": "form", "atoms": [1, 2]}],
            "composition_preserved": True,
            "charge": 0,
            "multiplicity": 1,
            "low_level_endpoint_distinct": True,
            "evidence_artifact_id": "attempt_same",
        },
    )
    manifest = Manifest.new(run_id="test", stage="dft-minima")
    manifest.extend(
        [candidate, reactant, reactant_calc, duplicate, duplicate_calc]
    )
    registered = MinimumRegistryStage().run(
        manifest,
        {"require_real_qm": True},
        StageContext(tmp_path / "registry_same", "test", {}),
    )
    planned = ReactionPlanStage().run(
        registered,
        {},
        StageContext(tmp_path / "plan_same", "test", {}),
    )

    assert len(planned.latest_artifacts("minimum_basin")) == 1
    assert planned.latest_artifacts("reaction") == []
    case = planned.latest_artifacts("reaction_case")[0]
    assert case.data["state"] == "same_basin"
    assert case.data["scientific_conclusion_supported"] is True

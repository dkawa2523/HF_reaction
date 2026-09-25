from pathlib import Path

import numpy as np

from hfauto.chemistry.connectivity import connectivity_changes
from hfauto.chemistry.xyz import XYZ, write_xyz
from hfauto.chemistry.xyz_trajectory import write_xyz_trajectory
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import StageContext
from hfauto.stages.preopt import PreoptStage
from hfauto.stages.relaxation_discovery import RelaxationDiscoveryStage


def _endpoint_pair() -> tuple[XYZ, XYZ]:
    initial = XYZ(
        ["N", "H", "F"],
        np.asarray([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.9, 0.0, 0.0]]),
    )
    final = XYZ(
        ["N", "H", "F"],
        np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.9, 0.0, 0.0]]),
        "Etot=-10.5",
    )
    return initial, final


def _source(tmp_path: Path, initial: XYZ) -> Artifact:
    path = write_xyz(initial, tmp_path / "source.xyz")
    return Artifact(
        artifact_id="spc_base_hf",
        artifact_type="species",
        paths={"xyz": str(path)},
        data={
            "species_id": "spc_base_hf",
            "state": "encounter_complex",
            "xyz_path": str(path),
            "charge": 0,
            "multiplicity": 1,
            "components": [
                {"component_id": "base", "atom_indices": [0]},
                {"component_id": "hf", "atom_indices": [1, 2]},
            ],
        },
    )


def test_connectivity_change_detects_proton_transfer_graph() -> None:
    initial, final = _endpoint_pair()

    result = connectivity_changes(initial, final)

    assert result["accepted"] is True
    assert result["topology_changed"] is True
    assert result["formed_bonds"] == [[0, 1]]
    assert result["broken_bonds"] == [[1, 2]]


def test_preopt_topology_change_becomes_candidate_without_barrier_claim(
    tmp_path: Path,
) -> None:
    initial, final = _endpoint_pair()
    source = _source(tmp_path, initial)
    final_path = write_xyz(final, tmp_path / "final.xyz")
    geometry = Artifact(
        artifact_id="preopt_geom_spc_base_hf",
        artifact_type="preopt_geometry",
        parents=[source.artifact_id, "calc_preopt"],
        paths={"xyz": str(final_path), "source_xyz": source.paths["xyz"]},
        data={"electronic_energy_hartree": -10.5},
    )
    manifest = Manifest.new(run_id="relax", stage="preopt")
    manifest.extend([source, geometry])

    result = RelaxationDiscoveryStage().run(
        manifest,
        {"include_crest_topology_stops": False},
        StageContext(tmp_path / "out", "relax", {"mode": "production"}),
    )

    candidate = result.latest_artifacts("reaction_candidate")[-1]
    trial = result.find(candidate.data["trial_id"])
    attempt = result.find(candidate.data["evidence_artifact_id"])
    product = result.find(candidate.data["product_species_id"])
    assert trial is not None and trial.data["driver_order"] == ["unbiased_preopt"]
    assert attempt is not None and attempt.data["barrier_claimed"] is False
    assert product is not None and product.data["state"] == "product_candidate"
    assert candidate.data["effectively_barrierless_low_level_hypothesis"] is True
    assert candidate.qc["activation_barrier_validated"] is False


def test_production_preopt_does_not_overwrite_a_connectivity_changed_state(
    tmp_path: Path, monkeypatch
) -> None:
    initial, final = _endpoint_pair()
    source = _source(tmp_path, initial)
    final_path = write_xyz(final, tmp_path / "optimized.xyz")

    class Engine:
        def optimize_frequency(self, _species, _method, _workdir):
            return Artifact(
                artifact_id="calc_preopt",
                artifact_type="calculation",
                paths={"final_xyz": str(final_path)},
                data={"electronic_energy_hartree": -10.5},
                method={"engine": "test"},
                qc={"geometry_sane": True, "fallback_dummy": False},
            )

    monkeypatch.setattr("hfauto.stages.preopt.get_qm_engine", lambda *_a, **_k: Engine())
    manifest = Manifest.new(run_id="preopt", stage="build-complexes")
    manifest.add_artifact(source)

    result = PreoptStage().run(
        manifest,
        {"states": ["encounter_complex"], "update_species": True},
        StageContext(tmp_path / "out", "preopt", {"mode": "production"}),
    )

    calculation = result.find("calc_preopt")
    geometry = result.find("preopt_geom_spc_base_hf")
    assert calculation is not None
    assert calculation.qc["connectivity_change"]["topology_changed"] is True
    assert calculation.qc["preopt_promotion_accepted"] is False
    assert geometry is not None and geometry.status.status == "partial"
    assert result.latest_artifacts("species_preopt") == []


def test_preopt_resumes_only_hash_verified_real_calculation(
    tmp_path: Path, monkeypatch
) -> None:
    initial, _final = _endpoint_pair()
    source = _source(tmp_path, initial)
    calls: list[Path] = []

    class Engine:
        def optimize_frequency(self, _species, _method, workdir):
            root = Path(workdir)
            calls.append(root)
            final_path = write_xyz(initial, root / "optimized.xyz")
            return Artifact(
                artifact_id="calc_preopt_resume",
                artifact_type="calculation",
                paths={"final_xyz": str(final_path)},
                data={"electronic_energy_hartree": -10.0},
                method={"engine": "xtb"},
                qc={
                    "geometry_sane": True,
                    "fallback_dummy": False,
                    "real_qm_executed": True,
                },
            )

    monkeypatch.setattr("hfauto.stages.preopt.get_qm_engine", lambda *_a, **_k: Engine())
    manifest = Manifest.new(run_id="resume", stage="build-complexes")
    manifest.add_artifact(source)
    context = StageContext(tmp_path / "out", "resume", {"mode": "production"})
    settings = {"engine": "xtb", "states": ["encounter_complex"]}

    PreoptStage().run(manifest, settings, context)
    result = PreoptStage().run(manifest, settings, context)

    assert [path.name for path in calls] == ["attempt_00"]
    assert result.metadata["preopt_resumed_jobs"] == 1
    calculation = result.find("calc_preopt_resume")
    assert calculation is not None
    assert calculation.provenance["reused_from_checkpoint"] is True

    Path(calculation.paths["final_xyz"]).write_text(
        "3\ncorrupted checkpoint\nN 0 0 0\nH 1.1 0 0\nF 2.9 0 0\n",
        encoding="utf-8",
    )
    PreoptStage().run(manifest, settings, context)
    assert [path.name for path in calls] == ["attempt_00", "attempt_01"]


def test_preopt_can_select_only_crest_nci_sources(
    tmp_path: Path, monkeypatch
) -> None:
    initial, _final = _endpoint_pair()
    seed = _source(tmp_path, initial)
    nci = seed.model_copy(deep=True)
    nci.artifact_id = "spc_base_hf_nci_0"
    nci.data.update(
        {
            "species_id": nci.artifact_id,
            "crest_nci_conformer_id": "spc_base_hf_conf0000",
        }
    )
    calls: list[str] = []

    class Engine:
        def optimize_frequency(self, species, _method, workdir):
            calls.append(species.artifact_id)
            final_path = write_xyz(initial, Path(workdir) / "optimized.xyz")
            return Artifact(
                artifact_id=f"calc_{species.artifact_id}",
                artifact_type="calculation",
                paths={"final_xyz": str(final_path)},
                method={"engine": "xtb"},
                qc={
                    "geometry_sane": True,
                    "fallback_dummy": False,
                    "real_qm_executed": True,
                },
            )

    monkeypatch.setattr(
        "hfauto.stages.preopt.get_qm_engine", lambda *_a, **_k: Engine()
    )
    manifest = Manifest.new(run_id="select_nci", stage="build-complexes")
    manifest.extend([seed, nci])

    result = PreoptStage().run(
        manifest,
        {
            "engine": "xtb",
            "states": ["encounter_complex"],
            "required_source_data_keys": ["crest_nci_conformer_id"],
        },
        StageContext(
            tmp_path / "selected",
            "select_nci",
            {"mode": "production"},
        ),
    )

    assert calls == ["spc_base_hf_nci_0"]
    assert result.metadata["preopt_required_source_data_keys"] == [
        "crest_nci_conformer_id"
    ]


def test_crest_topology_stop_trajectory_uses_the_same_candidate_contract(
    tmp_path: Path,
) -> None:
    initial, final = _endpoint_pair()
    source = _source(tmp_path, initial)
    raw = tmp_path / "crest"
    trajectory = write_xyz_trajectory(
        [initial, final], raw / "crestopt.log"
    )
    stdout = raw / "crest.stdout"
    stdout.write_text(
        "Geometry successfully optimized.\nChange in topology detected\n",
        encoding="utf-8",
    )
    failure = Artifact.failure(
        "crest_failed_spc_base_hf",
        "conformer",
        "CREST stopped after topology change",
        category="crest_failed",
        parents=[source.artifact_id],
        command_result={"cwd": str(raw), "stdout_path": str(stdout)},
    )
    manifest = Manifest.new(run_id="crest", stage="build-complexes")
    manifest.extend([source, failure])

    result = RelaxationDiscoveryStage().run(
        manifest,
        {"include_crest_topology_stops": True},
        StageContext(tmp_path / "out", "crest", {"mode": "production"}),
    )

    candidate = result.latest_artifacts("reaction_candidate")[-1]
    attempt = result.find(candidate.data["evidence_artifact_id"])
    assert attempt is not None
    assert attempt.data["backend"] == "crest_initial_optimization"
    assert attempt.paths["optimization_trajectory"] == str(trajectory)

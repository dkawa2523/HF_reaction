from pathlib import Path

from hfauto.chemistry.minimum_connections import assess_trial_between_minima
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.chemistry import ReactionTrialRecord
from hfauto.core.schemas.manifest import Manifest
from hfauto.stages.base import StageContext
from hfauto.stages.connect_minima import ConnectMinimaStage


def _write_hcn_hnc(tmp_path: Path) -> tuple[Path, Path]:
    hcn = tmp_path / "hcn.xyz"
    hnc = tmp_path / "hnc.xyz"
    hcn.write_text(
        "3\nHCN\nH 0.00 0 0\nC 1.06 0 0\nN 2.22 0 0\n",
        encoding="utf-8",
    )
    hnc.write_text(
        "3\nHNC\nH 3.22 0 0\nC 1.06 0 0\nN 2.22 0 0\n",
        encoding="utf-8",
    )
    return hcn, hnc


def _trial() -> ReactionTrialRecord:
    return ReactionTrialRecord(
        trial_id="trial_hcn_hnc",
        source_species_id="hcn",
        associations=[{"atoms": [0, 2], "label": "H-N"}],
        dissociations=[{"atoms": [0, 1], "label": "H-C"}],
        driver_order=["nt2"],
        mechanism_hint="isomerization",
        priority=10.0,
        charge=0,
        multiplicity=1,
    )


def test_minimum_connection_accepts_forward_and_reverse_trial(tmp_path: Path) -> None:
    hcn, hnc = _write_hcn_hnc(tmp_path)

    forward = assess_trial_between_minima(_trial(), hcn, hnc)
    reverse = assess_trial_between_minima(_trial(), hnc, hcn)

    assert forward["accepted"] is True
    assert forward["direction"] == "forward"
    assert reverse["accepted"] is True
    assert reverse["direction"] == "reverse"
    assert {change["kind"] for change in reverse["bond_changes"]} == {
        "form",
        "break",
    }


def test_threshold_crossing_without_large_distance_change_is_coordinate_only(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.xyz"
    target = tmp_path / "target.xyz"
    source.write_text(
        "3\nsource\nF 0 0 0\nH 1.04 0 0\nN 2.28 0 0\n",
        encoding="utf-8",
    )
    target.write_text(
        "3\ntarget\nF 0 0 0\nH 1.08 0 0\nN 2.28 0 0\n",
        encoding="utf-8",
    )
    trial = ReactionTrialRecord(
        trial_id="trial_small_reorganization",
        source_species_id="source",
        associations=[{"atoms": [1, 2], "label": "H-N"}],
        dissociations=[{"atoms": [0, 1], "label": "H-F"}],
        driver_order=["afir"],
        charge=0,
        multiplicity=1,
    )

    result = assess_trial_between_minima(trial, source, target)

    assert result["accepted"] is True
    assert result["connection_kind"] == "coordinate_reorganization"
    assert result["bond_changes"] == []
    assert 0.05 < result["forward_coordinate_progress_A"] < 0.15


def test_explicit_dihedral_trial_connects_conformational_minima(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.xyz"
    target = tmp_path / "target.xyz"
    source.write_text(
        "4\nammonia up\nN 0 0 0.36\nH 0.94 0 0\n"
        "H -0.47 0.81406388 0\nH -0.47 -0.81406388 0\n",
        encoding="utf-8",
    )
    target.write_text(
        "4\nammonia down\nN 0 0 -0.36\nH 0.94 0 0\n"
        "H -0.47 0.81406388 0\nH -0.47 -0.81406388 0\n",
        encoding="utf-8",
    )
    trial = ReactionTrialRecord(
        trial_id="trial_inversion",
        source_species_id="ammonia",
        driver_order=["nt2"],
        reaction_coordinate={
            "min_change": 0.2,
            "terms": [
                {
                    "kind": "dihedral",
                    "atoms": [0, 1, 2, 3],
                    "coefficient": 1.0,
                }
            ],
        },
    )

    result = assess_trial_between_minima(trial, source, target)

    assert result["accepted"] is True
    assert result["connection_kind"] == "coordinate_reorganization"
    assert result["bond_changes"] == []
    assert result["reaction_coordinate"]["terms"]


def test_connect_minima_publishes_a_dft_level_candidate(tmp_path: Path) -> None:
    hcn, hnc = _write_hcn_hnc(tmp_path)
    manifest = Manifest.new(run_id="connect", stage="minimum-registry")
    for species_id, path in (("hcn", hcn), ("hnc", hnc)):
        manifest.add_artifact(
            Artifact(
                artifact_id=f"opt_{species_id}",
                artifact_type="species_optimized",
                paths={"xyz": str(path)},
                data={
                    "species_id": species_id,
                    "state": "minimum",
                    "xyz_path": str(path),
                    "charge": 0,
                    "multiplicity": 1,
                },
            )
        )
    signature = {
        "composition": {"C": 1, "H": 1, "N": 1},
        "charge": 0,
        "multiplicity": 1,
        "method": {"method": {"engine": "nwchem"}},
    }
    for basin_id, species_id in (("basin_hcn", "hcn"), ("basin_hnc", "hnc")):
        manifest.add_artifact(
            Artifact(
                artifact_id=basin_id,
                artifact_type="minimum_basin",
                data={
                    "basin_id": basin_id,
                    "representative_species_id": species_id,
                    "signature": signature,
                    "unresolved_against": [],
                },
            )
        )
    manifest.add_artifact(
        Artifact(
            artifact_id="minimum_registry",
            artifact_type="minimum_registry",
            data={
                "species_to_basin": {"hcn": "basin_hcn", "hnc": "basin_hnc"}
            },
        )
    )
    manifest.add_artifact(
        Artifact(
            artifact_id="trial_hcn_hnc",
            artifact_type="reaction_trial",
            data=_trial().model_dump(mode="json"),
        )
    )

    output = ConnectMinimaStage().run(
        manifest,
        {"max_connections": 2},
        StageContext(out_dir=tmp_path / "out", run_id="connect", global_config={}),
    )
    candidates = output.latest_artifacts("reaction_candidate")

    assert len(candidates) == 1
    assert candidates[0].data["driver"] == "dft_minimum_ensemble"
    assert candidates[0].data["endpoint_evidence_level"] == "dft_minimum_ensemble"
    assert candidates[0].qc["distinct_registry_basins"] is True

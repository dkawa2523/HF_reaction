from pathlib import Path

import pandas as pd

from hfauto.core.config import load_yaml
from hfauto.core.io import read_manifest
from hfauto.reporting.html_report import latest_manifest_path
from hfauto.workflow.runner import run_pipeline


def test_dummy_pipeline(tmp_path):
    cfg = load_yaml(Path(__file__).parents[1] / "configs" / "pipelines" / "dummy.yaml")
    cfg["run_root"] = str(tmp_path / "runs")
    cfg["input"]["sdf"] = str(Path(__file__).parents[1] / "examples" / "candidates.sdf")
    run_dir = run_pipeline(cfg, "pytest_demo")
    manifest = read_manifest(latest_manifest_path(run_dir))
    assert manifest.stage == "rank"
    assert any(a.artifact_type == "ranking" for a in manifest.artifacts)
    assert any(a.artifact_type == "thermo" for a in manifest.artifacts)
    assert any(a.artifact_type == "reaction_validated" for a in manifest.artifacts)

    rank_dir = Path(latest_manifest_path(run_dir)).parent
    for name in ["rank_scavenger.csv", "rank_activation.csv", "candidate_summary.csv", "failure_report.csv"]:
        assert (rank_dir / name).exists(), name

    scav = pd.read_csv(rank_dir / "rank_scavenger.csv")
    assert {"delta_G_assoc_kcal_mol", "delta_G_ionpair_kcal_mol", "delta_G_act_kcal_mol", "quality_tier"}.issubset(scav.columns)
    assert scav["delta_G_assoc_kcal_mol"].notna().any()

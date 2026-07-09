from pathlib import Path

from hfauto.core.config import load_yaml
from hfauto.core.io import read_manifest
from hfauto.reporting.html_report import latest_manifest_path, render_report
from hfauto.workflow.runner import run_pipeline


def test_phase2_local_outputs(tmp_path):
    root = Path(__file__).parents[1]
    cfg = load_yaml(root / "configs" / "pipelines" / "phase2_local.yaml")
    cfg["run_root"] = str(tmp_path / "runs")
    cfg["input"]["sdf"] = str(root / "examples" / "candidates.sdf")
    run_dir = run_pipeline(cfg, "phase2_test")
    manifest = read_manifest(latest_manifest_path(run_dir))

    assert manifest.stage == "rank"
    ranking_ids = {a.artifact_id for a in manifest.iter_artifacts("ranking")}
    assert {"rank_scavenger", "rank_activation", "cluster_risk"} <= ranking_ids
    table_ids = {a.artifact_id for a in manifest.iter_artifacts("table")}
    assert {"candidate_summary", "reaction_results", "failure_report"} <= table_ids
    assert any(a.artifact_type == "thermo" and a.data.get("delta_G_assoc_kcal_mol") is not None for a in manifest.artifacts)
    assert any(a.artifact_type == "descriptor" and a.data.get("r_HF_A") is not None for a in manifest.artifacts)

    report_path = render_report(run_dir, tmp_path / "report.html")
    assert report_path.exists()
    assert "hfauto run report" in report_path.read_text(encoding="utf-8")

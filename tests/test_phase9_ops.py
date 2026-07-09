from __future__ import annotations

from pathlib import Path

import pandas as pd

from hfauto.core.config import load_yaml
from hfauto.core.io import read_manifest
from hfauto.reporting.html_report import latest_manifest_path
from hfauto.workflow.runner import run_pipeline
from hfauto_ops.bundle import build_operations_bundle
from hfauto_ops.core.compare import compare_runs
from hfauto_ops.core.resource import build_resource_plan, minutes_to_slurm_time


def test_phase9_resource_plan_contract(tmp_path: Path):
    root = Path(__file__).parents[1]
    cfg_path = root / "configs" / "pipelines" / "phase9_hpc_ops_offline.yaml"
    df = build_resource_plan(cfg_path, run_dir=None, ops_config={})
    assert {"stage", "ncores", "memory_gb", "time_slurm", "command"} <= set(df.columns)
    assert "dft-minima" in set(df["stage"])
    assert df.loc[df["stage"] == "ts-search", "memory_gb"].max() >= 32
    assert minutes_to_slurm_time(1500).startswith("1-")


def test_phase9_ops_bundle_on_existing_run(tmp_path: Path):
    root = Path(__file__).parents[1]
    cfg = load_yaml(root / "configs" / "pipelines" / "phase7_full_publicdb_rank_offline.yaml")
    cfg["run_root"] = str(tmp_path / "runs")
    run_dir = Path(run_pipeline(cfg, "ops_base"))
    out = Path(latest_manifest_path(run_dir)).parent
    paths = build_operations_bundle(run_dir, out, pipeline_config=root / "configs" / "pipelines" / "phase9_hpc_ops_offline.yaml", ops_config={"scheduler": "slurm"})
    assert paths["artifact_csv"].exists()
    assert paths["stage_csv"].exists()
    assert paths["resource_csv"].exists()
    assert paths["retry_csv"].exists()
    assert paths["slurm_submit_all"].exists()
    assert paths["operations_report"].exists()
    artifact_index = pd.read_csv(paths["artifact_csv"])
    assert {"artifact_id", "artifact_type", "status", "method_id"} <= set(artifact_index.columns)
    retry_plan = pd.read_csv(paths["retry_csv"])
    assert "retry_stage" in retry_plan.columns


def test_phase9_pipeline_ops_stage(tmp_path: Path):
    cfg = load_yaml("configs/pipelines/phase9_hpc_ops_offline.yaml")
    cfg["run_root"] = str(tmp_path / "runs")
    # Keep the integration test light while still exercising the stage bridge.
    for stage in cfg["stages"]:
        if stage.get("name") == "viz":
            stage["settings"] = {"max_molecule_dossiers": 1, "max_reaction_dossiers": 1, "top_energy_profiles": 1, "library_mode": "none", "asset_mode": "none"}
        if stage.get("name") == "ops":
            stage["pipeline_config"] = "configs/pipelines/phase9_hpc_ops_offline.yaml"
    run_dir = Path(run_pipeline(cfg, "phase9_ops"))
    manifest = read_manifest(latest_manifest_path(run_dir))
    assert manifest.stage == "ops"
    assert any(a.artifact_type == "operations_bundle" for a in manifest.artifacts)
    assert (Path(latest_manifest_path(run_dir)).parent / "operations_report.html").exists()
    assert (Path(latest_manifest_path(run_dir)).parent / "slurm" / "submit_all.sh").exists()


def test_phase9_compare_runs(tmp_path: Path):
    root = Path(__file__).parents[1]
    cfg = load_yaml(root / "configs" / "pipelines" / "phase7_full_publicdb_rank_offline.yaml")
    cfg["run_root"] = str(tmp_path / "runs")
    run_a = Path(run_pipeline(cfg, "run_a"))
    run_b = Path(run_pipeline(cfg, "run_b"))
    paths = compare_runs(run_a, run_b, tmp_path / "compare")
    assert paths["candidate_comparison"].exists()
    assert paths["reaction_comparison"].exists()
    assert paths["html"].exists()
    cand = pd.read_csv(paths["candidate_comparison"])
    assert "mol_id" in cand.columns

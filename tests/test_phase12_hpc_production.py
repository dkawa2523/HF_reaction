from pathlib import Path

import pandas as pd

from hfauto.core.config import load_yaml
from hfauto.core.io import read_manifest
from hfauto.reporting.html_report import latest_manifest_path
from hfauto.workflow.runner import run_pipeline
from hfauto_ops.schedulers.portable import render_scheduler_bundle


def test_phase12_pipeline_ops_outputs(tmp_path: Path):
    cfg = load_yaml("configs/pipelines/phase12_hpc_production_offline.yaml")
    cfg["run_root"] = str(tmp_path / "runs")
    for stage in cfg["stages"]:
        if stage.get("name") == "viz":
            stage["settings"] = {"max_molecule_dossiers": 1, "max_reaction_dossiers": 1, "top_energy_profiles": 1, "library_mode": "none", "asset_mode": "none"}
    run_dir = Path(run_pipeline(cfg, "phase12_ops"))
    manifest = read_manifest(latest_manifest_path(run_dir))
    assert manifest.stage == "ops"
    assert manifest.find("ops_phase12_bundle") is not None
    ops_dir = Path(latest_manifest_path(run_dir)).parent
    for rel in [
        "array_job_plan.csv",
        "array_items.jsonl",
        "scheduler/slurm/submit_all.sh",
        "scheduler/pbs/submit_all.sh",
        "scheduler/lsf/submit_all.sh",
        "slurm/submit_all.sh",
        "reuse_plan.csv",
        "qcarchive/qcarchive_records.jsonl",
        "backend_comparison_dashboard.html",
        "resource_autotune_report.html",
        "operations_report.html",
    ]:
        assert (ops_dir / rel).exists(), rel
    arr = pd.read_csv(ops_dir / "array_job_plan.csv")
    assert {"array_index", "stage", "target_artifact_id", "target_manifest", "cache_key"} <= set(arr.columns)
    assert len(arr) > 0


def test_phase12_portable_scheduler_bundle_contract(tmp_path: Path):
    plan = pd.DataFrame([
        {"stage": "dft-minima", "stage_index": 0, "command": "echo dft", "ncores": 16, "memory_gb": 32, "time_slurm": "24:00:00", "queue": "long"},
        {"stage": "ts-search", "stage_index": 1, "command": "echo ts", "ncores": 32, "memory_gb": 64, "time_slurm": "2-00:00:00", "queue": "long"},
    ])
    for scheduler in ["slurm", "pbs", "lsf", "local"]:
        paths = render_scheduler_bundle(plan, tmp_path / scheduler, scheduler=scheduler, run_id="demo", project_root="$PWD", prefix="stage")
        assert paths["submit_all"].exists()
        assert paths["status_command"].exists()

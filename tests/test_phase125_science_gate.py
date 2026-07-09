from pathlib import Path

import pandas as pd

from hfauto.core.config import load_yaml
from hfauto.core.qc import cap_confidence, production_rank_eligible, scientific_rank_eligible
from hfauto.reporting.html_report import latest_manifest_path
from hfauto.workflow.runner import run_pipeline


def test_confidence_caps_and_gate_helpers():
    assert cap_confidence(0.9, "Q1") == 0.30
    assert cap_confidence(0.9, "Q4", has_dummy_or_fallback=True) == 0.20
    assert scientific_rank_eligible("Q2") is True
    assert scientific_rank_eligible("Q1") is False
    assert production_rank_eligible("Q4", production_thermo_ready=True, real_irc_executed=True) is True
    assert production_rank_eligible("Q4", production_thermo_ready=False, real_irc_executed=True) is False


def test_phase125_rank_outputs_and_science_gate(tmp_path: Path):
    cfg = load_yaml("configs/pipelines/phase125_cleanup_offline.yaml")
    cfg["run_root"] = str(tmp_path / "runs")
    for stage in cfg["stages"]:
        if stage.get("name") == "viz":
            stage["settings"] = {"max_molecule_dossiers": 1, "max_reaction_dossiers": 1, "top_energy_profiles": 1, "library_mode": "none", "asset_mode": "none"}
    run_dir = Path(run_pipeline(cfg, "phase125_gate"))
    rank_dir = latest_manifest_path(run_dir).parent
    # latest manifest is ops, rank is two stages before viz/ops in the stable pipeline.
    rank_candidates = list(run_dir.glob("*_rank/candidate_summary.csv"))
    assert rank_candidates, "candidate_summary.csv not produced"
    rank_dir = rank_candidates[-1].parent
    for name in ["rank_screening.csv", "rank_scientific.csv", "rank_production.csv", "ranking_summary.json"]:
        assert (rank_dir / name).exists(), name
    summary = pd.read_csv(rank_dir / "candidate_summary.csv")
    assert {"scientific_rank_eligible", "production_rank_eligible", "science_next_action", "ops_next_action"} <= set(summary.columns)
    # Offline dummy run should be screening only, with capped confidence.
    assert not summary["production_rank_eligible"].fillna(False).astype(bool).any()
    assert summary["confidence_score"].max() <= 0.30 + 1e-12
    prod = pd.read_csv(rank_dir / "rank_production.csv")
    assert len(prod) == 0

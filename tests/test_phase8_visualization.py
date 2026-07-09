from __future__ import annotations

import json
from pathlib import Path

from hfauto.core.config import load_yaml
from hfauto.workflow.runner import run_pipeline
from hfauto_viz.core.run_loader import RunData, load_run
from hfauto_viz.reports.run_report import write_run_report, render_run_visual_report
from hfauto_viz.reports.molecule_dossier import render_molecule_dossier
from hfauto_viz.reports.reaction_dossier import render_reaction_dossier
from hfauto_viz.renderers.graphviz_network import render_reaction_network
from hfauto_viz.renderers.plotly_energy import energy_points


def test_phase8_energy_points_contract():
    row = {
        "delta_G_assoc_pressure_corrected_kcal_mol": -10.0,
        "delta_G_assoc_standard_kcal_mol": -15.0,
        "delta_G_act_kcal_mol": 4.0,
        "delta_G_ionpair_kcal_mol": -2.0,
    }
    pts = energy_points(row, mode="process")
    assert [p["state"] for p in pts] == ["B + (HF)n", "B···(HF)n", "TS", "BH+···F(HF)n-1−"]
    assert pts[1]["G"] == -10.0
    assert pts[2]["G"] == -6.0
    assert pts[3]["G"] == -12.0
    labels, values = energy_points(row, process=False)
    assert labels[1] == "B···(HF)n"
    assert values[1] == -15.0


def _make_run(tmp_path: Path, run_id: str = "viz_run") -> Path:
    cfg = load_yaml("configs/pipelines/phase7_full_publicdb_rank_offline.yaml")
    cfg["run_root"] = str(tmp_path / "runs")
    return Path(run_pipeline(cfg, run_id))


def test_phase8_visual_report_generates_expected_assets(tmp_path: Path):
    run_dir = _make_run(tmp_path, "viz_report")
    run = load_run(run_dir)
    out = run.run_dir / "14_viz"
    report = render_run_visual_report(run, out, top_n=2)
    assert report.exists()
    html = report.read_text(encoding="utf-8")
    assert "Phase 8 visualization report" in html
    assert (out / "screening" / "energy_profile_top_candidates.html").exists()
    assert (out / "screening" / "scavenger_vs_activation.html").exists()
    assert (out / "networks" / "reaction_network.dot").exists()
    assert (out / "chemiscope" / "hf_screening.json.gz").exists()
    vm = json.loads((out / "viz_manifest.json").read_text(encoding="utf-8"))
    assert vm["schema_version"] == "hfauto_viz.manifest.v1"
    assert any(v["viz_type"] == "run_report" for v in vm["visualizations"])


def test_phase8_write_run_report_compatibility_api(tmp_path: Path):
    run_dir = _make_run(tmp_path, "viz_compat")
    report = write_run_report(run_dir, out_dir=run_dir / "14_viz", top_n=2, asset_mode="none")
    assert report.exists()
    assert list((run_dir / "14_viz" / "molecules").glob("*/dossier.html"))
    assert list((run_dir / "14_viz" / "reactions").glob("*/dossier.html"))


def test_phase8_molecule_and_reaction_dossiers(tmp_path: Path):
    run_dir = _make_run(tmp_path, "viz_dossier")
    run = RunData(run_dir)
    candidate = run.read_table("candidate_summary")
    mol_id = str(candidate.iloc[0]["mol_id"])
    mol_paths = render_molecule_dossier(run, mol_id, run.run_dir / "14_viz" / "molecules")
    assert mol_paths["dossier"].exists()
    assert mol_paths["structure_2d"].exists()
    rr = run.read_table("reaction_results")
    reaction_id = str(rr.iloc[0]["reaction_id"])
    rxn_paths = render_reaction_dossier(run, reaction_id, run.run_dir / "14_viz" / "reactions")
    assert rxn_paths["dossier"].exists()
    assert rxn_paths["energy"].exists()


def test_phase8_graphviz_network_writes_dot(tmp_path: Path):
    run_dir = _make_run(tmp_path, "viz_network")
    run = RunData(run_dir)
    paths = render_reaction_network(run.read_table("candidate_summary"), run.read_table("reaction_results"), run.run_dir / "14_viz" / "networks", max_reactions=4)
    assert paths["dot"].exists()
    assert "digraph" in paths["dot"].read_text(encoding="utf-8")
    assert paths["html"].exists()


def test_phase8_pipeline_viz_stage_generates_bundle(tmp_path: Path):
    cfg = load_yaml("configs/pipelines/phase8_full_viz_offline.yaml")
    cfg["run_root"] = str(tmp_path / "runs")
    # Keep this integration test light while still exercising the stage bridge.
    for stage in cfg["stages"]:
        if stage.get("name") == "viz":
            stage["settings"] = {
                "max_molecule_dossiers": 1,
                "max_reaction_dossiers": 1,
                "top_energy_profiles": 1,
                "library_mode": "none",
            }
    run_dir = Path(run_pipeline(cfg, "viz_stage"))
    latest = json.loads((run_dir / "15_viz" / "manifest.json").read_text(encoding="utf-8"))
    assert latest["stage"] == "viz"
    assert (run_dir / "15_viz" / "report.html").exists()
    assert (run_dir / "15_viz" / "viz_manifest.json").exists()
    assert any(a["artifact_type"] == "visualization_bundle" for a in latest["artifacts"])

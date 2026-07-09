from pathlib import Path

from hfauto.backends.thermo.goodvibes import GoodVibesEngine
from hfauto.backends.kinetics.cantera import CanteraEngine
from hfauto.backends.thermo.arkane import ArkaneEngine
from hfauto.core.config import load_yaml
from hfauto.workflow.runner import run_pipeline
from hfauto.reporting.html_report import latest_manifest_path
from hfauto.core.io import read_manifest


def test_goodvibes_parser_accepts_qh_g_columns(tmp_path):
    csv_path = tmp_path / "Goodvibes_test.csv"
    csv_path.write_text("Structure,qh-G(T),H(T)\nfoo.out,-100.123456,-100.100000\n", encoding="utf-8")
    engine = GoodVibesEngine()
    table = engine.parse_goodvibes_csv(csv_path)
    assert "foo.out" in table.rows
    assert table.rows["foo.out"]["goodvibes_gibbs_hartree"] == -100.123456


def test_cantera_validation_is_structured_when_not_requested(tmp_path):
    engine = CanteraEngine(allow_external=False)
    path = tmp_path / "mechanism.yaml"
    path.write_text("phases: []\n", encoding="utf-8")
    result = engine.validate_mechanism(path)
    assert result["cantera_validation_status"] == "not_requested"
    assert result["production_cantera_ready"] is False


def test_arkane_export_summary_is_written(tmp_path):
    engine = ArkaneEngine(run_external=False)
    recs = [{"reaction_id": "rxn1", "T_K": 298.15, "delta_G_act_kcal_mol": 5.0, "quality_tier": "Q4"}]
    meta = engine.export_and_maybe_run(recs, tmp_path)
    assert Path(meta["skeleton_path"]).exists()
    assert meta["arkane_status"] == "not_requested"


def test_phase10_offline_pipeline_connector_audit(tmp_path):
    cfg = load_yaml("configs/pipelines/phase10_production_connectors_offline.yaml")
    cfg["run_root"] = str(tmp_path / "runs")
    run_dir = run_pipeline(cfg, "phase10_test")
    manifest = read_manifest(latest_manifest_path(run_dir))
    audit = manifest.find("production_connector_audit")
    assert audit is not None
    assert Path(audit.paths["csv"]).exists()
    assert manifest.find("reference_residuals") is not None

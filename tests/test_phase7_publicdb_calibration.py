from pathlib import Path

import pandas as pd

from hfauto.backends.db.pubchem import PubChemProvider
from hfauto.core.config import load_yaml
from hfauto.core.io import read_manifest
from hfauto.core.schemas.artifact import Artifact
from hfauto.reporting.html_report import latest_manifest_path
from hfauto.workflow.runner import run_pipeline


def test_pubchem_fixture_provider_enriches_without_network():
    provider = PubChemProvider(fixture_path="examples/db_fixtures/public_db_fixtures.json", allow_network=False)
    mol = Artifact(
        artifact_id="mol_ammonia",
        artifact_type="molecule",
        data={"mol_id": "mol00000", "name": "ammonia", "canonical_smiles": "N", "inchikey": "QGZKDVFQNNGYKY-UHFFFAOYSA-N"},
    )
    data = provider.enrich(mol)
    assert data["identity"]["pubchem_cid"] == 222
    assert data["public_data"]["pubchem"]["from_fixture"] is True
    assert data["db_provider_status"]["pubchem"]["matched"] is True


def test_phase7_publicdb_calibration_pipeline(tmp_path):
    cfg = load_yaml("configs/pipelines/phase7_publicdb_calibration_offline.yaml")
    cfg["run_root"] = str(tmp_path / "runs")
    run_dir = run_pipeline(cfg, "phase7_test")
    manifest = read_manifest(latest_manifest_path(run_dir))
    assert manifest.stage == "calibrate"
    assert any(a.artifact_type == "calibration" for a in manifest.artifacts)
    assert any(a.artifact_id == "method_validation_metrics" for a in manifest.artifacts)
    cov = next(a for a in manifest.artifacts if a.artifact_id == "public_data_coverage")
    df_cov = pd.read_csv(cov.paths["csv"])
    assert "public_db_support_score" in df_cov.columns
    assert df_cov["public_db_support_score"].max() > 0
    metrics = next(a for a in manifest.artifacts if a.artifact_id == "method_validation_metrics")
    df_metrics = pd.read_csv(metrics.paths["csv"])
    assert "basicity_rank_spearman_proxy" in set(df_metrics["metric"])
    report = next(a for a in manifest.artifacts if a.artifact_id == "method_validation_report")
    assert Path(report.paths["html"]).exists()

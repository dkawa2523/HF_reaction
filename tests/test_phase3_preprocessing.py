from pathlib import Path

from hfauto.backends.qm.xtb import parse_xtb_energy, parse_xtb_convergence
from hfauto.core.config import load_yaml
from hfauto.core.io import read_manifest
from hfauto.reporting.html_report import latest_manifest_path
from hfauto.workflow.runner import run_pipeline


def test_xtb_parsers():
    text = """
      | TOTAL ENERGY             -12.345678 Eh
      geometry optimization converged
      normal termination of xtb
    """
    assert parse_xtb_energy(text) == -12.345678
    qc = parse_xtb_convergence(text)
    assert qc["geometry_converged"] is True
    assert qc["scf_converged"] is True


def test_phase3_xtb_crest_fallback_pipeline(tmp_path):
    root = Path(__file__).parents[1]
    cfg = load_yaml(root / "configs" / "pipelines" / "phase3_xtb_crest.yaml")
    cfg["run_root"] = str(tmp_path / "runs")
    cfg["input"]["sdf"] = str(root / "examples" / "candidates.sdf")
    # Make the test deterministic even if crest/xtb happen to be installed.
    for stage in cfg["stages"]:
        if stage["name"] == "conformers":
            stage["allow_subprocess"] = False
        if stage["name"] == "preopt":
            stage["settings"]["allow_subprocess"] = False
    run_dir = run_pipeline(cfg, "phase3_fallback_test")
    manifest = read_manifest(latest_manifest_path(run_dir))
    assert manifest.stage == "rank"
    assert any(a.artifact_type == "conformer" and a.qc.get("fallback_reason") for a in manifest.artifacts)
    assert any(a.artifact_type == "species" and a.qc.get("preoptimized") for a in manifest.artifacts)
    assert any(a.artifact_type == "calculation" and a.qc.get("fallback_dummy") for a in manifest.artifacts)
    # The latest species revision should point to a preopt-updated xyz path.
    species = next(a for a in manifest.latest_artifacts("species") if a.data.get("state") == "reactant_complex")
    assert "preopt" in species.data

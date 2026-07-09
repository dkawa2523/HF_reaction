from __future__ import annotations

import math
from pathlib import Path

from hfauto.core.config import load_yaml
from hfauto.core.io import read_manifest
from hfauto.core.thermochemistry import eyring_rate, wigner_tunneling_factor
from hfauto.workflow.runner import run_pipeline


def test_wigner_factor_and_eyring_rate_are_positive():
    assert wigner_tunneling_factor(-1000.0, 423.15) > 1.0
    assert math.isclose(wigner_tunneling_factor(None, 423.15), 1.0)
    assert eyring_rate(5.0, 423.15) > 0.0


def test_phase6_offline_pipeline_exports_thermo_kinetics(tmp_path: Path):
    cfg = load_yaml(Path("configs/pipelines/phase6_thermo_kinetics_offline.yaml"))
    cfg["run_root"] = str(tmp_path / "runs")
    run_dir = run_pipeline(cfg, "phase6_test")
    manifest = read_manifest(Path((Path(run_dir) / "manifest.path").read_text()))
    thermo = list(manifest.iter_artifacts("thermo"))
    kinetics = list(manifest.iter_artifacts("kinetics"))
    assert thermo
    assert kinetics
    first = thermo[0].data
    assert "delta_G_assoc_standard_kcal_mol" in first
    assert "delta_G_assoc_pressure_corrected_kcal_mol" in first
    assert "K_assoc_standard" in first
    assert "quasi_rrho_applied" in first
    assert (Path(run_dir) / "12_kinetics" / "cantera_mechanism.yaml").exists()
    assert (Path(run_dir) / "12_kinetics" / "arkane" / "arkane_tst_input.py").exists()

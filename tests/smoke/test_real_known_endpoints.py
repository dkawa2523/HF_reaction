"""A real known_endpoints run of HCN -> HNC (WSL, HFAUTO_REAL=1; about a minute of QM) meets its
validation case (validation/cases.yaml ``hcn``): the outcome, the end-fragment equation and
dG_eff through reporting.validation.compare."""

import os
from pathlib import Path

import pytest

from hfauto.chemistry.xyz import read_xyz
from hfauto.pipeline.config import load
from hfauto.pipeline.runner import run_pipeline
from hfauto.reporting.validation import compare, load_cases

pytestmark = pytest.mark.real
REPO = Path(__file__).resolve().parents[2]


def test_hcn_known_endpoints_meets_its_case(real_site, tmp_path, monkeypatch):
    monkeypatch.chdir(REPO)  # methods and the site path resolve as for the CLI
    case = load_cases(REPO / "validation" / "cases.yaml")["hcn"]
    (step,) = case.chain
    resolved = load(REPO / "configs" / "pipelines" / f"{step.pipeline}.yaml",
                    REPO / "configs" / "systems" / f"{step.system}.yaml",
                    Path(os.environ["HFAUTO_SITE"]))
    layout = run_pipeline(resolved, tmp_path / "hcn")
    assert {s.status for s in layout.read_state()} == {"done"}
    view = layout.view()
    assert compare(case, view, lambda g: read_xyz(layout.run_dir / g.file.path)) == []

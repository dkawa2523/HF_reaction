"""Stage integration fixtures: the pipeline Runtime with fake engines (design §10.1)."""

from dataclasses import replace

import pytest

from hfauto.core.method import EngineSite
from hfauto.pipeline.config import PipelineConfig, ResolvedConfig, SiteConfig
from hfauto.pipeline.layout import RunLayout
from hfauto.pipeline.runner import build_runtime


@pytest.fixture
def fake_runtime(tmp_path):
    """``fake_runtime(system, {(capability, name): fake}, methods=None, stage_id="stage")``:
    a StageRuntime of the run ``tmp_path / "run"`` (``tmp_run``) built by build_runtime, whose
    ``create`` returns the fakes."""

    def make(system, engines, methods=None, stage_id="stage"):
        site = SiteConfig(site="fake", scratch_root=str(tmp_path / "scratch"), cores=4,
                          engines={n: EngineSite(version="0") for _, n in engines})
        resolved = ResolvedConfig(pipeline=PipelineConfig(pipeline_id="fake", stages=[]),
                                  system=system, site=site, methods=dict(methods or {}),
                                  code_version="test")
        layout = RunLayout(tmp_path / "run")
        layout.stage_dir(stage_id).mkdir(parents=True, exist_ok=True)
        runtime = replace(build_runtime(resolved, layout),
                          create=lambda capability, name, **_: engines[(capability, name)])
        return runtime.bind(stage_id, layout.stage_dir(stage_id))

    return make

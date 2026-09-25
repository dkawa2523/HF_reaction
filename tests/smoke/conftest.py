"""Real-engine smoke fixtures (marker ``real``): skipped unless HFAUTO_REAL=1 and HFAUTO_SITE."""

import os
from pathlib import Path

import pytest

from hfauto.pipeline.config import SiteConfig, load_site


@pytest.fixture
def real_site() -> SiteConfig:
    if os.environ.get("HFAUTO_REAL") != "1" or not os.environ.get("HFAUTO_SITE"):
        pytest.skip("real engines need HFAUTO_REAL=1 and HFAUTO_SITE=<site yaml>")
    return load_site(Path(os.environ["HFAUTO_SITE"]))


@pytest.fixture
def real_engine(real_site, tmp_path):
    """``real_engine(capability, name)``: a registry engine; jobs go to ``tmp_path/run/jobs``."""
    from hfauto.backends.engines import create
    from hfauto.execution.jobs import JobRunner
    from hfauto.execution.jobstore import JobStore

    jobs = JobRunner(JobStore(tmp_path / "run" / "jobs"), cores=real_site.cores)
    return lambda cap, name: create(cap, name, jobs=jobs, site=real_site.engines[name])

from __future__ import annotations

import platform
import sys

from hfauto.backends.protocols import Capability, Requirements
from hfauto.core.method import EngineSite
from hfauto.core.system import SystemConfig
from hfauto.pipeline.config import PipelineConfig, ResolvedConfig, SiteConfig
from hfauto.pipeline.preflight import check_site, engine_uses, preflight


def site(scratch: str, **engines: EngineSite) -> SiteConfig:
    return SiteConfig(site="s", scratch_root=scratch, cores=4, engines=engines)


def test_scratch_on_a_windows_drive_is_reported():
    ok = ResolvedConfig(pipeline=PipelineConfig(pipeline_id="p", stages=[]), methods={},
                        system=SystemConfig(system_id="s", species=[]),
                        site=site("/home/u/scratch"), code_version="t")
    assert preflight(ok, dry_run=True) == []
    mnt = site("/mnt/c/scratch", xtb=EngineSite(version="1", scratch_dir="/mnt/d/x"))
    found = preflight(ok.model_copy(update={"site": mnt}), dry_run=True)
    assert len(found) == 2 and all("ext4" in p for p in found)


def test_executables_version_pins_and_worker_modules(tmp_path):
    py = EngineSite(version=platform.python_version(), executables={"python": sys.executable})
    ok = Requirements(executables=("python",), python_modules=("json",),
                      version_command=("python", "--version"))
    bad = Requirements(executables=("no_such_exe_hfauto",), python_modules=("no_such_mod_hfauto",))
    assert check_site(site(str(tmp_path), py=py), {"py": ok}) == []
    pinned = site(str(tmp_path), py=py.model_copy(update={"version": "0.0.1"}), bad=py)
    found = check_site(pinned, {"py": ok, "bad": bad, "gone": ok})
    assert [p.split()[:2] for p in found] == [
        ["gone:", "not"], ["bad:", "executable"], ["bad:", "python"], ["py:", "version"]]
    assert check_site(pinned, {"bad": bad}, dry_run=True) == []


def test_engine_uses_pair_engines_with_methods():
    paths = {"id": "paths", "stage": "reaction-paths", "method": "pbe0",
             "engines": {"qm": "nwchem"}, "screen": {"method": "gfn2", "path": "pysis_neb"}}
    thermo = {"id": "thermo", "stage": "thermo", "engine": "goodvibes"}
    uses = engine_uses(PipelineConfig(pipeline_id="p", stages=[paths, thermo]))
    assert {(u.capability, u.name, u.method) for u in uses} == {(None, "goodvibes", None),
        (Capability.QM, "nwchem", "pbe0"), (Capability.PATH, "pysis_neb", "gfn2")}

from __future__ import annotations

import platform
import sys
from importlib.metadata import PackageNotFoundError

import pytest

from hfauto.backends.protocols import Capability, Requirements
from hfauto.core.method import EngineSite
from hfauto.core.system import SystemConfig
from hfauto.pipeline import preflight as pf
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
    crest = {"id": "conformers", "stage": "conformers", "engine": "crest"}
    uses = engine_uses(PipelineConfig(pipeline_id="p", stages=[paths, crest]))
    assert {(u.capability, u.name, u.method) for u in uses} == {(None, "crest", None),
        (Capability.QM, "nwchem", "pbe0"), (Capability.PATH, "pysis_neb", "gfn2")}


@pytest.mark.parametrize(("found", "problem"), [("4.3.0", False), ("4.2.1", True), (None, True)])
def test_a_thermo_stage_needs_goodvibes_4_3_0_in_this_interpreter(monkeypatch, found, problem):
    def version(name):
        if found is None:
            raise PackageNotFoundError(name)
        return found

    monkeypatch.setattr(pf, "version", version)
    thermo = PipelineConfig(pipeline_id="p", stages=[{"id": "thermo", "stage": "thermo"}])
    run = ResolvedConfig(pipeline=thermo, methods={}, system=SystemConfig(system_id="s", species=[]),
                         site=site("/home/u/scratch"), code_version="t")
    assert bool(preflight(run)) is problem and preflight(run, dry_run=True) == []
    assert pf.goodvibes_problems(thermo.model_copy(update={"stages": []})) == []

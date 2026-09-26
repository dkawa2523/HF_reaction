from __future__ import annotations

from types import SimpleNamespace

import pytest
import yaml

import hfauto.execution.jobs as jobs_module
from hfauto.backends import engines
from hfauto.backends.protocols import Capability
from hfauto.core.evidence import Failure, FailureKind, Geometry
from hfauto.core.manifest import Artifact
from hfauto.core.method import EngineSite
from hfauto.core.records import ArtifactType, ReportRecord
from hfauto.core.system import SpeciesInput, SystemConfig
from hfauto.pipeline.config import PipelineConfig, ResolvedConfig, SiteConfig
from hfauto.pipeline.layout import RunLayout
from hfauto.pipeline.runner import build_runtime, run_pipeline
from hfauto.stages import catalog
from hfauto.stages.spec import StageConfig, StageSpec

REPORT = ArtifactType.REPORT
CALLS: list[str] = []
E, C = {"id": "e", "stage": "emit"}, {"id": "c", "stage": "collect"}


class Ids(StageConfig):
    ids: tuple[str, ...] = ("a",)


class Emit:  # emits <stage>.<id>; Collect emits <stage>.n<number of artifacts in its view>
    spec = StageSpec("emit", Ids, consumes=(), produces=(REPORT,))

    def run(self, inputs, config, rt):
        CALLS.append(rt.stage_id)
        ids = config.ids if isinstance(config, Ids) else [f"n{len(inputs.artifacts)}"]
        rec = ReportRecord(rows=(), tables={})
        return [Artifact(artifact_id=f"{rt.stage_id}.{i}", type=REPORT, payload=rec) for i in ids]


class Collect(Emit):
    spec = StageSpec("collect", StageConfig, consumes=(REPORT,), produces=(REPORT,))


class Bad(Emit):
    spec = StageSpec("bad", Ids, consumes=(), produces=(ArtifactType.SPECIES,))


class Fail(Emit):  # leaves only failed artifacts
    spec = StageSpec("fail", Ids, consumes=(), produces=(REPORT,))

    def run(self, inputs, config, rt):
        why = Failure(kind=FailureKind.INPUT_INVALID, reason="reaction iso: atom index out of range")
        return [Artifact(artifact_id=f"{rt.stage_id}.{i}", type=REPORT, status="failed",
                         failure=why) for i in config.ids]


@pytest.fixture(autouse=True)
def dummy_stages():
    CALLS.clear()
    with catalog.override("emit", Emit), catalog.override("collect", Collect):  # noqa: SIM117
        with catalog.override("bad", Bad), catalog.override("fail", Fail):
            yield


def resolved(tmp_path, pipeline_id: str, *stages: dict) -> ResolvedConfig:
    site = SiteConfig(site="t", scratch_root=str(tmp_path / "scratch"), cores=2,
                      engines={"xtb": EngineSite(version="6.7.1")})
    system = SystemConfig(system_id="s", species=[SpeciesInput(id="h2", smiles="[H][H]")])
    return ResolvedConfig(pipeline=PipelineConfig(pipeline_id=pipeline_id, stages=list(stages)),
                          system=system, site=site, methods={}, code_version="test")


def test_resume_from_config_change_and_appending_pipeline(tmp_path):
    run, main = tmp_path / "run", resolved(tmp_path, "main", E, C)
    layout = run_pipeline(main, run)
    run_pipeline(main, run)  # resume: everything is fresh
    assert CALLS == ["e", "c"] and layout.state("c").n_ok == 1 and layout.state("c").input_sha
    assert [p.action for p in run_pipeline(main, run, start="e", dry_run=True)] == ["run", "run"]
    run_pipeline(main, run, start="e", stop="e")  # --from e --to e: c is left stale
    assert [(s.stage_id, s.status) for s in layout.read_state()] == [("c", "stale"), ("e", "done")]
    run_pipeline(main, run)
    run_pipeline(resolved(tmp_path, "main", {**E, "ids": ["a", "b"]}, C), run)  # emit changed
    assert CALLS == ["e", "c", "e", "c", "e", "c"]
    run_pipeline(resolved(tmp_path, "panel", {**C, "id": "pc"}), run)  # appends to the run
    assert [a.artifact_id for a in layout.view().artifacts] == ["e.a", "e.b", "c.n2", "pc.n3"]
    with pytest.raises(ValueError, match="already belongs"):
        run_pipeline(resolved(tmp_path, "other", E), run)
    records = yaml.safe_load(layout.resolved_config_path.read_text(encoding="utf-8"))
    assert list(records) == ["main", "panel"] and records["main"]["code_version"] == "test"
    assert records["panel"]["pipeline"]["stages"] == [{"id": "pc", "stage": "collect"}]


def test_consumes_and_produces_are_checked(tmp_path):
    with pytest.raises(ValueError, match="no input of type"):
        run_pipeline(resolved(tmp_path, "p", C), tmp_path / "r1")
    assert RunLayout(tmp_path / "r1").state("c").status == "failed" and CALLS == []
    with pytest.raises(ValueError, match="produced"):
        run_pipeline(resolved(tmp_path, "p", {"id": "b", "stage": "bad"}), tmp_path / "r2")


def test_missing_input_names_up_to_three_upstream_failures(tmp_path):
    fail = {"id": "f", "stage": "fail", "ids": ["a", "b", "c", "d"]}
    with pytest.raises(ValueError, match="upstream failures: f.a: reaction iso: atom") as exc:
        run_pipeline(resolved(tmp_path, "p", fail, C), tmp_path / "r")
    assert "f.c: reaction iso: atom index out of range" in str(exc.value)
    assert "f.d" not in str(exc.value)


class StubQM:  # structurally a QMEngine
    name = "xtb"
    requirements = supports = energy = optimize = frequencies = lambda self, *a, **k: True

    def __init__(self, jobs, site):
        self.jobs, self.site = jobs, site


def test_build_runtime_with_stub_runner_and_engine(tmp_path, monkeypatch):
    monkeypatch.setattr(jobs_module, "JobRunner",
                        lambda store, *, cores: SimpleNamespace(store=store, cores=cores))
    layout = RunLayout(tmp_path)
    with engines.override(Capability.QM, "xtb", StubQM):
        rt = build_runtime(resolved(tmp_path, "p"), layout).bind("s", layout.stage_dir("s"))
        qm = rt.engine(Capability.QM, "xtb")
        assert rt.engine(Capability.QM, "xtb") is qm and qm.site.version == "6.7.1"
        assert qm.jobs.cores == 2 and qm.jobs.store.root == layout.jobs_dir
    with pytest.raises(KeyError):
        rt.engine(Capability.QM, "nwchem")  # not configured in the site
    (xyz := tmp_path / "h.xyz").write_text("1\n\nH 0 0 0\n", encoding="utf-8")
    ref = rt.file_ref(xyz)
    assert ref.path == "h.xyz" and rt.stage_id == "s" and rt.run_id == tmp_path.name
    assert rt.load_xyz(Geometry(file=ref, fingerprint="f", symbols=("H",))).symbols == ["H"]
    assert rt.thread_map(lambda x: 2 * x, [1, 2, 3], threads_per_item=1) == [2, 4, 6]

from __future__ import annotations

from types import SimpleNamespace

import pytest
import yaml

from hfauto.backends import engines
from hfauto.backends.protocols import Capability
from hfauto.core.evidence import Failure, FailureKind, Geometry
from hfauto.core.manifest import Artifact, load_manifest
from hfauto.core.method import EngineSite
from hfauto.core.records import ArtifactType, ReportRecord
from hfauto.core.system import SpeciesInput, SystemConfig
from hfauto.execution.jobs import JobStats
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


class Contain(Emit):  # item "x" raises: contained to that item (a failed artifact)
    spec = StageSpec("contain", Ids, consumes=(), produces=(REPORT,))

    def run(self, inputs, config, rt):
        CALLS.append(rt.stage_id)

        def item(i):
            if i == "x":
                raise KeyError("Te")
            return ReportRecord(rows=(), tables={})

        got = {i: rt.contain(i, lambda i=i: item(i)) for i in config.ids}
        return [Artifact(artifact_id=f"{rt.stage_id}.{i}", type=REPORT, status="failed", failure=r)
                if isinstance(r, Failure) else
                Artifact(artifact_id=f"{rt.stage_id}.{i}", type=REPORT, payload=r)
                for i, r in got.items()]


class Fail(Emit):  # leaves only failed artifacts
    spec = StageSpec("fail", Ids, consumes=(), produces=(REPORT,))

    def run(self, inputs, config, rt):
        why = Failure(kind=FailureKind.INPUT_INVALID, reason="reaction iso: atom index out of range")
        return [Artifact(artifact_id=f"{rt.stage_id}.{i}", type=REPORT, status="failed",
                         failure=why) for i in config.ids]


@pytest.fixture(autouse=True)
def dummy_stages(monkeypatch):  # the runner finds these stage names through catalog.get
    CALLS.clear()
    stages = {"emit": Emit, "collect": Collect, "bad": Bad, "fail": Fail, "contain": Contain}
    get = catalog.get
    monkeypatch.setattr(catalog, "get", lambda name: stages.get(name) or get(name))


def resolved(tmp_path, pipeline_id: str, *stages: dict) -> ResolvedConfig:
    site = SiteConfig(site="t", scratch_root=str(tmp_path / "scratch"), cores=2,
                      engines={"xtb": EngineSite(version="6.7.1")})
    system = SystemConfig(system_id="s",
                          species=[SpeciesInput(id="h2", smiles="[H][H]", multiplicity=1)])
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
    layout = run_pipeline(resolved(tmp_path, "p", C, E), tmp_path / "r1")  # c has no input
    assert CALLS == ["e"] and load_manifest(layout.manifest_path("c")).artifacts == []
    assert (layout.state("c").status, layout.state("c").n_ok) == ("done", 0)
    with pytest.raises(ValueError, match="produced"):
        run_pipeline(resolved(tmp_path, "p", {"id": "b", "stage": "bad"}), tmp_path / "r2")
    assert RunLayout(tmp_path / "r2").state("b").status == "failed"


def test_an_item_that_raises_is_contained_counted_and_run_again(tmp_path, monkeypatch, caplog):
    """X7-2: the item, not the stage, fails (error:<type>); the stage is done with n_errors 1,
    which a rerun runs again. HFAUTO_STRICT=1 re-raises: the stage fails."""
    contain, run = {"id": "k", "stage": "contain", "ids": ["a", "x"]}, tmp_path / "r"
    with pytest.raises(KeyError):
        run_pipeline(resolved(tmp_path, "p", contain), run)
    assert RunLayout(run).state("k").status == "failed"
    monkeypatch.delenv("HFAUTO_STRICT")
    layout = run_pipeline(resolved(tmp_path, "p", contain), run)
    state = layout.state("k")
    assert (state.status, state.n_ok, state.n_failed, state.n_errors) == ("done", 1, 1, 1)
    failure = load_manifest(layout.manifest_path("k")).get("k.x").failure
    assert failure == Failure(kind=FailureKind.ERROR, reason="error:KeyError: 'Te'")
    assert "item 'x' raised" in caplog.text and state.unfinished()
    run_pipeline(resolved(tmp_path, "p", contain), run)  # not fresh: runs again
    assert CALLS == ["k", "k", "k"]


def test_a_stage_whose_job_failed_for_its_environment_runs_again(tmp_path):
    """X7-3: a missing executable is not stored, and its stage is not fresh: after a site fix
    a rerun runs it, with no --retry-failed."""
    run, main = tmp_path / "run", resolved(tmp_path, "main", E)
    layout = run_pipeline(main, run)
    assert [p.action for p in run_pipeline(main, run, dry_run=True)] == ["skip"]
    layout.update("e", jobs=JobStats(failures_by_kind={"executable_missing": 1}))
    assert [p.action for p in run_pipeline(main, run, dry_run=True)] == ["run"]


def test_missing_input_logs_up_to_three_upstream_failures(tmp_path, caplog):
    fail = {"id": "f", "stage": "fail", "ids": ["a", "b", "c", "d"]}
    run_pipeline(resolved(tmp_path, "p", fail, C), tmp_path / "r")
    assert "upstream failures: f.a: reaction iso: atom" in caplog.text
    assert "f.c: reaction iso: atom index out of range" in caplog.text
    assert "f.d" not in caplog.text


def test_build_runtime_with_stub_engine(tmp_path, monkeypatch):
    monkeypatch.setattr(engines, "create", lambda capability, name, *, jobs, site: (
        SimpleNamespace(jobs=jobs, site=site)))
    layout = RunLayout(tmp_path)
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
    assert rt.thread_map(lambda x: 2 * x, [1, 2, 3]) == [2, 4, 6]

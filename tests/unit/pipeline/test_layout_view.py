from hfauto.core.manifest import Artifact, Manifest, save_manifest
from hfauto.core.records import ArtifactType, ReportRecord
from hfauto.pipeline.layout import RunLayout


def finish(layout: RunLayout, stage_id: str, pipeline_id: str, *ids: str) -> None:
    layout.begin(stage_id, pipeline_id)
    payload = ReportRecord(rows=(), tables={})
    arts = [Artifact(artifact_id=i, type=ArtifactType.REPORT, payload=payload) for i in ids]
    save_manifest(Manifest(run_id="r", stage_id=stage_id, created_at="t", artifacts=arts),
                  layout.manifest_path(stage_id))
    layout.update(stage_id, status="done")


def test_view_is_the_union_of_earlier_done_stages(tmp_path):
    layout = RunLayout(tmp_path / "run")
    finish(layout, "s1", "main", "a", "b")
    finish(layout, "s2", "main", "b", "c")
    finish(layout, "panel", "panel", "d")  # another pipeline appends to the run
    view = {s: [a.artifact_id for a in layout.view(s).artifacts] for s in ("s1", "s2", "panel")}
    assert view == {"s1": [], "s2": ["a", "b"], "panel": ["a", "b", "c"]}
    assert [a.artifact_id for a in layout.view().artifacts] == ["a", "b", "c", "d"]


def test_run_state_transitions(tmp_path):
    layout = RunLayout(tmp_path / "run")
    assert layout.ensure("s0", "p").status == "pending"
    finish(layout, "s1", "p", "a")
    finish(layout, "s2", "p", "b")
    layout.begin("s1", "p")  # a done stage re-runs: moved last, stages after it go stale
    assert [(s.stage_id, s.status) for s in layout.read_state()] == [
        ("s0", "pending"), ("s2", "stale"), ("s1", "running")]
    layout.update("s1", status="done")
    layout.mark_stale("s0")  # --from a stage that never ran: nothing derives from it
    layout.mark_stale("s2")  # --from s2: s2 and everything executed after it
    assert [s.status for s in layout.read_state()] == ["pending", "stale", "stale"]

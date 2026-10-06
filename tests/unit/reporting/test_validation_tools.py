"""The validation harness must distinguish evidence, execution and unfinished evaluation."""

import gzip
import importlib.util
from pathlib import Path

import pytest

from hfauto.chemistry.xyz import read_xyz
from hfauto.pipeline.layout import RunLayout, StageState
from hfauto.reporting.validation import Case, Step, load_cases, load_references

REPO = Path(__file__).resolve().parents[3]
VALIDATION = REPO / "validation"


def _tool(name):
    spec = importlib.util.spec_from_file_location(name, VALIDATION / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_a_matched_regression_with_a_reference_deviation_does_not_pass(tmp_path, monkeypatch,
                                                                    capsys):
    check = _tool("check")
    source = REPO / "tests/golden/data/manifests/h_c2h4"
    run = tmp_path / "h_c2h4"
    layout = RunLayout(run)
    layout.write_state([StageState(stage_id=s, pipeline_id="known_endpoints", status="done")
                        for s in ("dft", "paths", "thermo")])
    for stage in ("dft", "paths", "thermo"):
        path = source / stage / "manifest.json"
        data = path.read_bytes() if path.is_file() else gzip.decompress(
            path.with_suffix(".json.gz").read_bytes())
        layout.stage_dir(stage).mkdir()
        layout.manifest_path(stage).write_bytes(data)
    monkeypatch.setattr(check, "_execution_problems", lambda *_: [])
    monkeypatch.setattr(check, "read_xyz", lambda path: read_xyz(source / path.relative_to(run)))
    case = load_cases(VALIDATION / "cases.yaml")["h_c2h4"]
    references = load_references(VALIDATION / "bh76/subset.yaml")
    assert check.check(run, run.name, case, references) is False
    output = capsys.readouterr().out
    assert output.startswith("DEVIATION h_c2h4") and "case expectations matched" in output
    assert "BH76 row 33" in output and "PASS" not in output


@pytest.mark.parametrize("changes", [None, {"status": "incomplete"}, {"status": "running"},
                                     {"status": "done", "n_errors": 1}])
def test_matching_records_cannot_hide_an_unfinished_required_stage(tmp_path, monkeypatch,
                                                                 capsys, changes):
    check = _tool("check")
    pipeline = tmp_path / "pipeline.yaml"
    pipeline.write_text("stages:\n  - {id: screen}\n  - {id: report}\n")
    case = Case(chain=(Step(pipeline=pipeline.as_posix(), system="unused"),))
    layout = RunLayout(tmp_path / "run")
    states = [StageState(stage_id="screen", pipeline_id="test", status="done")]
    if changes is not None:
        states.append(StageState(stage_id="report", pipeline_id="test", **changes))
    layout.write_state(states)
    monkeypatch.setattr(check, "compare", lambda *_: pytest.fail("unfinished records compared"))
    assert check.check(layout.run_dir, "run", case, {}) is False
    assert capsys.readouterr().out.startswith("INCOMPLETE run\n    report:")


def test_intentional_stage_bounds_and_superseded_evidence_are_respected(tmp_path):
    check = _tool("check")
    pipeline = tmp_path / "pipeline.yaml"
    pipeline.write_text("stages:\n  - {id: structures}\n  - {id: screen}\n  - {id: report}\n")
    case = Case(chain=(Step(pipeline=pipeline.as_posix(), system="unused",
                           args=("--from", "screen", "--to", "screen")),))
    layout = RunLayout(tmp_path / "run")
    layout.write_state([StageState(stage_id="screen", pipeline_id="test", status="done")])
    assert check._execution_problems(layout, case) == []
    archived = Case(chain=(Step(pipeline="removed", system="unused"),), superseded="replaced")
    assert check._execution_problems(layout, archived) == []


def test_complete_fresh_validation_requires_the_second_run(tmp_path, monkeypatch, capsys):
    check = _tool("check")
    case = Case(chain=(), twice=True)
    run = tmp_path / "boundary"
    run.mkdir()
    (run / "run_state.json").write_text("[]")
    monkeypatch.setattr(check, "load_cases", lambda *_: {"boundary": case})
    monkeypatch.setattr(check, "load_references", lambda *_: {})
    monkeypatch.setattr(check, "check", lambda *_: True)
    assert check.main([str(tmp_path), "boundary"]) == 0
    assert capsys.readouterr().out == "NOT_CHECKED boundary.2 (second fresh run)\n"
    assert check.main([str(tmp_path), "boundary", "--require-complete"]) == 1
    assert capsys.readouterr().out == "MISSING boundary.2 (second fresh run)\n"
    second = tmp_path / "boundary.2"
    second.mkdir()
    (second / "run_state.json").write_text("[]")
    assert check.main([str(tmp_path), "boundary", "--require-complete"]) == 0


def test_ran_records_execution_without_claiming_validation_and_detects_stopped_chains(monkeypatch):
    replay = _tool("replay")
    case = Case(chain=(Step(pipeline="one", system="unused"),
                       Step(pipeline="two", system="unused")))
    assert replay.chain_problems(case, Path("run"), [0, 1]) == []  # chemistry rc 1 is expected
    monkeypatch.setattr(replay, "run_chain", lambda *_: [124])
    result = replay.fresh(case, [Path("run")], Path("unused.log"), 1)
    assert result["verdict"] == "RAN"
    assert result["problems"] == ["run: ran 1 of 2 pipelines", "run: pipeline 1 exited 124"]


def test_strict_fresh_execution_exits_nonzero_when_its_chain_stops(tmp_path, monkeypatch, capsys):
    replay = _tool("replay")
    monkeypatch.setattr(replay, "ROOT", tmp_path)
    monkeypatch.setattr(replay, "run_chain", lambda *_: [124])  # no process or QM is started
    assert replay.main(["out", "--runs", "hcn", "--strict"]) == 1
    assert "RAN is execution only" in capsys.readouterr().out

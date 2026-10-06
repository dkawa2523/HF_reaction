"""CLI (design §7.4): five commands, a dry run on the shipped configs and status."""

import signal
import sys
from pathlib import Path

import pytest
import typer
from typer.testing import CliRunner

from hfauto.cli.main import app, stop_on_signals
from hfauto.execution import process
from hfauto.execution.jobs import JobStats
from hfauto.pipeline.layout import RunLayout, StageState

cli = CliRunner()


def test_help_lists_exactly_the_five_commands():
    assert cli.invoke(app, ["--help"]).exit_code == 0
    assert sorted(typer.main.get_command(app).commands) == [  # type: ignore[attr-defined]
        "case", "doctor", "report", "run", "status"]


@pytest.mark.parametrize("extra, code", [([], 0), (["--from", "nope"], 2),
                                         (["--retry-failed", "timeout,nope"], 2)])
def test_run_dry_run_resolves_names_and_writes_nothing(tmp_path, monkeypatch, extra, code):
    monkeypatch.chdir(Path(__file__).resolve().parents[3])
    result = cli.invoke(app, ["run", "known_endpoints", "--system", "hcn", "--site", "wsl_local",
                              "--run-dir", str(tmp_path / "run"), "--dry-run", *extra])
    assert result.exit_code == code, result.output
    assert code or all(s in result.output for s in ("structures", "dft", "paths", "report"))
    assert not (tmp_path / "run").exists()


def test_status_exit_code_is_the_runs_completeness(tmp_path):
    """X7-2: 0 when every stage is done, failed artifacts and chemical job failures included
    (coverage only); 1 when a stage is stopped, an item raised or a job failed for its
    environment (a rerun resumes it); 2 when a stage failed. As for run."""
    dft = StageState(stage_id="dft", pipeline_id="p", status="done", n_ok=3, n_failed=2,
                     jobs=JobStats(hits=1, misses=3, failures_by_kind={"scf_not_converged": 2}))
    paths = {"stage_id": "paths", "pipeline_id": "p", "status": "done",
             "jobs": JobStats(hits=1, failures_by_kind={"nonzero_exit": 1})}
    for change, code in (({}, 0), ({"status": "incomplete"}, 1), ({"n_errors": 1}, 1),
                         ({"jobs": JobStats(failures_by_kind={"timeout": 1})}, 1),
                         ({"status": "failed"}, 2)):
        RunLayout(tmp_path).write_state([dft, StageState.model_validate(paths | change)])
        result = cli.invoke(app, ["status", str(tmp_path)])
        assert result.exit_code == code, (change, result.output)
    assert "job failures: nonzero_exit=1, scf_not_converged=2" in result.output
    assert "job reuse: 2/5 (40%)" in result.output
    assert cli.invoke(app, ["status", str(tmp_path / "missing")]).exit_code == 2


@pytest.mark.skipif(sys.platform == "win32", reason="POSIX signals")
def test_a_signal_stops_the_external_programs_and_the_run_as_ctrl_c(monkeypatch):
    """SIGTERM stops the run as Ctrl-C does (KeyboardInterrupt): the running stage is left
    incomplete and run exits 1 (tests/integration/test_resume_e2e)."""
    handlers, stopped = {}, []
    monkeypatch.setattr(signal, "signal", lambda signum, handler: handlers.update({signum: handler}))
    monkeypatch.setattr(process, "stop_all", lambda: stopped.append(True))
    stop_on_signals()
    assert handlers[signal.SIGINT] is handlers[signal.SIGTERM] is handlers[signal.SIGHUP]
    with pytest.raises(KeyboardInterrupt, match="SIGTERM"):
        handlers[signal.SIGTERM](signal.SIGTERM, None)
    assert stopped == [True]

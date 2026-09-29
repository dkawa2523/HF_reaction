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


def test_status_sums_failures_by_kind_and_job_reuse(tmp_path):
    RunLayout(tmp_path).write_state([
        StageState(stage_id="dft", pipeline_id="p", status="done", n_ok=3,
                   jobs=JobStats(hits=1, misses=3, failures_by_kind={"timeout": 2})),
        StageState(stage_id="paths", pipeline_id="p", status="failed",
                   jobs=JobStats(hits=1, failures_by_kind={"timeout": 1, "nonzero_exit": 1}))])
    result = cli.invoke(app, ["status", str(tmp_path)])
    assert result.exit_code == 0, result.output
    assert "job failures: nonzero_exit=1, timeout=3" in result.output
    assert "job reuse: 2/5 (40%)" in result.output
    assert cli.invoke(app, ["status", str(tmp_path / "missing")]).exit_code == 2


@pytest.mark.skipif(sys.platform == "win32", reason="POSIX signals")
def test_ctrl_c_stops_the_external_programs_and_exits_130(monkeypatch):
    handlers, stopped = {}, []
    monkeypatch.setattr(signal, "signal", lambda signum, handler: handlers.update({signum: handler}))
    monkeypatch.setattr(process, "stop_all", lambda: stopped.append(True))
    stop_on_signals()
    assert handlers[signal.SIGINT] is handlers[signal.SIGTERM] is handlers[signal.SIGHUP]
    with pytest.raises(SystemExit) as exit_:
        handlers[signal.SIGINT](signal.SIGINT, None)
    assert exit_.value.code == 130 and stopped == [True]

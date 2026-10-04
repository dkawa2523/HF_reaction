import json
import sys
import time
from pathlib import Path

import pytest

from hfauto.execution import process
from hfauto.execution.lock import pid_alive
from hfauto.execution.process import Command, resolve_executable, run_command

# A command that starts a grandchild (like an MPI rank), writes its pid and sleeps.
_TREE = (
    "import os, subprocess, sys, time\n"
    "p = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(60)'])\n"
    "open('pid.tmp', 'w').write(str(p.pid)); os.replace('pid.tmp', 'child.pid')\n"
    "{then}\n"
    "time.sleep(60)"
)


def _py(tmp_path, code, **kwargs):
    return Command(argv=(sys.executable, "-c", code), cwd=tmp_path, **kwargs)


def _gone(pid_file: Path) -> bool:
    child = int(pid_file.read_text())  # a grandchild: only a tree kill stops it
    end = time.monotonic() + 5
    while pid_alive(child) and time.monotonic() < end:
        time.sleep(0.05)
    return not pid_alive(child)


def test_success_nonzero_env_stdin_and_sidecar(tmp_path):
    code = "import sys; print('out'); print('err', file=sys.stderr)"
    res = run_command(_py(tmp_path, code), timeout_s=30)
    assert (res.returncode, res.timed_out) == (0, False)
    assert res.stdout.read_text().strip() == "out" and res.stderr.read_text().strip() == "err"
    sidecar = json.loads((tmp_path / "command_result.json").read_text())
    assert sidecar["returncode"] == 0 and sidecar["argv"] == [sys.executable, "-c", code]
    (tmp_path / "in.txt").write_text("7")
    code = "import os, sys; sys.exit(int(sys.stdin.read()) + int(os.environ['HFAUTO_T']))"
    cmd = _py(tmp_path, code, env={"HFAUTO_T": "2"}, stdin=tmp_path / "in.txt")
    assert run_command(cmd, timeout_s=30).returncode == 9


def test_timeout_kills_the_process_tree(tmp_path):
    res = run_command(_py(tmp_path, _TREE.format(then="")), timeout_s=1)
    assert (res.timed_out, res.returncode) == (True, None) and res.duration_s < 10
    sidecar = json.loads((tmp_path / "command_result.json").read_text())
    assert sidecar["timed_out"] and "stopped" not in sidecar
    assert _gone(tmp_path / "child.pid")


def test_a_command_started_while_stopping_dies_and_is_no_result(tmp_path, monkeypatch):
    monkeypatch.setattr(process, "_stopping", True)
    start = time.monotonic()
    with pytest.raises(SystemExit):
        run_command(_py(tmp_path, "import time; time.sleep(60)"), timeout_s=30)
    assert time.monotonic() - start < 10 and not process._GROUPS


@pytest.mark.skipif(sys.platform == "win32", reason="POSIX signals and process groups")
def test_sigterm_kills_the_running_groups_and_stops_as_ctrl_c(tmp_path):
    # The command sends SIGTERM to hfauto (its parent), as `timeout` or a scheduler would; the
    # handler raises KeyboardInterrupt, which ``hfauto run`` turns into exit code 1.
    child = _TREE.format(then="import signal; os.kill(os.getppid(), signal.SIGTERM)")
    harness = (
        "import sys\nfrom pathlib import Path\n"
        "from hfauto.cli.main import stop_on_signals\n"
        "from hfauto.execution.process import Command, run_command\n"
        "stop_on_signals()\n"
        "try:\n"
        f"    run_command(Command(argv=(sys.executable, '-c', {child!r}), cwd=Path('job')),"
        " timeout_s=60)\n"
        "except KeyboardInterrupt as stop:\n"
        "    sys.exit(f'stopped by {stop}')\n"
    )
    root = str(Path(__file__).resolve().parents[3])
    res = run_command(_py(tmp_path, harness, env={"PYTHONPATH": root}), timeout_s=30)
    stderr = res.stderr.read_text()
    assert res.returncode == 1 and "stopped by SIGTERM" in stderr and res.duration_s < 20, stderr
    assert _gone(tmp_path / "job" / "child.pid")


def test_resolve_executable(tmp_path):
    assert resolve_executable("python", sys.executable) == sys.executable
    assert resolve_executable("python", str(tmp_path / "missing")) is None
    assert resolve_executable("hfauto-no-such-program") is None

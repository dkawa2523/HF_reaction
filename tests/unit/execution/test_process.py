import json
import sys
import time

from hfauto.execution.lock import pid_alive
from hfauto.execution.process import Command, resolve_executable, run_command


def _py(tmp_path, code, **kwargs):
    return Command(argv=(sys.executable, "-c", code), cwd=tmp_path, **kwargs)


def test_success_nonzero_env_stdin_and_sidecar(tmp_path):
    code = "import sys; print('out'); print('err', file=sys.stderr)"
    res = run_command(_py(tmp_path, code), timeout_s=30)
    assert (res.returncode, res.timed_out, res.stopped) == (0, False, None)
    assert res.stdout.read_text().strip() == "out" and res.stderr.read_text().strip() == "err"
    sidecar = json.loads((tmp_path / "command_result.json").read_text())
    assert sidecar["returncode"] == 0 and sidecar["argv"] == [sys.executable, "-c", code]
    (tmp_path / "in.txt").write_text("7")
    code = "import os, sys; sys.exit(int(sys.stdin.read()) + int(os.environ['HFAUTO_T']))"
    cmd = _py(tmp_path, code, env={"HFAUTO_T": "2"}, stdin=tmp_path / "in.txt")
    assert run_command(cmd, timeout_s=30).returncode == 9


def test_timeout(tmp_path):
    res = run_command(_py(tmp_path, "import time; time.sleep(30)"), timeout_s=0.5)
    assert (res.timed_out, res.returncode, res.stopped) == (True, None, None)
    assert res.duration_s < 8
    assert json.loads((tmp_path / "command_result.json").read_text())["timed_out"]


def test_monitor_stops_the_process_tree(tmp_path):
    code = (
        "import os, subprocess, sys, time\n"
        "p = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(60)'])\n"
        "open('pid.tmp', 'w').write(str(p.pid)); os.replace('pid.tmp', 'child.pid')\n"
        "time.sleep(60)"
    )

    def monitor(workdir):  # stops only if it is handed the command's working directory
        return "stagnated" if (workdir / "child.pid").exists() else None

    res = run_command(_py(tmp_path, code), timeout_s=30, monitor=monitor, poll_s=0.1)
    assert (res.stopped, res.timed_out, res.returncode) == ("stagnated", False, None)
    child = int((tmp_path / "child.pid").read_text())
    end = time.monotonic() + 5
    while pid_alive(child) and time.monotonic() < end:
        time.sleep(0.05)
    assert not pid_alive(child)


def test_resolve_executable(tmp_path):
    assert resolve_executable("python", sys.executable) == sys.executable
    assert resolve_executable("python", str(tmp_path / "missing")) is None
    assert resolve_executable("hfauto-no-such-program") is None

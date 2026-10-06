import json
import os
import sys
import time
from pathlib import Path

import pytest

from hfauto.execution.jobs import thread_map
from hfauto.execution.lock import SITE_LOCK_NAME, SiteLock
from hfauto.execution.process import Command, run_command


def echo(job, workdir):  # worker target: the job and where it ran
    return {"job": job, "cwd": str(Path.cwd()), "workdir": str(workdir)}


def test_site_lock_refuses_second_run_and_takes_over_dead_holders(tmp_path):
    scratch, lock = tmp_path / "scratch", tmp_path / "scratch" / SITE_LOCK_NAME
    with (
        SiteLock(scratch, tmp_path / "run_a"),
        pytest.raises(RuntimeError, match="run_a"),
        SiteLock(scratch, tmp_path / "run_b"),
    ):
        pass
    assert not lock.exists()
    code = "import os; print(os.getpid())"
    dead_pid = int(run_command(Command((sys.executable, "-c", code), tmp_path), timeout_s=30)
                   .stdout.read_text())
    # A dead pid, and this pid left by an earlier process (another token), are both stale.
    for stale in ({"pid": dead_pid}, {"pid": os.getpid(), "token": "earlier process"}):
        lock.write_text(json.dumps(stale | {"run_dir": "old", "time": "t"}))
        with SiteLock(scratch, tmp_path / "run_c"):
            holder = json.loads(lock.read_text())
            assert holder["pid"] == os.getpid() and holder["run_dir"].endswith("run_c")
    # U9-I14: a holder on another host is alive, whatever its pid is here
    lock.write_text(json.dumps({"pid": dead_pid, "host": "another-node", "run_dir": "old"}))
    with pytest.raises(RuntimeError, match="another-node"), SiteLock(scratch, tmp_path / "d"):
        pass


def test_worker_echo_runs_in_the_job_directory(tmp_path):
    job_dir = tmp_path / "attempt_00"
    job_dir.mkdir()
    (job_dir / "job.json").write_text(json.dumps({"n": 3}))
    argv = (sys.executable, "-m", "hfauto.execution.worker", "test_lock_worker:echo",
            str(job_dir / "job.json"))
    env = {"PYTHONPATH": str(Path(__file__).parent)}
    res = run_command(Command(argv, tmp_path, env), timeout_s=60)
    assert res.returncode == 0, res.stderr.read_text()
    out = json.loads((job_dir / "result.json").read_text())
    assert out["job"] == {"n": 3} and os.path.samefile(out["cwd"], job_dir)


def test_thread_map_keeps_input_order():
    def slow_identity(i):
        time.sleep(0.02 * (5 - i))
        return i

    assert thread_map(slow_identity, list(range(5)), workers=4) == list(range(5))
    assert thread_map(slow_identity, [3], workers=1) == [3]


def test_thread_map_cancels_the_items_not_started_when_one_raises():
    started = []

    def item(i):
        started.append(i)
        if i == 0:
            raise ValueError("item 0")
        time.sleep(0.2)
        return i

    with pytest.raises(ValueError, match="item 0"):
        thread_map(item, list(range(8)), workers=2)
    assert len(started) < 8

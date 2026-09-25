"""The only module that starts subprocesses (design §7.1).

``run_command`` streams stdout / stderr into files in the working directory, stops the
whole process tree on timeout or when a monitor asks for it, and writes
``command_result.json`` next to the output.
"""

from __future__ import annotations

import contextlib
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO

STDOUT_NAME = "stdout.txt"
STDERR_NAME = "stderr.txt"
RESULT_NAME = "command_result.json"
_TERMINATION_GRACE_S = 2.0

Monitor = Callable[[Path], str | None]


@dataclass(frozen=True)
class Command:
    argv: tuple[str, ...]
    cwd: Path
    env: Mapping[str, str] = field(default_factory=dict)  # overrides the inherited environment
    stdin: Path | None = None


@dataclass(frozen=True)
class CommandResult:
    returncode: int | None  # None when the process was killed (timeout or monitor stop)
    timed_out: bool
    stopped: str | None  # reason returned by the monitor
    duration_s: float
    stdout: Path
    stderr: Path


def resolve_executable(name: str, explicit: str | None = None) -> str | None:
    """The site's explicit executable (a path or a PATH name), else ``name`` on PATH.

    An explicit value that cannot be found gives None; there is no fallback to ``name``.
    """
    if explicit:
        path = Path(os.path.expandvars(explicit)).expanduser()
        if path.is_file() and os.access(path, os.X_OK):
            return str(path)
        return shutil.which(explicit)
    return shutil.which(name)


def run_command(
    cmd: Command,
    *,
    timeout_s: float,
    monitor: Monitor | None = None,
    poll_s: float = 10.0,
) -> CommandResult:
    """Run ``cmd`` in ``cmd.cwd`` (created if missing) and write ``command_result.json``.

    ``monitor(cwd)`` is called every ``poll_s`` seconds; a non-None reason stops the
    process tree and is recorded in ``stopped``. A missing executable raises
    FileNotFoundError (PermissionError when it is not executable).
    """
    cwd = Path(cmd.cwd)
    cwd.mkdir(parents=True, exist_ok=True)
    stdout, stderr = cwd / STDOUT_NAME, cwd / STDERR_NAME
    start = time.monotonic()
    with (
        stdout.open("wb") as out,
        stderr.open("wb") as err,
        _stdin(cmd.stdin) as inp,
    ):
        process = subprocess.Popen(
            list(cmd.argv),
            cwd=cwd,
            env={**os.environ, **cmd.env},
            stdin=inp,
            stdout=out,
            stderr=err,
            # A group / session of its own, so that the whole tree can be stopped.
            creationflags=subprocess.CREATE_NEW_PROCESS_GROUP if sys.platform == "win32" else 0,
            start_new_session=sys.platform != "win32",
        )
        stopped, timed_out = _supervise(process, cwd, start + timeout_s, monitor, poll_s)
    killed = timed_out or stopped is not None
    result = CommandResult(
        returncode=None if killed else process.returncode,
        timed_out=timed_out,
        stopped=stopped,
        duration_s=time.monotonic() - start,
        stdout=stdout,
        stderr=stderr,
    )
    _write_sidecar(cmd, result)
    return result


def _stdin(path: Path | None) -> contextlib.AbstractContextManager[IO[bytes] | int]:
    if path is None:
        return contextlib.nullcontext(subprocess.DEVNULL)
    return Path(path).open("rb")


def _supervise(
    process: subprocess.Popen[bytes],
    cwd: Path,
    end: float,
    monitor: Monitor | None,
    poll_s: float,
) -> tuple[str | None, bool]:
    """Wait for the process; returns (stop reason, timed out). Kills the tree on any exit path."""
    try:
        while True:
            remaining = end - time.monotonic()
            if remaining <= 0:
                _kill_tree(process)
                return None, True
            wait_s = remaining if monitor is None else min(poll_s, remaining)
            try:
                process.wait(timeout=wait_s)
                return None, False
            except subprocess.TimeoutExpired:
                pass
            reason = monitor(cwd) if monitor is not None else None
            if reason is not None:
                _kill_tree(process)
                return reason, False
    except BaseException:
        _kill_tree(process)
        raise


def _kill_tree(process: subprocess.Popen[bytes]) -> None:
    """Stop the process and its descendants (MPI ranks, workers); every wait is bounded."""
    if process.poll() is not None:
        return
    if sys.platform == "win32":
        with contextlib.suppress(OSError, subprocess.SubprocessError):
            subprocess.run(
                ["taskkill", "/PID", str(process.pid), "/T", "/F"],
                check=False,
                capture_output=True,
                timeout=_TERMINATION_GRACE_S,
            )
    else:
        # The command is a session leader, so this never reaches the caller's group.
        with contextlib.suppress(OSError):
            os.killpg(process.pid, signal.SIGKILL)
    if process.poll() is None:
        with contextlib.suppress(OSError):
            process.kill()
    with contextlib.suppress(subprocess.TimeoutExpired):
        process.wait(timeout=_TERMINATION_GRACE_S)


def _write_sidecar(cmd: Command, result: CommandResult) -> None:
    payload = {
        "argv": list(cmd.argv),
        "cwd": str(cmd.cwd),
        "returncode": result.returncode,
        "timed_out": result.timed_out,
        "stopped": result.stopped,
        "duration_s": result.duration_s,
        "stdout": str(result.stdout),
        "stderr": str(result.stderr),
    }
    path = Path(cmd.cwd) / RESULT_NAME
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    os.replace(tmp, path)

"""The only module that starts subprocesses (design §8).

``run_command`` streams stdout / stderr into files in the working directory, stops the
whole process tree on timeout, and writes ``command_result.json`` next to the output.
``stop_all`` kills every running command's group (the CLI's SIGTERM / SIGHUP handler).
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
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO

from hfauto.core.files import write_atomic

STDOUT_NAME = "stdout.txt"
STDERR_NAME = "stderr.txt"
RESULT_NAME = "command_result.json"
_TERMINATION_GRACE_S = 2.0
_GROUPS: set[int] = set()  # process groups of the running commands (the leader's pid)
_stopping = False


@dataclass(frozen=True)
class Command:
    argv: tuple[str, ...]
    cwd: Path
    env: Mapping[str, str] = field(default_factory=dict)  # overrides the inherited environment
    stdin: Path | None = None


@dataclass(frozen=True)
class CommandResult:
    returncode: int | None  # None when the process was killed on timeout
    timed_out: bool
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


def run_command(cmd: Command, *, timeout_s: float) -> CommandResult:
    """Run ``cmd`` in ``cmd.cwd`` (created if missing) and write ``command_result.json``.

    After ``timeout_s`` the whole process tree is killed and ``timed_out`` is set. A missing
    executable raises FileNotFoundError (PermissionError when it is not executable); after
    ``stop_all`` the call raises SystemExit.
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
            # A group / session of its own, so that the whole tree can be killed.
            creationflags=subprocess.CREATE_NEW_PROCESS_GROUP if sys.platform == "win32" else 0,
            start_new_session=sys.platform != "win32",
        )
        _GROUPS.add(process.pid)
        try:
            if _stopping:  # stop_all ran before the group was registered
                _kill_tree(process)
            timed_out = _supervise(process, timeout_s)
        finally:
            _GROUPS.discard(process.pid)
    if _stopping:  # a killed job is no result: no thread may record it
        raise SystemExit("stopped by a signal")
    result = CommandResult(
        returncode=None if timed_out else process.returncode,
        timed_out=timed_out,
        duration_s=time.monotonic() - start,
        stdout=stdout,
        stderr=stderr,
    )
    _write_sidecar(cmd, result)
    return result


def stop_all() -> None:
    """Kill the process group of every running command (POSIX); later commands die at once.

    Worker threads of ``thread_map`` then see their command end and raise SystemExit."""
    global _stopping
    _stopping = True
    if sys.platform != "win32":
        for group in list(_GROUPS):
            with contextlib.suppress(OSError):
                os.killpg(group, signal.SIGKILL)


def _stdin(path: Path | None) -> contextlib.AbstractContextManager[IO[bytes] | int]:
    if path is None:
        return contextlib.nullcontext(subprocess.DEVNULL)
    return Path(path).open("rb")


def _supervise(process: subprocess.Popen[bytes], timeout_s: float) -> bool:
    """Wait for the process; True when it timed out. Kills the tree on any abnormal exit."""
    try:
        process.wait(timeout=max(timeout_s, 0.0))
    except subprocess.TimeoutExpired:
        _kill_tree(process)
        return True
    except BaseException:
        _kill_tree(process)
        raise
    return False


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
        "duration_s": result.duration_s,
        "stdout": str(result.stdout),
        "stderr": str(result.stderr),
    }
    write_atomic(Path(cmd.cwd) / RESULT_NAME, json.dumps(payload, indent=1))

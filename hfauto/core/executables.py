"""Small subprocess/executable utilities used by simulation backends.

The utilities are deliberately minimal so backend behavior is easy to audit:
- an executable is either explicitly configured or found on PATH;
- every command writes stdout/stderr and a JSON sidecar;
- timeouts and non-zero exit codes are represented as structured results.
"""

from __future__ import annotations

import os
import shutil
import signal
import subprocess
import sys
import time
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from hfauto.core.io import ensure_dir, write_json_atomic

TIMEOUT_RETURNCODE = 124
_TERMINATION_GRACE_S = 2.0


@dataclass
class CommandResult:
    command: list[str]
    cwd: str
    returncode: int | None
    stdout_path: str
    stderr_path: str
    duration_s: float
    timed_out: bool = False
    executable: str | None = None

    @property
    def ok(self) -> bool:
        return self.returncode == 0 and not self.timed_out

    def to_dict(self) -> dict:
        return asdict(self)


def _popen_process_group_kwargs() -> dict[str, Any]:
    """Start each command in a killable group owned by this invocation."""

    if os.name == "posix":
        return {"start_new_session": True}
    if os.name == "nt":
        flag = getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0)
        return {"creationflags": flag} if flag else {}
    return {}


def _signal_posix_process_group(process: subprocess.Popen[str], sig: int) -> None:
    """Signal only the session/process group created for ``process``."""

    try:
        os.killpg(process.pid, sig)
    except ProcessLookupError:
        pass
    except OSError:
        if process.poll() is None:
            try:
                process.send_signal(sig)
            except OSError:
                pass


def _terminate_process_tree(
    process: subprocess.Popen[str],
) -> None:
    """Terminate the launcher and ranks, then reap the launcher.

    POSIX commands are session leaders, so signalling their process group also
    reaches MPI ranks and other descendants.  On Windows, ``taskkill /T`` is
    used for the exact PID created here; direct termination remains a fallback.
    Every wait is bounded to avoid replacing one timeout with another hang.
    """

    if os.name == "posix":
        # The process was started as a new session leader, so this cannot
        # target the caller's process group.
        _signal_posix_process_group(process, signal.SIGKILL)
    elif os.name == "nt" and process.poll() is None:
        try:
            subprocess.run(
                ["taskkill", "/PID", str(process.pid), "/T", "/F"],
                check=False,
                capture_output=True,
                timeout=_TERMINATION_GRACE_S,
            )
        except (OSError, subprocess.SubprocessError):
            pass
    if os.name != "posix" and process.poll() is None:
        try:
            process.kill()
        except OSError:
            pass

    try:
        process.wait(timeout=_TERMINATION_GRACE_S)
    except subprocess.TimeoutExpired:
        try:
            process.kill()
        except OSError:
            pass
        try:
            process.wait(timeout=_TERMINATION_GRACE_S)
        except (OSError, subprocess.TimeoutExpired):
            pass


def resolve_executable(name: str, explicit: str | None = None) -> str | None:
    """Return an executable path from explicit config or PATH, otherwise None."""
    env_name = f"HFAUTO_{name.upper().replace('-', '_')}_EXECUTABLE"
    explicit = explicit or os.environ.get(env_name)
    if explicit:
        p = Path(os.path.expandvars(explicit)).expanduser()
        if p.exists() and os.access(str(p), os.X_OK):
            return str(p)
        found = shutil.which(str(explicit))
        return found
    found = shutil.which(name)
    if found:
        return found
    scripts_dir = Path(sys.executable).parent
    candidates = [scripts_dir / name, scripts_dir / f"{name}.exe"]
    for candidate in candidates:
        if candidate.exists() and os.access(str(candidate), os.X_OK):
            return str(candidate)
    return None


def run_command(
    command: Sequence[str],
    cwd: str | Path,
    timeout_s: int = 3600,
    env: Mapping[str, str] | None = None,
    stdout_name: str = "stdout.txt",
    stderr_name: str = "stderr.txt",
) -> CommandResult:
    wd = ensure_dir(cwd)
    stdout_path = wd / stdout_name
    stderr_path = wd / stderr_name
    t0 = time.time()
    merged_env = os.environ.copy()
    if env:
        merged_env.update({str(k): str(v) for k, v in env.items()})
    command_list = list(command)
    # Stream long external calculations directly into their evidence files.
    # Keeping stdout/stderr in PIPE buffers until completion hides progress and
    # makes memory consumption grow with a multi-hour NEB/IRC calculation.
    with (
        stdout_path.open("w", encoding="utf-8") as stdout_stream,
        stderr_path.open("w", encoding="utf-8") as stderr_stream,
    ):
        process = subprocess.Popen(
            command_list,
            cwd=wd,
            text=True,
            encoding="utf-8",
            errors="replace",
            stdout=stdout_stream,
            stderr=stderr_stream,
            env=merged_env,
            **_popen_process_group_kwargs(),
        )
        try:
            process.wait(timeout=float(timeout_s))
            returncode = int(process.returncode)
            timed_out = False
        except subprocess.TimeoutExpired:
            _terminate_process_tree(process)
            returncode = TIMEOUT_RETURNCODE
            timed_out = True

    result = CommandResult(
        command=command_list,
        cwd=str(wd),
        returncode=returncode,
        stdout_path=str(stdout_path),
        stderr_path=str(stderr_path),
        duration_s=time.time() - t0,
        timed_out=timed_out,
        executable=command_list[0] if command_list else None,
    )
    write_json_atomic(wd / "command_result.json", result.to_dict())
    return result

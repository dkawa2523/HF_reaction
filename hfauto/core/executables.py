from __future__ import annotations

"""Small subprocess/executable utilities used by simulation backends.

The utilities are deliberately minimal so backend behavior is easy to audit:
- an executable is either explicitly configured or found on PATH;
- every command writes stdout/stderr and a JSON sidecar;
- timeouts and non-zero exit codes are represented as structured results.
"""

from dataclasses import dataclass, asdict
from pathlib import Path
import os
import shutil
import subprocess
import time
from typing import Mapping, Sequence

from hfauto.core.io import ensure_dir, write_json


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


def resolve_executable(name: str, explicit: str | None = None) -> str | None:
    """Return an executable path from explicit config or PATH, otherwise None."""
    if explicit:
        p = Path(explicit).expanduser()
        if p.exists() and os.access(str(p), os.X_OK):
            return str(p)
        found = shutil.which(str(explicit))
        return found
    return shutil.which(name)


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
    try:
        proc = subprocess.run(
            list(command),
            cwd=wd,
            text=True,
            capture_output=True,
            timeout=int(timeout_s),
            env=merged_env,
        )
        stdout_path.write_text(proc.stdout or "", encoding="utf-8")
        stderr_path.write_text(proc.stderr or "", encoding="utf-8")
        result = CommandResult(
            command=list(command),
            cwd=str(wd),
            returncode=int(proc.returncode),
            stdout_path=str(stdout_path),
            stderr_path=str(stderr_path),
            duration_s=time.time() - t0,
            timed_out=False,
            executable=list(command)[0] if command else None,
        )
    except subprocess.TimeoutExpired as exc:
        stdout_path.write_text(exc.stdout or "", encoding="utf-8")
        stderr_path.write_text(exc.stderr or "", encoding="utf-8")
        result = CommandResult(
            command=list(command),
            cwd=str(wd),
            returncode=None,
            stdout_path=str(stdout_path),
            stderr_path=str(stderr_path),
            duration_s=time.time() - t0,
            timed_out=True,
            executable=list(command)[0] if command else None,
        )
    write_json(wd / "command_result.json", result.to_dict())
    return result

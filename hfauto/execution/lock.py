"""Lock files created with O_EXCL; a lock whose holder pid is dead is taken over.

The same primitive serves ``SiteLock`` (one run per scratch root, design §7.1) and the
per-job ``.lock`` of the JobStore; plain files, so it works on Windows and POSIX alike.
"""

from __future__ import annotations

import json
import os
import sys
import time
import uuid
from collections.abc import Mapping
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Self

SITE_LOCK_NAME = ".hfauto_site.lock"
_UNREADABLE_STALE_S = 10.0  # a holder writes its JSON right after creating the file
# Tells this process apart from an earlier, dead one that had the same pid.
_TOKEN = uuid.uuid4().hex


def pid_alive(pid: object) -> bool:
    if not isinstance(pid, int) or pid <= 0:
        return False
    if pid == os.getpid():
        return True
    if sys.platform == "win32":
        return _windows_pid_alive(pid)
    try:
        os.kill(pid, 0)  # never on Windows: there os.kill terminates the process
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _windows_pid_alive(pid: int) -> bool:
    import ctypes
    from ctypes import wintypes

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.OpenProcess.restype = wintypes.HANDLE
    kernel32.OpenProcess.argtypes = (wintypes.DWORD, wintypes.BOOL, wintypes.DWORD)
    kernel32.GetExitCodeProcess.argtypes = (wintypes.HANDLE, ctypes.POINTER(wintypes.DWORD))
    kernel32.CloseHandle.argtypes = (wintypes.HANDLE,)
    handle = kernel32.OpenProcess(0x1000, False, pid)  # PROCESS_QUERY_LIMITED_INFORMATION
    if not handle:
        return ctypes.get_last_error() == 5  # ERROR_ACCESS_DENIED: exists, not ours
    try:
        code = wintypes.DWORD()
        ok = kernel32.GetExitCodeProcess(handle, ctypes.byref(code))
        return bool(ok) and code.value == 259  # STILL_ACTIVE
    finally:
        kernel32.CloseHandle(handle)


def holder_alive(holder: Mapping[str, Any]) -> bool:
    pid = holder.get("pid")
    if pid == os.getpid():
        return holder.get("token") == _TOKEN
    return pid_alive(pid)


def try_acquire(path: Path, holder: Mapping[str, Any]) -> dict[str, Any] | None:
    """Create ``path`` exclusively and write ``holder`` (plus this process's token) into it.

    Returns None on success, otherwise the current (live) holder. A lock whose holder is
    dead, or that stays unreadable, is removed and creation is retried.
    """
    while True:
        try:
            fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            current = _read_holder(path)
            if current is not None and holder_alive(current):
                return current
            if current is None and not _stale_unreadable(path):
                return {"pid": None, "unreadable": str(path)}
            release(path)
            continue
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump({**holder, "token": _TOKEN}, stream)
        return None


def release(path: Path) -> None:
    path.unlink(missing_ok=True)


def _read_holder(path: Path) -> dict[str, Any] | None:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def _stale_unreadable(path: Path) -> bool:
    try:
        return time.time() - path.stat().st_mtime > _UNREADABLE_STALE_S
    except OSError:
        return True  # vanished meanwhile


class SiteLock:
    """``<scratch_root>/.hfauto_site.lock``, held for a whole ``hfauto run``.

    Entering raises RuntimeError naming the holder when another live process (or this
    one) already holds the lock; a lock left by a dead pid is taken over.
    """

    def __init__(self, scratch_root: Path, run_dir: Path) -> None:
        self.path = Path(scratch_root) / SITE_LOCK_NAME
        self.run_dir = Path(run_dir)

    def __enter__(self) -> Self:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        holder = {
            "pid": os.getpid(),
            "run_dir": str(self.run_dir.resolve()),
            "time": datetime.now(UTC).isoformat(timespec="seconds"),
        }
        current = try_acquire(self.path, holder)
        if current is not None:
            raise RuntimeError(
                f"site lock {self.path} is held by pid {current.get('pid')} "
                f"(run {current.get('run_dir')}, since {current.get('time')})"
            )
        return self

    def __exit__(self, *exc_info: object) -> None:
        release(self.path)

"""The one atomic file write (standard library only: the worker process imports it)."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path


def write_atomic(path: str | Path, text: str) -> Path:
    """Write ``text`` through a temporary file in the same directory and ``os.replace``, so
    that a reader (or a resumed run after a kill) never sees a partial file."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=target.parent, prefix=f".{target.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp, target)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise
    return target

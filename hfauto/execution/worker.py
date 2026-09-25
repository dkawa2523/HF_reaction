"""Isolated in-process library runs.

Invoked as ``python -m hfauto.execution.worker <module:function> <job.json>``.

The working directory becomes the directory of ``job.json`` (the attempt directory),
``function(job, workdir)`` is called, and the returned dict is written to ``result.json``
there. SCINE, pysisyphus and GoodVibes are imported only inside such functions, so only
this child process ever loads them. Standard library only.
"""

from __future__ import annotations

import importlib
import json
import os
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

RESULT_NAME = "result.json"
_USAGE = "usage: python -m hfauto.execution.worker <module:function> <job.json>"

WorkerFunction = Callable[[dict[str, Any], Path], dict[str, Any]]


def _resolve(spec: str) -> WorkerFunction:
    module, _, name = spec.partition(":")
    if not module or not name:
        raise SystemExit(f"{_USAGE}\ninvalid target {spec!r}")
    return getattr(importlib.import_module(module), name)


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 2:
        print(_USAGE, file=sys.stderr)
        return 2
    function = _resolve(args[0])
    job_path = Path(args[1]).resolve()
    workdir = job_path.parent
    os.chdir(workdir)
    job = json.loads(job_path.read_text(encoding="utf-8"))
    result = function(job, workdir)
    if not isinstance(result, dict):
        raise TypeError(f"{args[0]} returned {type(result).__name__}, expected dict")
    tmp = workdir / (RESULT_NAME + ".tmp")
    tmp.write_text(json.dumps(result, indent=1), encoding="utf-8")
    os.replace(tmp, workdir / RESULT_NAME)  # a killed worker never leaves a partial result
    return 0


def _echo(job: dict[str, Any], workdir: Path) -> dict[str, Any]:
    """Test target: returns the job and where it ran."""
    return {"job": job, "cwd": str(Path.cwd()), "workdir": str(workdir)}


if __name__ == "__main__":
    sys.exit(main())

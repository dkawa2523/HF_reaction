"""Content-addressed job cache of one run: ``<run>/jobs/<k[:2]>/<key>/`` (design §7.1).

A job directory holds ``job.json`` (the key document), ``result.json`` and the
``attempt_NN/`` working directories. ``result.json`` is
``{"kind", "data": model_dump(mode="json"), "files": {relpath: sha256}}`` where ``files``
lists every FileRef inside the result (paths relative to the run directory).
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
import time
from collections.abc import Generator, Iterable, Iterator, Mapping
from contextlib import contextmanager
from pathlib import Path
from typing import TYPE_CHECKING, Any, TypeVar

from pydantic import BaseModel, ValidationError

from hfauto.core.evidence import Failure, FailureKind, FileRef
from hfauto.core.hashing import sha256_file
from hfauto.execution.lock import release, try_acquire

if TYPE_CHECKING:
    from hfauto.execution.jobs import Task

T = TypeVar("T", bound=BaseModel)

# Failures that a rerun of the same request cannot fix; remembered unless retried explicitly.
TERMINAL = frozenset(
    {
        FailureKind.INPUT_INVALID,
        FailureKind.METHOD_MISMATCH,
        FailureKind.EXECUTABLE_MISSING,
        FailureKind.INCOMPLETE_OUTPUT,
    }
)
_FAILURE_KIND = "failure"
_LOCK_POLL_S = 0.05


def key_document(task: Task) -> dict[str, Any]:
    """What identifies a job; ExecutionSpec and render inputs are not part of it."""
    return {
        "engine": task.engine,
        "version_pin": task.version_pin,
        "kind": task.kind,
        "key_payload": task.key_payload,
    }


def canonical_json(data: Any) -> str:
    return json.dumps(
        data, sort_keys=True, separators=(",", ":"), ensure_ascii=True, default=_jsonable
    )


def _jsonable(value: Any) -> Any:
    if isinstance(value, BaseModel):
        return value.model_dump(mode="json")
    if isinstance(value, Mapping):
        return dict(value)
    if isinstance(value, Path):
        return value.as_posix()
    raise TypeError(f"not JSON serializable in a job key: {type(value).__name__}")


def write_json_atomic(path: Path, data: Any) -> None:
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=path.name, suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            stream.write(json.dumps(data, indent=1, sort_keys=True, default=_jsonable))
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def file_refs(value: Any) -> Iterator[FileRef]:
    """Every FileRef reachable from a result model."""
    if isinstance(value, FileRef):
        yield value
    elif isinstance(value, BaseModel):
        for name in type(value).model_fields:
            yield from file_refs(getattr(value, name))
    elif isinstance(value, list | tuple):
        for item in value:
            yield from file_refs(item)
    elif isinstance(value, dict):
        for item in value.values():
            yield from file_refs(item)


class JobStore:
    """Cache under ``root`` (always ``<run>/jobs``); FileRef paths are relative to its parent."""

    def __init__(self, root: Path, *, retry_failed: Iterable[FailureKind | str] = ()) -> None:
        self.root = Path(root)
        self.run_dir = self.root.parent
        self.retry_failed = frozenset(FailureKind(kind) for kind in retry_failed)

    def key(self, task: Task) -> str:
        return hashlib.sha256(canonical_json(key_document(task)).encode("ascii")).hexdigest()

    def job_dir(self, key: str) -> Path:
        return self.root / key[:2] / key

    def attempt_dir(self, key: str, index: int) -> Path:
        return self.job_dir(key) / f"attempt_{index:02d}"

    @contextmanager
    def lock(self, key: str) -> Generator[None, None, None]:
        """Per-job exclusion via ``<job>/.lock``; waits while a live process holds it."""
        path = self.job_dir(key) / ".lock"
        path.parent.mkdir(parents=True, exist_ok=True)
        while try_acquire(path, {"pid": os.getpid()}) is not None:
            time.sleep(_LOCK_POLL_S)
        try:
            yield
        finally:
            release(path)

    def load(self, key: str, result_type: type[T]) -> T | Failure | None:
        """The stored result, a remembered terminal failure, or None (cache miss).

        Changed or missing files and data that no longer validate are cache misses.
        """
        path = self.job_dir(key) / "result.json"
        if not path.is_file():
            return None
        record = json.loads(path.read_text(encoding="utf-8"))
        try:
            if record["kind"] == _FAILURE_KIND:
                failure = Failure.model_validate(record["data"])
                return None if failure.kind in self.retry_failed else failure
            if not self._files_intact(record["files"]):
                return None
            return result_type.model_validate(record["data"])
        except ValidationError:
            return None

    def begin(self, key: str, task: Task) -> int:
        """Write ``job.json``; returns the index of the first new attempt directory."""
        job_dir = self.job_dir(key)
        job_dir.mkdir(parents=True, exist_ok=True)
        write_json_atomic(job_dir / "job.json", key_document(task))
        return len(list(job_dir.glob("attempt_*")))

    def save(self, key: str, result: BaseModel) -> None:
        """Store a result or a terminal failure; any other failure clears the old record."""
        path = self.job_dir(key) / "result.json"
        if isinstance(result, Failure) and result.kind not in TERMINAL:
            path.unlink(missing_ok=True)
            return
        kind = _FAILURE_KIND if isinstance(result, Failure) else getattr(result, "kind", None)
        files = {ref.path: ref.sha256 for ref in file_refs(result)}
        record = {"kind": kind, "data": result.model_dump(mode="json"), "files": files}
        path.parent.mkdir(parents=True, exist_ok=True)
        write_json_atomic(path, record)

    def file_ref(self, path: Path) -> FileRef:
        """FileRef of a file below the run directory (POSIX path relative to it)."""
        rel = Path(path).resolve().relative_to(self.run_dir.resolve())
        return FileRef(path=rel.as_posix(), sha256=sha256_file(path))

    def _files_intact(self, files: dict[str, str]) -> bool:
        for rel, sha in files.items():
            path = self.run_dir / rel
            if not path.is_file() or sha256_file(path) != sha:
                return False
        return True

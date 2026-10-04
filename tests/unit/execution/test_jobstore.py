import json
from types import MappingProxyType
from typing import Literal

from pydantic import BaseModel

from hfauto.core.evidence import Failure, FailureKind, FileRef, Geometry
from hfauto.core.method import ExecutionSpec
from hfauto.execution.jobs import Task
from hfauto.execution.jobstore import JobStore


class Out(BaseModel):
    kind: Literal["out"] = "out"
    files: tuple[FileRef, ...] = ()


def _task(**changes):
    fields = {"engine": "e", "version_pin": "1.0", "kind": "optimize",
              "key_payload": {"x": [1, 2.5]}, "execution": ExecutionSpec()}
    return Task(**(fields | changes))


def test_key_covers_request_and_pin_but_not_execution(tmp_path):
    store = JobStore(tmp_path / "jobs")
    key = store.key(_task())
    assert len(key) == 64 and key == store.key(_task(key_payload={"x": (1, 2.5)}))
    assert key == store.key(_task(key_payload=MappingProxyType({"x": [1, 2.5]})))
    assert key == store.key(_task(execution=ExecutionSpec(ranks=4, timeout_s=5), inputs={"a": 1}))
    for other in ({"version_pin": "1.1"}, {"kind": "string"}, {"key_payload": {"x": [1, 2.6]}}):
        assert store.key(_task(**other)) != key
    assert store.attempt_dir(key, 3) == tmp_path / "jobs" / key[:2] / key / "attempt_03"


def test_restores_result_type_and_checks_file_shas(tmp_path):
    store = JobStore(tmp_path / "jobs")
    key = store.key(_task())
    assert store.begin(key, _task()) == 0
    out_file = store.attempt_dir(key, 0) / "out.txt"
    out_file.parent.mkdir(parents=True)
    out_file.write_text("42")
    ref = store.file_ref(out_file)
    assert ref.path == f"jobs/{key[:2]}/{key}/attempt_00/out.txt"
    store.save(key, Out(files=(ref,)))
    record = json.loads((store.job_dir(key) / "result.json").read_text())
    assert record["kind"] == "out" and record["files"] == {ref.path: ref.sha256}
    assert json.loads((store.job_dir(key) / "job.json").read_text())["version_pin"] == "1.0"
    assert store.load(key, Out) == Out(files=(ref,)) and store.begin(key, _task()) == 1
    assert store.load(key, Failure) is None  # data that does not validate is a miss
    out_file.write_text("43")
    assert store.load(key, Out) is None


def test_failures_are_remembered_unless_retried(tmp_path):
    store = JobStore(tmp_path / "jobs")
    key = store.key(_task())
    for kind in set(FailureKind) - {FailureKind.TIMEOUT}:
        store.save(key, Failure(kind=kind, reason="r"))
        assert store.load(key, Out) == Failure(kind=kind, reason="r")
    frame = store.attempt_dir(key, 2) / "last.xyz"  # a saddle search stopped at maxiter
    frame.parent.mkdir(parents=True)
    frame.write_text("1\n\nH 0 0 0\n")
    final = Geometry(file=store.file_ref(frame), fingerprint="f", symbols=("H",))
    maxiter = Failure(kind=FailureKind.GEOMETRY_MAXITER, reason="maxiter", final=final)
    store.save(key, maxiter)
    assert store.load(key, Out) == maxiter
    assert JobStore(tmp_path / "jobs", retry_failed=["timeout"]).load(key, Out) == maxiter
    assert JobStore(tmp_path / "jobs", retry_failed=["geometry_maxiter"]).load(key, Out) is None
    frame.unlink()  # a deleted attempt directory
    assert store.load(key, Out) is None


def test_timeouts_are_never_stored(tmp_path):
    store = JobStore(tmp_path / "jobs")
    key = store.key(_task())
    store.save(key, Failure(kind=FailureKind.INPUT_INVALID, reason="autoz"))
    store.save(key, Failure(kind=FailureKind.TIMEOUT, reason="timeout"))  # clears the old record
    assert store.load(key, Out) is None
    assert not (store.job_dir(key) / "result.json").exists()


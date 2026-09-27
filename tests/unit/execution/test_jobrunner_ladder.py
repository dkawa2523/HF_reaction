import sys
from dataclasses import replace
from typing import Literal

import pytest
from pydantic import BaseModel

from hfauto.core.evidence import Failure, FailureKind, FileRef
from hfauto.core.method import Deadline, ExecutionSpec
from hfauto.execution.jobs import JobRunner, Task, thread_map
from hfauto.execution.jobstore import JobStore
from hfauto.execution.process import Command

OK = "open('out.txt', 'w').write('done')"
SLOW = "import time; time.sleep(30)"
BAD = "import sys; sys.exit(3)"


class Out(BaseModel):
    kind: Literal["out"] = "out"
    value: str
    file: FileRef
    job_key: str = ""


class ScriptAdapter:
    """Runs inputs['codes'][0] after dropping a unique marker file; continuation drops that code."""

    result_type = Out

    def __init__(self, store, marker):
        self.store, self.marker = store, marker

    def prepare(self, task, workdir):
        mark = f"import uuid; open({str(self.marker)!r} + '/' + uuid.uuid4().hex, 'w').close()\n"
        code = mark + task.inputs["codes"][0]
        return Command(argv=(task.inputs.get("exe", sys.executable), "-c", code), cwd=workdir)

    def parse(self, task, workdir, result):
        if result.timed_out:
            return Failure(kind=FailureKind.TIMEOUT, reason="timeout")
        if result.returncode != 0:
            kind = FailureKind.INPUT_INVALID if result.returncode == 3 else FailureKind.NONZERO_EXIT
            return Failure(kind=kind, reason=f"rc {result.returncode}")
        out = workdir / "out.txt"
        return Out(value=out.read_text(), file=self.store.file_ref(out))

    def continuation(self, task, workdir, failure):
        codes = task.inputs["codes"][1:]
        return replace(task, inputs={"codes": codes}, execution=ExecutionSpec()) if codes else None


def _setup(tmp_path, cores=1, **store_options):
    store = JobStore(tmp_path / "run" / "jobs", **store_options)
    (tmp_path / "marker").mkdir(exist_ok=True)
    return JobRunner(store, cores=cores), ScriptAdapter(store, tmp_path / "marker")


def _task(*codes, name="a", timeout_s=30.0, ranks=1):
    execution = ExecutionSpec(timeout_s=timeout_s, ranks=ranks)
    return Task("script", "1", "run", {"name": name}, execution, {"codes": codes})


def _runs(tmp_path):
    marker = tmp_path / "marker"
    return len(list(marker.iterdir())) if marker.exists() else 0


def test_timeout_continues_then_result_is_reused(tmp_path):
    runner, adapter = _setup(tmp_path)
    task = _task(SLOW, OK, timeout_s=0.5)
    out = runner.run(task, adapter)
    key = runner.store.key(task)
    assert isinstance(out, Out) and out.value == "done" and out.job_key == key
    assert runner.store.attempt_dir(key, 1).exists()
    runs = _runs(tmp_path)
    assert runner.run(task, adapter) == out
    fresh, _ = _setup(tmp_path)  # a new process restores through result_type
    assert fresh.run(task, adapter) == out and _runs(tmp_path) == runs
    assert (runner.stats().hits, runner.stats().misses, fresh.stats().misses) == (1, 1, 0)


def test_ladder_limits_and_terminal_failure_memory(tmp_path):
    runner, adapter = _setup(tmp_path)
    first = runner.run(_task(BAD, BAD, OK), adapter)  # INPUT_INVALID is continued only once
    assert first.kind is FailureKind.INPUT_INVALID and first.job_key and _runs(tmp_path) == 2
    assert runner.run(_task(BAD, BAD, OK), adapter) == first and _runs(tmp_path) == 2
    assert runner.run(_task(BAD, name="b"), adapter).kind is FailureKind.INPUT_INVALID
    assert _runs(tmp_path) == 3  # no continuation offered: a single attempt
    nonzero = _task("import sys; sys.exit(1)", OK, name="c")
    assert runner.run(nonzero, adapter).kind is FailureKind.NONZERO_EXIT  # never continued
    assert runner.run(nonzero, adapter).kind is FailureKind.NONZERO_EXIT  # nor remembered
    assert _runs(tmp_path) == 5
    assert runner.stats().failures_by_kind == {"input_invalid": 3, "nonzero_exit": 2}
    retry, _ = _setup(tmp_path, retry_failed=["input_invalid"])
    assert isinstance(retry.run(_task(BAD, BAD, OK), adapter), Failure) and _runs(tmp_path) == 7


def test_budget_and_missing_executable(tmp_path):
    runner, adapter = _setup(tmp_path)
    task = _task(OK)
    out = runner.run(task, adapter, deadline=Deadline.after(30))
    assert out.kind is FailureKind.BUDGET_EXHAUSTED and _runs(tmp_path) == 0
    assert isinstance(runner.run(task, adapter, deadline=Deadline.after(3600)), Out)
    missing = replace(_task(OK, name="m"), inputs={"codes": (OK,), "exe": "hfauto-no-such-exe"})
    missing = runner.run(missing, adapter)
    assert missing.kind is FailureKind.EXECUTABLE_MISSING
    assert (runner.store.job_dir(missing.job_key) / "result.json").exists()


def test_parse_exception_is_a_remembered_incomplete_output(tmp_path, monkeypatch):
    runner, adapter = _setup(tmp_path)

    def parse(task, workdir, result):
        raise ValueError("no energy line\nin the output")

    adapter.parse = parse
    task = _task(OK, OK, name="p")  # a continuation exists but is not offered
    with pytest.raises(ValueError, match="no energy line"):  # HFAUTO_STRICT=1 re-raises
        runner.run(task, adapter)
    monkeypatch.delenv("HFAUTO_STRICT")
    out = runner.run(task, adapter)
    assert out == Failure(kind=FailureKind.INCOMPLETE_OUTPUT, job_key=runner.store.key(task),
                          reason="parse:ValueError: no energy line")
    fresh, _ = _setup(tmp_path)  # replayed from the JobStore on disk
    assert fresh.run(task, adapter) == out and _runs(tmp_path) == 2
    retry, _ = _setup(tmp_path, retry_failed=["incomplete_output"])
    assert retry.run(task, adapter) == out and _runs(tmp_path) == 3


def _meet(runner, adapter, room, wait_s, ranks):
    """Two jobs each touch a file in ``room`` and wait for the other; returns what each saw."""
    code = (f"import pathlib, time; d = pathlib.Path({str(room)!r}); d.mkdir(exist_ok=True)\n"
            f"(d / pathlib.Path.cwd().parent.name).touch(); end = time.time() + {wait_s}\n"
            "while len(list(d.iterdir())) < 2 and time.time() < end: time.sleep(0.02)\n"
            "open('out.txt', 'w').write(str(len(list(d.iterdir()))))")
    tasks = [_task(code, name=f"{room.name}{i}", ranks=ranks) for i in range(2)]
    outs = thread_map(lambda t: runner.run(t, adapter), tasks, workers=2)
    assert [o.job_key for o in outs] == [runner.store.key(t) for t in tasks]
    return sorted(out.value for out in outs)


def test_core_semaphore_serializes_and_job_lock_dedupes(tmp_path):
    runner, adapter = _setup(tmp_path, cores=2)
    assert _meet(runner, adapter, tmp_path / "wide", 0.3, 2) == ["1", "2"]  # one at a time
    assert _meet(runner, adapter, tmp_path / "narrow", 10, 1) == ["2", "2"]  # concurrent
    same = thread_map(lambda t: runner.run(t, adapter), [_task(OK, name="s")] * 2, workers=2)
    assert same[0] == same[1] and _runs(tmp_path) == 5

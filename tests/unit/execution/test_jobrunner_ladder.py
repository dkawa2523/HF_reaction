import sys
from dataclasses import replace
from typing import Literal

import pytest
from pydantic import BaseModel

from hfauto.core.evidence import Failure, FailureKind, FileRef
from hfauto.core.method import ExecutionSpec
from hfauto.execution.jobs import JobRunner, JobStats, Task, thread_map
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
    """Runs inputs['codes'][0] after dropping a unique marker file; continuation drops that code.
    ``seen`` holds the execution of every attempt, ``continued`` (attempt, failure kind) of
    every continuation asked."""

    result_type = Out

    def __init__(self, store, marker):
        self.store, self.marker, self.seen, self.continued = store, marker, [], []

    def prepare(self, task, workdir):
        self.seen.append(task.execution)
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
        self.continued.append((workdir.name, failure.kind))
        codes = task.inputs["codes"][1:]
        return replace(task, inputs={"codes": codes}, execution=ExecutionSpec()) if codes else None


def _setup(tmp_path, cores=1, **store_options):
    store = JobStore(tmp_path / "run" / "jobs", **store_options)
    (tmp_path / "marker").mkdir(exist_ok=True)
    return JobRunner(store, cores=cores), ScriptAdapter(store, tmp_path / "marker")


def _task(*codes, name="a", timeout_s=30.0, ranks=1, threads=1):
    execution = ExecutionSpec(timeout_s=timeout_s, ranks=ranks, threads=threads)
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


def test_ladder_limits_and_failure_memory(tmp_path):
    runner, adapter = _setup(tmp_path)
    first = runner.run(_task(BAD, BAD, OK), adapter)  # INPUT_INVALID is continued only once
    assert first.kind is FailureKind.INPUT_INVALID and first.job_key and _runs(tmp_path) == 2
    assert runner.run(_task(BAD, BAD, OK), adapter) == first and _runs(tmp_path) == 2
    assert runner.run(_task(BAD, name="b"), adapter).kind is FailureKind.INPUT_INVALID
    assert _runs(tmp_path) == 3  # no continuation offered: a single attempt
    before, nonzero = runner.stats(), _task("import sys; sys.exit(1)", OK, name="c")
    assert runner.run(nonzero, adapter).kind is FailureKind.NONZERO_EXIT  # never continued
    assert runner.run(nonzero, adapter).kind is FailureKind.NONZERO_EXIT  # but remembered
    assert _runs(tmp_path) == 4
    assert runner.stats().failures_by_kind == {"input_invalid": 3, "nonzero_exit": 2}
    assert runner.stats().since(before) == JobStats(hits=1, misses=1,
                                                    failures_by_kind={"nonzero_exit": 2})
    retry, _ = _setup(tmp_path, retry_failed=["input_invalid"])
    assert isinstance(retry.run(_task(BAD, BAD, OK), adapter), Failure) and _runs(tmp_path) == 6


def test_environment_failures_are_not_stored_and_rerun_after_a_site_fix(tmp_path):
    """X7-3: a timeout, a missing executable and a job out of memory are never stored."""
    runner, adapter = _setup(tmp_path)
    before, slow = runner.stats(), _task(SLOW, name="slow", timeout_s=0.5)  # no continuation
    assert [runner.run(slow, adapter).kind for _ in range(2)] == [FailureKind.TIMEOUT] * 2
    assert runner.stats().since(before) == JobStats(misses=2, failures_by_kind={"timeout": 2})
    missing = replace(_task(OK, name="m"), inputs={"codes": (OK,), "exe": "hfauto-no-such-exe"})
    assert runner.run(missing, adapter).kind is FailureKind.EXECUTABLE_MISSING
    assert not runner.store.has_result(runner.store.key(missing))
    fixed = runner.run(replace(missing, inputs={"codes": (OK,)}), adapter)  # the same key
    assert isinstance(fixed, Out) and fixed.job_key == runner.store.key(missing)
    oom = "import sys; sys.stderr.write('ga_create failed: Cannot allocate memory'); sys.exit(1)"
    for code, kind in ((oom, FailureKind.OUT_OF_MEMORY),
                       ("import sys; sys.exit(1)", FailureKind.NONZERO_EXIT)):
        task = _task(code, name=code)
        assert runner.run(task, adapter).kind is kind
        assert runner.store.has_result(runner.store.key(task)) is (kind is FailureKind.NONZERO_EXIT)


def test_a_parse_exception_is_the_items_and_stores_nothing(tmp_path):
    """A parser that raises fails the item that ran the job (StageRuntime.contain), not the
    job: nothing is stored, so a fixed parser reads the job again."""
    runner, adapter = _setup(tmp_path)
    parse = adapter.parse

    def broken(task, workdir, result):
        raise ValueError("no energy line")

    adapter.parse = broken
    task = _task(OK, OK, name="p")  # a continuation exists but is not offered
    with pytest.raises(ValueError, match="no energy line"):
        runner.run(task, adapter)
    assert not runner.store.has_result(runner.store.key(task)) and not adapter.continued
    adapter.parse = parse
    assert isinstance(runner.run(task, adapter), Out)


def test_an_attempt_that_left_no_result_is_continued_once(tmp_path):
    """U9-P3: a stopped run left attempt_00 without a result; the next run continues it as
    after a timeout, in attempt_01; with no continuation it starts afresh."""
    runner, adapter = _setup(tmp_path)
    for name, codes in (("stopped", (BAD, OK)), ("one", (OK,))):
        task = _task(*codes, name=name)
        key = runner.store.key(task)
        runner.store.begin(key, task)
        runner.store.attempt_dir(key, 0).mkdir()  # killed before it wrote anything
        adapter.continued.clear()
        assert isinstance(runner.run(task, adapter), Out)  # BAD is never run
        assert adapter.continued == [("attempt_00", FailureKind.TIMEOUT)]
        assert (runner.store.attempt_dir(key, 1) / "out.txt").is_file()
    assert _runs(tmp_path) == 2


def _meet(runner, adapter, room, wait_s, **execution):
    """Two jobs each touch a file in ``room`` and wait for the other; returns what each saw."""
    code = (f"import pathlib, time; d = pathlib.Path({str(room)!r}); d.mkdir(exist_ok=True)\n"
            f"(d / pathlib.Path.cwd().parent.name).touch(); end = time.time() + {wait_s}\n"
            "while len(list(d.iterdir())) < 2 and time.time() < end: time.sleep(0.02)\n"
            "open('out.txt', 'w').write(str(len(list(d.iterdir()))))")
    tasks = [_task(code, name=f"{room.name}{i}", **execution) for i in range(2)]
    outs = thread_map(lambda t: runner.run(t, adapter), tasks, workers=2)
    assert [o.job_key for o in outs] == [runner.store.key(t) for t in tasks]
    return sorted(out.value for out in outs)


def test_core_semaphore_serializes_and_job_lock_dedupes(tmp_path):
    runner, adapter = _setup(tmp_path, cores=2)
    assert _meet(runner, adapter, tmp_path / "wide", 0.3, threads=2) == ["1", "2"]  # in turn
    assert _meet(runner, adapter, tmp_path / "narrow", 10) == ["2", "2"]  # concurrent
    assert _meet(runner, adapter, tmp_path / "ranks", 10, ranks=2) == ["2", "2"]  # 1 rank each
    same = thread_map(lambda t: runner.run(t, adapter), [_task(OK, name="s")] * 2, workers=2)
    assert same[0] == same[1] and _runs(tmp_path) == 7


def test_rank_share_applies_only_inside_thread_map_and_never_to_the_key(tmp_path):
    runner, adapter = _setup(tmp_path, cores=4)

    def run(task):
        return runner.run(task, adapter)

    three = thread_map(run, [_task(OK, name=f"t{i}", ranks=4) for i in range(3)], workers=4)
    two = thread_map(run, [_task(OK, name=f"d{i}", ranks=4) for i in range(2)], workers=4)
    thread_map(run, [_task(OK, name="one", ranks=4)], workers=4)
    run(_task(OK, name="serial", ranks=4))
    # ranks 4 // (items at once); the timeout keeps the job's core-seconds
    assert [(e.ranks, e.timeout_s) for e in adapter.seen] == [
        (1, 120.0)] * 3 + [(2, 60.0)] * 2 + [(4, 30.0)] * 2
    full = _task(OK, name="t0", ranks=4)
    assert runner.store.key(full) == runner.store.key(replace(full, execution=ExecutionSpec()))
    assert run(full) == three[0] and two[0].job_key and len(adapter.seen) == 7  # a hit
    thread_map(run, [_task(OK, name=f"w{i}", ranks=4) for i in range(5)], workers=4)
    assert sorted(e.ranks for e in adapter.seen[7:]) == [1, 1, 1, 1, 4]  # waves of 4 and 1

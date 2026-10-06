"""Tasks, adapters and the JobRunner: cache, attempt ladder, core semaphore and the rank share
of concurrent items (design §7.1)."""

from __future__ import annotations

import re
import threading
from collections import Counter
from collections.abc import Callable, Generator, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Protocol, TypeVar

from pydantic import BaseModel, ConfigDict

from hfauto.core.evidence import Failure, FailureKind
from hfauto.core.method import ExecutionSpec
from hfauto.execution.jobstore import JobStore
from hfauto.execution.process import Command, CommandResult, run_command

T = TypeVar("T", bound=BaseModel)
ItemT = TypeVar("ItemT")
ResultT = TypeVar("ResultT")

# The only automatic retries: failure kind -> continuations allowed (each needs a non-None
# adapter.continuation). Every other failure is returned at once.
LADDER: dict[FailureKind, int] = {
    FailureKind.TIMEOUT: 2,
    FailureKind.GEOMETRY_MAXITER: 2,
    FailureKind.INPUT_INVALID: 1,
    FailureKind.SCF_NOT_CONVERGED: 1,
}
# Set by thread_map in its workers: the MPI ranks one of several concurrent items may use.
_RANK_SHARE: ContextVar[int | None] = ContextVar("rank_share", default=None)
# A nonzero exit that is an allocation failure (OOM, MPI / GA / MA, malloc) in the output's tail.
_OUT_OF_MEMORY = re.compile(
    r"out of memory|cannot allocate memory|bad_alloc|MemoryError|insufficient memory"
    r"|not enough memory|(?:ga_create|ma_push_get|ma_alloc_get|armci_malloc)\b.{0,40}fail"
    r"|signal 9 \(Killed\)", re.IGNORECASE)
_TAIL_BYTES = 20_000


@dataclass(frozen=True)
class Task:
    engine: str
    version_pin: str  # EngineSite.version, known before running
    kind: str  # "optimize", "frequencies", "string", ...
    key_payload: Mapping[str, Any]  # method.signature(), fingerprints, parameters, input shas
    execution: ExecutionSpec
    inputs: Mapping[str, Any] = field(default_factory=dict)  # for prepare only; not in the key


class Adapter(Protocol[T]):
    result_type: type[T]  # pydantic model used to restore results from the JobStore

    def prepare(self, task: Task, workdir: Path) -> Command: ...

    def parse(self, task: Task, workdir: Path, result: CommandResult) -> T | Failure: ...

    def continuation(self, task: Task, workdir: Path, failure: Failure) -> Task | None: ...


class JobStats(BaseModel):
    """Cache hits, misses and failures by kind of a JobRunner, or of one stage (``since``)."""

    model_config = ConfigDict(frozen=True, extra="forbid")
    hits: int = 0
    misses: int = 0
    failures_by_kind: dict[str, int] = {}

    def since(self, before: JobStats) -> JobStats:
        """The counts added after ``before`` (the counters only grow)."""
        failures = {kind: n - before.failures_by_kind.get(kind, 0)
                    for kind, n in self.failures_by_kind.items()}
        return JobStats(hits=self.hits - before.hits, misses=self.misses - before.misses,
                        failures_by_kind={kind: n for kind, n in failures.items() if n})


class _CoreSemaphore:
    """Reserves cores so that the ranks x threads of running jobs never exceed the total."""

    def __init__(self, total: int) -> None:
        self._total = total
        self._free = total
        self._cond = threading.Condition()

    @contextmanager
    def reserve(self, n: int) -> Generator[None, None, None]:
        n = min(max(n, 1), self._total)
        with self._cond:
            self._cond.wait_for(lambda: self._free >= n)
            self._free -= n
        try:
            yield
        finally:
            with self._cond:
                self._free += n
                self._cond.notify_all()


class JobRunner:
    """Runs tasks through the JobStore: cached result, or attempts along LADDER."""

    def __init__(self, store: JobStore, *, cores: int) -> None:
        if cores < 1:
            raise ValueError(f"cores must be >= 1, got {cores}")
        self.store = store
        self.cores = cores
        self._semaphore = _CoreSemaphore(cores)
        self._lock = threading.Lock()
        self._hits = 0
        self._misses = 0
        self._failures: Counter[str] = Counter()

    def run(self, task: Task, adapter: Adapter[T]) -> T | Failure:
        """Result or Failure of ``task``; ``job_key`` is stamped on models that carry it."""
        key = self.store.key(task)
        with self.store.lock(key):
            result = self.store.load(key, adapter.result_type)
            hit = result is not None
            if result is None:
                result = _stamp(self._execute(key, task, adapter), key)
                self.store.save(key, result)
        self._record(hit, result)
        return result

    def stats(self) -> JobStats:
        with self._lock:
            return JobStats(
                hits=self._hits, misses=self._misses, failures_by_kind=dict(self._failures)
            )

    def _execute(self, key: str, task: Task, adapter: Adapter[T]) -> T | Failure:
        """Attempts along LADDER. An earlier attempt that left no result (a stopped run, a
        timeout, an environment failure) is continued once as after a timeout (from its last
        frame, when the adapter has one); otherwise the task starts afresh."""
        index = self.store.begin(key, task)
        used: Counter[FailureKind] = Counter()
        current = task
        if index and not self.store.has_result(key):
            stopped = Failure(kind=FailureKind.TIMEOUT, reason="no result")
            current = adapter.continuation(task, self.store.attempt_dir(key, index - 1),
                                           stopped) or task
        while True:
            current = replace(current, execution=_shared(current.execution))
            workdir = self.store.attempt_dir(key, index)
            outcome = self._attempt(current, adapter, workdir)
            if not isinstance(outcome, Failure):
                return outcome
            following = _continue(current, adapter, workdir, outcome, used)
            if following is None:
                return outcome
            current, index = following, index + 1

    def _attempt(self, task: Task, adapter: Adapter[T], workdir: Path) -> T | Failure:
        workdir.mkdir(parents=True, exist_ok=True)
        cmd = adapter.prepare(task, workdir)
        with self._semaphore.reserve(task.execution.ranks * task.execution.threads):
            try:
                result = run_command(cmd, timeout_s=task.execution.timeout_s)
            except (FileNotFoundError, PermissionError) as exc:
                return Failure(kind=FailureKind.EXECUTABLE_MISSING, reason=f"{cmd.argv[0]}: {exc}")
        outcome = adapter.parse(task, workdir, result)  # an exception is the item's (contain)
        if isinstance(outcome, Failure) and outcome.kind is FailureKind.NONZERO_EXIT and (
                _out_of_memory(result)):
            return outcome.model_copy(update={"kind": FailureKind.OUT_OF_MEMORY})
        return outcome

    def _record(self, hit: bool, result: BaseModel) -> None:
        with self._lock:
            if hit:
                self._hits += 1
            else:
                self._misses += 1
            if isinstance(result, Failure):
                self._failures[str(result.kind)] += 1


def _out_of_memory(result: CommandResult) -> bool:
    """Killed by SIGKILL (the kernel's OOM killer, a scheduler's memory limit) or an allocation
    failure reported at the end of stdout or stderr."""
    if result.returncode in (-9, 137):
        return True
    for path in (result.stdout, result.stderr):
        with path.open("rb") as stream:
            stream.seek(max(path.stat().st_size - _TAIL_BYTES, 0))
            if _OUT_OF_MEMORY.search(stream.read().decode("utf-8", errors="replace")):
                return True
    return False


def _shared(execution: ExecutionSpec) -> ExecutionSpec:
    """At most _RANK_SHARE ranks, with the timeout stretched by the same factor so that a job
    that ends within its timeout at full ranks also ends at the share (same core-seconds)."""
    share = _RANK_SHARE.get()
    if share is None or share >= execution.ranks:
        return execution
    return execution.model_copy(
        update={"ranks": share, "timeout_s": execution.timeout_s * execution.ranks / share})


def _continue(
    task: Task, adapter: Adapter[T], workdir: Path, failure: Failure, used: Counter[FailureKind]
) -> Task | None:
    if used[failure.kind] >= LADDER.get(failure.kind, 0):
        return None
    following = adapter.continuation(task, workdir, failure)
    if following is not None:
        used[failure.kind] += 1
    return following


def _stamp(result: T | Failure, key: str) -> T | Failure:
    model: type[BaseModel] = type(result)
    if "job_key" in model.model_fields:
        return result.model_copy(update={"job_key": key})
    return result


def thread_map(
    fn: Callable[[ItemT], ResultT], items: Sequence[ItemT], *, workers: int
) -> list[ResultT]:
    """``[fn(x) for x in items]`` on up to ``workers`` (the cores) threads, in input order.

    The items run in waves of up to ``workers``, and the jobs of an item run with at most
    workers // (items in its wave) MPI ranks: independent small jobs run side by side instead
    of each on all ranks, and a short last wave (9 items on 4 cores: 4, 4, 1) gets the cores
    it leaves free. Serial code keeps full ranks. An exception (or a stop) cancels the items
    not yet started; the running ones end first."""
    if workers <= 1 or len(items) <= 1:
        return [fn(item) for item in items]
    at_once = min(workers, len(items))

    def shared(indexed: tuple[int, ItemT]) -> ResultT:
        i, item = indexed  # set in the pool's threads, which end with this call
        _RANK_SHARE.set(workers // min(at_once, len(items) - i // at_once * at_once))
        return fn(item)

    pool = ThreadPoolExecutor(max_workers=at_once)
    try:
        return list(pool.map(shared, enumerate(items)))
    finally:
        pool.shutdown(cancel_futures=True)

"""Tasks, adapters and the JobRunner: cache, attempt ladder, core semaphore and the rank share
of concurrent items (design §7.1)."""

from __future__ import annotations

import os
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
from hfauto.core.method import MIN_ATTEMPT_S, Deadline, ExecutionSpec
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

    def run(
        self, task: Task, adapter: Adapter[T], *, deadline: Deadline | None = None
    ) -> T | Failure:
        """Result or Failure of ``task``; ``job_key`` is stamped on models that carry it."""
        key = self.store.key(task)
        with self.store.lock(key):
            result = self.store.load(key, adapter.result_type)
            hit = result is not None
            if result is None:
                result = _stamp(self._execute(key, task, adapter, deadline), key)
                self.store.save(key, result)
        self._record(hit, result)
        return result

    def stats(self) -> JobStats:
        with self._lock:
            return JobStats(
                hits=self._hits, misses=self._misses, failures_by_kind=dict(self._failures)
            )

    def _execute(
        self, key: str, task: Task, adapter: Adapter[T], deadline: Deadline | None
    ) -> T | Failure:
        index = self.store.begin(key, task)
        used: Counter[FailureKind] = Counter()
        current = task
        while True:
            current = replace(current, execution=_shared(current.execution))
            timeout_s = _attempt_timeout(current.execution, deadline)
            if timeout_s is None:
                return Failure(
                    kind=FailureKind.BUDGET_EXHAUSTED,
                    reason=f"less than {MIN_ATTEMPT_S:.0f} s of the deadline left",
                )
            workdir = self.store.attempt_dir(key, index)
            outcome = self._attempt(current, adapter, workdir, timeout_s)
            if not isinstance(outcome, Failure):
                return outcome
            following = _continue(current, adapter, workdir, outcome, used)
            if following is None:
                return outcome
            current, index = following, index + 1

    def _attempt(
        self, task: Task, adapter: Adapter[T], workdir: Path, timeout_s: float
    ) -> T | Failure:
        workdir.mkdir(parents=True, exist_ok=True)
        cmd = adapter.prepare(task, workdir)
        with self._semaphore.reserve(task.execution.ranks * task.execution.threads):
            try:
                result = run_command(cmd, timeout_s=timeout_s)
            except (FileNotFoundError, PermissionError) as exc:
                return Failure(kind=FailureKind.EXECUTABLE_MISSING, reason=f"{cmd.argv[0]}: {exc}")
        try:
            return adapter.parse(task, workdir, result)
        except Exception as exc:  # an unforeseen output fails this job only
            if os.environ.get("HFAUTO_STRICT") == "1":
                raise
            first = (str(exc).splitlines() or [""])[0][:200]
            return Failure(kind=FailureKind.INCOMPLETE_OUTPUT,
                           reason=f"parse:{type(exc).__name__}: {first}")

    def _record(self, hit: bool, result: BaseModel) -> None:
        with self._lock:
            if hit:
                self._hits += 1
            else:
                self._misses += 1
            if isinstance(result, Failure):
                self._failures[str(result.kind)] += 1


def _shared(execution: ExecutionSpec) -> ExecutionSpec:
    """At most _RANK_SHARE ranks, with the timeout stretched by the same factor so that a job
    that ends within its timeout at full ranks also ends at the share (same core-seconds)."""
    share = _RANK_SHARE.get()
    if share is None or share >= execution.ranks:
        return execution
    return execution.model_copy(
        update={"ranks": share, "timeout_s": execution.timeout_s * execution.ranks / share})


def _attempt_timeout(execution: ExecutionSpec, deadline: Deadline | None) -> float | None:
    if deadline is None:
        return execution.timeout_s
    if deadline.expired():
        return None
    return min(execution.timeout_s, deadline.remaining())


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
    it leaves free. Serial code keeps full ranks."""
    if workers <= 1 or len(items) <= 1:
        return [fn(item) for item in items]
    at_once = min(workers, len(items))

    def shared(indexed: tuple[int, ItemT]) -> ResultT:
        i, item = indexed  # set in the pool's threads, which end with this call
        _RANK_SHARE.set(workers // min(at_once, len(items) - i // at_once * at_once))
        return fn(item)

    with ThreadPoolExecutor(max_workers=at_once) as pool:
        return list(pool.map(shared, enumerate(items)))

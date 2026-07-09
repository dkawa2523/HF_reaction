from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any

from hfauto.core.executables import resolve_executable, run_command
from hfauto.core.io import ensure_dir, write_json


@dataclass
class SchedulerResult:
    action: str
    scheduler: str
    command: list[str]
    dry_run: bool
    ok: bool
    returncode: int | None = None
    stdout_path: str | None = None
    stderr_path: str | None = None
    parsed_job_id: str | None = None
    note: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


class SchedulerAdapter:
    name = "base"
    submit_executable = ""
    status_executable = ""
    cancel_executable = ""

    def submit_command(self, script: str | Path) -> list[str]:
        raise NotImplementedError

    def status_command(self, job_id: str | None = None, user: str | None = None) -> list[str]:
        raise NotImplementedError

    def cancel_command(self, job_id: str) -> list[str]:
        raise NotImplementedError

    def parse_submit_job_id(self, stdout: str) -> str | None:
        return stdout.strip().splitlines()[-1].strip() if stdout.strip() else None

    def run_scheduler_command(
        self,
        action: str,
        command: list[str],
        workdir: str | Path,
        dry_run: bool = True,
        allow_execute: bool = False,
        timeout_s: int = 120,
    ) -> SchedulerResult:
        ensure_dir(workdir)
        if dry_run or not allow_execute:
            return SchedulerResult(action=action, scheduler=self.name, command=command, dry_run=True, ok=True, note="dry_run_no_scheduler_command_executed")
        exe = resolve_executable(command[0])
        if not exe:
            return SchedulerResult(action=action, scheduler=self.name, command=command, dry_run=False, ok=False, note=f"executable_not_found:{command[0]}")
        result = run_command(command, cwd=workdir, timeout_s=timeout_s, stdout_name=f"{action}_stdout.txt", stderr_name=f"{action}_stderr.txt")
        stdout = Path(result.stdout_path).read_text(encoding="utf-8") if Path(result.stdout_path).exists() else ""
        parsed = self.parse_submit_job_id(stdout) if action == "submit" else None
        return SchedulerResult(action=action, scheduler=self.name, command=command, dry_run=False, ok=result.ok, returncode=result.returncode, stdout_path=result.stdout_path, stderr_path=result.stderr_path, parsed_job_id=parsed)


def write_scheduler_result(path: str | Path, result: SchedulerResult) -> Path:
    return write_json(path, result.to_dict())

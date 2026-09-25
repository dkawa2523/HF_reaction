from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

from hfauto.core import executables


def test_command_output_is_streamed_to_evidence_files(tmp_path: Path) -> None:
    result = executables.run_command(
        [
            sys.executable,
            "-c",
            "import sys; print('stdout evidence'); print('stderr evidence', file=sys.stderr)",
        ],
        tmp_path,
    )

    assert result.ok is True
    assert (tmp_path / "stdout.txt").read_text(encoding="utf-8").strip() == (
        "stdout evidence"
    )
    assert (tmp_path / "stderr.txt").read_text(encoding="utf-8").strip() == (
        "stderr evidence"
    )


class _TimeoutProcess:
    def __init__(self, stdout: Any = None, stderr: Any = None) -> None:
        self.stdout = stdout
        self.stderr = stderr
        self.pid = 987654
        self.returncode = None

    def wait(self, timeout: float | None = None):
        raise subprocess.TimeoutExpired(
            cmd=["fake", "arg"],
            timeout=timeout,
        )


def test_timeout_payloads_are_persisted_with_structured_sidecar(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    process = _TimeoutProcess()
    popen_calls: list[tuple[list[str], dict[str, Any]]] = []
    terminated: list[_TimeoutProcess] = []

    def fake_popen(command: list[str], **kwargs: Any) -> _TimeoutProcess:
        popen_calls.append((command, kwargs))
        kwargs["stdout"].write("partial stdout")
        kwargs["stderr"].write("partial stderr")
        return process

    def fake_terminate(target: _TimeoutProcess) -> None:
        terminated.append(target)

    monkeypatch.setattr(executables.subprocess, "Popen", fake_popen)
    monkeypatch.setattr(executables, "_terminate_process_tree", fake_terminate)

    result = executables.run_command(
        ["fake", "arg"],
        tmp_path,
        timeout_s=7,
        stdout_name="captured.out",
        stderr_name="captured.err",
    )

    assert result.timed_out is True
    assert result.returncode == executables.TIMEOUT_RETURNCODE == 124
    assert result.ok is False
    assert terminated == [process]
    assert (tmp_path / "captured.out").read_text(encoding="utf-8") == "partial stdout"
    assert (tmp_path / "captured.err").read_text(encoding="utf-8") == "partial stderr"

    sidecar = json.loads((tmp_path / "command_result.json").read_text(encoding="utf-8"))
    assert sidecar == result.to_dict()
    assert sidecar["command"] == ["fake", "arg"]
    assert sidecar["returncode"] == 124
    assert sidecar["timed_out"] is True
    assert sidecar["stdout_path"] == str(tmp_path / "captured.out")
    assert sidecar["stderr_path"] == str(tmp_path / "captured.err")
    assert sidecar["executable"] == "fake"

    assert len(popen_calls) == 1
    command, kwargs = popen_calls[0]
    assert command == ["fake", "arg"]
    assert kwargs["text"] is True
    assert kwargs["encoding"] == "utf-8"
    assert kwargs["errors"] == "replace"
    if os.name == "posix":
        assert kwargs["start_new_session"] is True
    elif os.name == "nt" and hasattr(subprocess, "CREATE_NEW_PROCESS_GROUP"):
        assert kwargs["creationflags"] == subprocess.CREATE_NEW_PROCESS_GROUP


def test_timeout_termination_does_not_overwrite_streamed_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    process = _TimeoutProcess()

    def fake_popen(*_args: Any, **kwargs: Any) -> _TimeoutProcess:
        kwargs["stdout"].write("prefix")
        kwargs["stderr"].write("warning")
        return process

    monkeypatch.setattr(executables.subprocess, "Popen", fake_popen)
    monkeypatch.setattr(
        executables,
        "_terminate_process_tree",
        lambda _process: None,
    )

    executables.run_command(["fake"], tmp_path, timeout_s=1)

    assert (tmp_path / "stdout.txt").read_text(encoding="utf-8") == "prefix"
    assert (tmp_path / "stderr.txt").read_text(encoding="utf-8") == "warning"


@pytest.mark.skipif(os.name != "posix", reason="POSIX process-group contract")
def test_posix_timeout_cleanup_signals_owned_group_and_reaps_launcher(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class Process:
        pid = 24680
        returncode: int | None = None
        stdout = None
        stderr = None

        def wait(self, timeout: float | None = None):
            assert timeout == executables._TERMINATION_GRACE_S
            self.returncode = -9
            return self.returncode

    process = Process()
    signals: list[int] = []
    monkeypatch.setattr(
        executables,
        "_signal_posix_process_group",
        lambda target, sig: signals.append(sig) if target is process else None,
    )

    result = executables._terminate_process_tree(process)  # type: ignore[arg-type]

    assert result is None
    assert signals == [signal.SIGKILL]

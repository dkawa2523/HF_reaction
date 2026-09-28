"""Shared fixtures (design §10.1). Test doubles live in tests/fakes.py (pythonpath = tests)."""

import pytest


@pytest.fixture(autouse=True)
def strict(monkeypatch):  # unforeseen exceptions re-raise; a containment test deletes it
    monkeypatch.setenv("HFAUTO_STRICT", "1")


@pytest.fixture
def tmp_run(tmp_path):  # the run directory; FileRefs of the fake engines are relative to it
    (tmp_path / "run").mkdir()
    return tmp_path / "run"

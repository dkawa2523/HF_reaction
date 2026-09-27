"""Shared fixtures (design §10.1). Test doubles live in tests/fakes.py (pythonpath = tests)."""

from contextlib import ExitStack

import pytest

from hfauto.backends import engines


@pytest.fixture(autouse=True)
def strict(monkeypatch):  # unforeseen exceptions re-raise; a containment test deletes it
    monkeypatch.setenv("HFAUTO_STRICT", "1")


@pytest.fixture
def override_engine():  # override_engine(capability, name, engine) until the test ends
    with ExitStack() as stack:
        yield lambda cap, name, engine: stack.enter_context(
            engines.override(cap, name, lambda **_: engine))


@pytest.fixture
def tmp_run(tmp_path):  # the run directory; FileRefs of the fake engines are relative to it
    (tmp_path / "run").mkdir()
    return tmp_path / "run"

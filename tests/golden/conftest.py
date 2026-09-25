"""The `golden` fixture: read excerpts of real engine outputs from tests/golden/data (design §10.2)."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

DATA = Path(__file__).resolve().parent / "data"


class Golden:
    """Resolve data-relative paths; `relpath` may name a stored `.gz` with or without the suffix."""

    def __init__(self, tmp: Path) -> None:
        self._tmp = tmp

    def _stored(self, relpath: str) -> Path:
        plain = DATA / relpath
        return plain if plain.is_file() else DATA / f"{relpath.removesuffix('.gz')}.gz"

    def path(self, relpath: str) -> Path:
        """Return a readable file; a gzipped excerpt is expanded into a temporary directory."""
        stored = self._stored(relpath)
        if stored.suffix != ".gz":
            return stored
        out = self._tmp / stored.relative_to(DATA).with_suffix("")
        if not out.is_file():
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_bytes(gzip.decompress(stored.read_bytes()))
        return out

    def text(self, relpath: str) -> str:
        return self.path(relpath).read_text(encoding="utf-8")


@pytest.fixture(scope="session")
def golden(tmp_path_factory: pytest.TempPathFactory) -> Golden:
    return Golden(tmp_path_factory.mktemp("golden"))

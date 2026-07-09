from __future__ import annotations

from typing import Protocol

from hfauto.core.schemas.artifact import Artifact


class DBProvider(Protocol):
    name: str

    def enrich(self, molecule_artifact: Artifact) -> dict:
        ...

from __future__ import annotations

from pathlib import Path
from typing import Protocol

from hfauto.backends.ts.base import IRCResult, TSSearchResult
from hfauto.core.schemas.artifact import Artifact


class QMBackend(Protocol):
    name: str

    def optimize_frequency(self, species: Artifact, method: dict, workdir: str) -> Artifact:
        ...

    def single_point(self, species: Artifact, method: dict, workdir: str) -> Artifact:
        ...


class TSBackend(Protocol):
    name: str

    def search_ts(
        self,
        reaction: Artifact,
        reactant: Artifact,
        product: Artifact,
        method: dict,
        workdir: str | Path,
    ) -> TSSearchResult:
        ...

    def run_irc(
        self,
        reaction: Artifact,
        ts_species: Artifact,
        reactant: Artifact,
        product: Artifact,
        method: dict,
        workdir: str | Path,
    ) -> IRCResult:
        ...

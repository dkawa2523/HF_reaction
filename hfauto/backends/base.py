from __future__ import annotations

from typing import Protocol

from hfauto.core.schemas.artifact import Artifact


class QMBackend(Protocol):
    name: str

    def optimize_frequency(self, species: Artifact, method: dict, workdir: str) -> Artifact:
        ...

    def single_point(self, species: Artifact, method: dict, workdir: str) -> Artifact:
        ...


class TSBackend(Protocol):
    name: str

    def search(self, reaction: Artifact, reactant: Artifact, product: Artifact, method: dict, workdir: str) -> list[Artifact]:
        ...

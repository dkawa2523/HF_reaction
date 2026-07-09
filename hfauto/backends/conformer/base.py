from __future__ import annotations

from typing import Protocol

from hfauto.core.schemas.artifact import Artifact


class ConformerBackend(Protocol):
    name: str

    def generate(self, molecule: Artifact, config: dict, workdir: str) -> list[Artifact]:
        ...

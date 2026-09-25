from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from hfauto.core.schemas.manifest import Manifest


@dataclass
class StageContext:
    out_dir: Path
    run_id: str
    global_config: dict[str, Any]


class Stage:
    name: str
    accepts_empty_manifest: bool = False

    def run(self, manifest: Manifest | None, config: dict[str, Any], context: StageContext) -> Manifest:
        raise NotImplementedError

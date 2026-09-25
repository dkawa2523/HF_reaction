from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class VizRecord:
    artifact_id: str
    viz_type: str
    path: str
    parents: list[str] = field(default_factory=list)
    renderer: str | None = None
    data: dict[str, Any] = field(default_factory=dict)

@dataclass
class VizManifest:
    run_id: str
    out_dir: Path
    records: list[VizRecord] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    def add(self, record: VizRecord) -> None: self.records.append(record)
    def warn(self, message: str) -> None: self.warnings.append(message)
    def write(self) -> Path:
        self.out_dir.mkdir(parents=True, exist_ok=True)
        path = self.out_dir / "viz_manifest.json"
        path.write_text(json.dumps({"schema_version":"hfauto_viz.manifest.v1","run_id":self.run_id,"warnings":self.warnings,"visualizations":[asdict(r) for r in self.records]}, ensure_ascii=False, indent=2), encoding="utf-8")
        return path

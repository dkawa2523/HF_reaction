from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, Any

from hfauto.core.schemas.manifest import Manifest


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def read_manifest(path: str | Path) -> Manifest:
    return Manifest.model_validate_json(Path(path).read_text(encoding="utf-8"))


def write_manifest(manifest: Manifest, out_dir: str | Path) -> Path:
    out = ensure_dir(out_dir) / "manifest.json"
    out.write_text(manifest.model_dump_json(indent=2), encoding="utf-8")
    return out


def write_json(path: str | Path, data: Any) -> Path:
    p = Path(path)
    ensure_dir(p.parent)
    p.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    return p


def read_json(path: str | Path) -> Any:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_jsonl(records: Iterable[dict], path: str | Path) -> Path:
    p = Path(path)
    ensure_dir(p.parent)
    with p.open("w", encoding="utf-8") as f:
        for record in records:
            f.write(json.dumps(record, ensure_ascii=False) + "\n")
    return p


def read_jsonl(path: str | Path) -> list[dict]:
    p = Path(path)
    if not p.exists():
        return []
    return [json.loads(line) for line in p.read_text(encoding="utf-8").splitlines() if line.strip()]

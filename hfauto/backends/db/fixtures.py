from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def load_fixture(fixture_path: str | None) -> dict[str, Any]:
    if not fixture_path:
        return {}
    p = Path(fixture_path)
    if not p.exists():
        return {}
    return json.loads(p.read_text(encoding="utf-8"))


def fixture_by_name(fixture: dict[str, Any], provider: str, name: str | None) -> dict[str, Any] | None:
    if not fixture or not name:
        return None
    block = fixture.get(provider, {}) if isinstance(fixture.get(provider), dict) else {}
    by_name = block.get("by_name", {}) if isinstance(block.get("by_name"), dict) else {}
    return by_name.get(str(name).strip().lower())

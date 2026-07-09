from __future__ import annotations

import json
import sqlite3
import time
from pathlib import Path
from typing import Any, Optional


class SQLiteJSONCache:
    """Tiny JSON cache for public database responses.

    The cache is deliberately simple so it can be reviewed and replaced easily.
    Keys are provider-specific strings such as ``pubchem:inchikey:XXXX``.
    """

    def __init__(self, path: str | Path | None):
        self.path = Path(path) if path else None
        if self.path:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            with sqlite3.connect(self.path) as con:
                con.execute(
                    "CREATE TABLE IF NOT EXISTS cache (provider TEXT, key TEXT, value TEXT, created_at REAL, PRIMARY KEY(provider, key))"
                )

    def get(self, provider: str, key: str) -> Optional[dict[str, Any]]:
        if not self.path:
            return None
        with sqlite3.connect(self.path) as con:
            row = con.execute(
                "SELECT value FROM cache WHERE provider=? AND key=?", (provider, key)
            ).fetchone()
        if row is None:
            return None
        return json.loads(row[0])

    def put(self, provider: str, key: str, value: dict[str, Any]) -> None:
        if not self.path:
            return
        with sqlite3.connect(self.path) as con:
            con.execute(
                "INSERT OR REPLACE INTO cache(provider, key, value, created_at) VALUES (?, ?, ?, ?)",
                (provider, key, json.dumps(value, ensure_ascii=False), time.time()),
            )

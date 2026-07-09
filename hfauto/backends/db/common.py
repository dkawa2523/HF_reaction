from __future__ import annotations

import json
import re
import time
from pathlib import Path
from typing import Any

from hfauto.backends.db.cache import SQLiteJSONCache


def norm_key(value: Any) -> str:
    return re.sub(r"\s+", "_", str(value or "").strip().lower())


def load_fixture(fixture_path: str | Path | None) -> dict[str, Any]:
    if not fixture_path:
        return {}
    path = Path(fixture_path)
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def lookup_fixture(provider: str, molecule_data: dict[str, Any], fixture: dict[str, Any]) -> dict[str, Any] | None:
    section = fixture.get(provider, {}) if fixture else {}
    if not isinstance(section, dict):
        return None
    keys = []
    for field in ["inchikey", "canonical_smiles", "name"]:
        if molecule_data.get(field):
            keys.append(("by_" + field, norm_key(molecule_data.get(field))))
    # User fixtures often key names as they appear in SDF, with spaces or underscores.
    if molecule_data.get("name"):
        n = str(molecule_data["name"])
        keys.extend([("by_name", n.strip().lower()), ("by_name", norm_key(n)), ("by_name", n.replace("_", " ").lower())])
    for bucket, key in keys:
        table = section.get(bucket, {})
        if key in table:
            result = dict(table[key])
            result.setdefault("matched", True)
            result.setdefault("from_fixture", True)
            return result
        # try normalized view of fixture keys too
        for raw_key, value in table.items():
            if norm_key(raw_key) == key:
                result = dict(value)
                result.setdefault("matched", True)
                result.setdefault("from_fixture", True)
                return result
    return None


class CachedHTTPMixin:
    """Small helper for optional cache-first public DB providers.

    Providers are offline by default.  Network access must be explicitly enabled
    with allow_network=True so tests and reviews are reproducible.
    """

    name: str

    def _setup_http(self, cache_path: str | None = None, fixture_path: str | None = None, allow_network: bool = False, rate_limit_per_sec: float = 2.0, **_: Any) -> None:
        self.cache = SQLiteJSONCache(cache_path)
        self.fixture = load_fixture(fixture_path)
        self.allow_network = bool(allow_network)
        self.rate_limit_per_sec = float(rate_limit_per_sec)
        self._last_request = 0.0
        self.fixture_path = fixture_path

    def _cache_get(self, key: str) -> dict[str, Any] | None:
        return self.cache.get(self.name, key) if getattr(self, "cache", None) else None

    def _cache_put(self, key: str, value: dict[str, Any]) -> None:
        if getattr(self, "cache", None):
            self.cache.put(self.name, key, value)

    def _get_json(self, url: str, params: dict[str, Any] | None = None, headers: dict[str, str] | None = None, timeout: int = 45) -> dict[str, Any]:
        if not getattr(self, "allow_network", False):
            raise RuntimeError(f"{self.name} network access disabled; set allow_network: true")
        try:
            import requests
        except Exception as exc:  # pragma: no cover
            raise RuntimeError("requests is not installed") from exc
        wait = max(0.0, (1.0 / max(self.rate_limit_per_sec, 0.1)) - (time.time() - self._last_request))
        if wait:
            time.sleep(wait)
        self._last_request = time.time()
        r = requests.get(url, params=params, headers=headers, timeout=timeout)
        r.raise_for_status()
        return r.json()

    def _get_text(self, url: str, params: dict[str, Any] | None = None, headers: dict[str, str] | None = None, timeout: int = 45) -> str:
        if not getattr(self, "allow_network", False):
            raise RuntimeError(f"{self.name} network access disabled; set allow_network: true")
        try:
            import requests
        except Exception as exc:  # pragma: no cover
            raise RuntimeError("requests is not installed") from exc
        wait = max(0.0, (1.0 / max(self.rate_limit_per_sec, 0.1)) - (time.time() - self._last_request))
        if wait:
            time.sleep(wait)
        self._last_request = time.time()
        r = requests.get(url, params=params, headers=headers, timeout=timeout)
        r.raise_for_status()
        return r.text


def merge_provider_status(data: dict[str, Any], provider: str, status: dict[str, Any]) -> dict[str, Any]:
    data.setdefault("db_provider_status", {})[provider] = status
    return data

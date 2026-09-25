from __future__ import annotations

"""Small, reviewable HTTP/cache utilities for public-data providers.

This module is intentionally tiny. Network access is opt-in through either
``allow_network`` or
``network_enabled``; otherwise cache misses are returned as data instead of
exceptions so offline pipelines remain reproducible.
"""

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from hfauto.backends.db.cache import SQLiteJSONCache


@dataclass
class HTTPResult:
    ok: bool
    status_code: int | None = None
    text: str | None = None
    json_data: dict[str, Any] | list[Any] | None = None
    from_cache: bool = False
    error: str | None = None
    url: str | None = None

    @property
    def reason(self) -> str | None:  # compatibility alias
        return self.error


class CachedHTTPClient:
    def __init__(
        self,
        provider: str,
        cache_path: str | Path | None = None,
        allow_network: bool | None = None,
        rate_limit_per_sec: float = 2.0,
        timeout_s: float = 30.0,
        user_agent: str = "hfauto/0.14 public-data-provider",
        network_enabled: bool | None = None,
        **_: Any,
    ) -> None:
        self.provider = provider
        self.cache = SQLiteJSONCache(cache_path)
        self.allow_network = bool(network_enabled if network_enabled is not None else allow_network)
        self.rate_limit_per_sec = float(rate_limit_per_sec)
        self.timeout_s = float(timeout_s)
        self.user_agent = user_agent
        self._last_request = 0.0

    @staticmethod
    def _full_url(url: str, params: dict[str, Any] | None = None) -> str:
        if not params:
            return url
        return f"{url}?{urlencode({k: v for k, v in params.items() if v is not None})}"

    def get_text(self, url: str, params: dict[str, Any] | None = None, cache_key: str | None = None) -> HTTPResult:
        full_url = self._full_url(url, params)
        key = cache_key or full_url
        cached = self.cache.get(self.provider, key)
        if cached is not None:
            return HTTPResult(True, cached.get("status_code"), cached.get("text"), cached.get("json"), True, None, full_url)
        if not self.allow_network:
            return HTTPResult(False, None, None, None, False, "network_disabled_and_cache_miss", full_url)
        wait = max(0.0, (1.0 / max(self.rate_limit_per_sec, 0.1)) - (time.time() - self._last_request))
        if wait:
            time.sleep(wait)
        self._last_request = time.time()
        try:
            req = Request(full_url, headers={"User-Agent": self.user_agent})
            with urlopen(req, timeout=self.timeout_s) as resp:
                raw = resp.read().decode("utf-8", errors="replace")
                status = getattr(resp, "status", 200)
            self.cache.put(self.provider, key, {"status_code": status, "text": raw, "json": None, "url": full_url})
            return HTTPResult(True, status, raw, None, False, None, full_url)
        except Exception as exc:  # pragma: no cover - depends on network
            return HTTPResult(False, None, None, None, False, str(exc), full_url)

    def get_json(self, url: str, params: dict[str, Any] | None = None, cache_key: str | None = None) -> HTTPResult:
        full_url = self._full_url(url, params)
        key = cache_key or full_url
        cached = self.cache.get(self.provider, key)
        if cached is not None:
            return HTTPResult(True, cached.get("status_code"), cached.get("text"), cached.get("json"), True, None, full_url)
        text_result = self.get_text(url, params=params, cache_key=key)
        if not text_result.ok or text_result.text is None:
            return text_result
        try:
            data = json.loads(text_result.text)
        except Exception as exc:
            return HTTPResult(False, text_result.status_code, text_result.text, None, text_result.from_cache, f"json_parse_failed: {exc}", full_url)
        self.cache.put(self.provider, key, {"status_code": text_result.status_code, "text": text_result.text, "json": data, "url": full_url})
        return HTTPResult(True, text_result.status_code, text_result.text, data, text_result.from_cache, None, full_url)

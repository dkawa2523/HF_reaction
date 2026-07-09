from __future__ import annotations

import re
from typing import Any

from hfauto.backends.db.common import load_fixture, lookup_fixture
from hfauto.backends.db.http import CachedHTTPClient
from hfauto.core.reference_data import load_reference_set, match_reference


def _parse_first_number_after(label: str, text: str) -> float | None:
    pattern = rf"{re.escape(label)}[^\n\r]*?([-+]?\d+(?:\.\d+)?)"
    m = re.search(pattern, text, flags=re.IGNORECASE)
    if not m:
        return None
    try:
        return float(m.group(1))
    except Exception:
        return None


def parse_nist_ion_energetics_text(text: str) -> dict[str, Any]:
    flat = re.sub(r"<[^>]+>", " ", text)
    flat = re.sub(r"\s+", " ", flat)
    out: dict[str, Any] = {}
    pa = _parse_first_number_after("Proton affinity", flat)
    gb = _parse_first_number_after("Gas basicity", flat)
    if pa is not None:
        out["proton_affinity_kj_mol"] = pa
    if gb is not None:
        out["gas_basicity_kj_mol"] = gb
    out["has_ir"] = "Infrared Spectrum" in text or "Gas Phase IR Spectrum" in text
    out["has_ion_energetics"] = "Proton affinity" in flat or "Gas basicity" in flat or "Ion Energetics" in flat
    return out


class NISTWebBookProvider:
    """NIST Chemistry WebBook PA/GB/IR provider.

    The WebBook is HTML rather than a stable JSON API, so online parsing is
    conservative and cached.  Offline fixture/bundled reference-set matches are
    supported for reproducible method-validation tests.
    """

    name = "nist_webbook"
    BASE = "https://webbook.nist.gov/cgi/cbook.cgi"

    def __init__(
        self,
        cache_path: str | None = None,
        fixture_path: str | None = None,
        network_enabled: bool = False,
        allow_network: bool | None = None,
        rate_limit_per_sec: float = 1.0,
        timeout_s: int = 30,
        reference_set_path: str | None = None,
        use_builtin_references: bool = True,
        **_: Any,
    ):
        if allow_network is not None:
            network_enabled = bool(allow_network)
        self.client = CachedHTTPClient(self.name, cache_path or ".hfauto_cache/nist_webbook.sqlite", network_enabled=network_enabled, rate_limit_per_sec=rate_limit_per_sec, timeout_s=timeout_s)
        self.fixture = load_fixture(fixture_path)
        self.reference_set = load_reference_set(reference_set_path) if use_builtin_references else {"references": []}

    def _fixture_or_reference(self, data: dict[str, Any]) -> dict[str, Any] | None:
        hit = lookup_fixture(self.name, data, self.fixture) if self.fixture else None
        if hit:
            result = dict(hit)
            result.setdefault("matched", True)
            result.setdefault("provider_status", "success_fixture")
            result.setdefault("source", "offline_fixture")
            return result
        ref = match_reference(data, self.reference_set)
        if ref:
            props = ref.get("properties") or {}
            if any(k in props for k in ["proton_affinity_kj_mol", "gas_basicity_kj_mol", "hf_stretch_cm1", "dipole_D"]):
                return {
                    "matched": True,
                    "provider_status": "builtin_reference_match",
                    "match_type": ref.get("match_type"),
                    "match_confidence": ref.get("match_confidence"),
                    "reference_id": ref.get("ref_id"),
                    "name": ref.get("name"),
                    **props,
                    "sources": ref.get("sources", []),
                    "source": "bundled_hf_amine_reference_set",
                    "note": "Use production NIST retrieval or curated project references for final reporting.",
                }
        return None

    def _query(self, data: dict[str, Any]) -> dict[str, Any]:
        params = {"Units": "SI", "Mask": "20"}
        if data.get("name"):
            params["Name"] = data.get("name")
        elif data.get("formula"):
            params["Formula"] = data.get("formula")
        else:
            return {"matched": False, "provider_status": "skipped", "reason": "missing_name_or_formula"}
        r = self.client.get_text(self.BASE, params=params)
        if not r.ok or not r.text:
            return {"matched": False, "provider_status": "skipped_or_failed", "reason": r.reason, "url": r.url}
        parsed = parse_nist_ion_energetics_text(r.text)
        return {"matched": bool(parsed), "provider_status": "success", "url": r.url, **parsed, "source": "NIST Chemistry WebBook"}

    def enrich(self, molecule_artifact):
        data = molecule_artifact.data.copy()
        data.setdefault("public_data", {})
        result = self._fixture_or_reference(data)
        if result is None:
            result = self._query(data)
        if not result.get("matched"):
            result.setdefault("intended_use", "PA/GB, IR, small-molecule method validation")
        data["public_data"]["nist_webbook"] = result
        if result.get("matched"):
            refs = data["public_data"].setdefault("reference_values", {})
            for k in ["proton_affinity_kj_mol", "gas_basicity_kj_mol", "hf_stretch_cm1", "dipole_D"]:
                if result.get(k) is not None:
                    refs[k] = result.get(k)
        data.setdefault("db_provider_status", {})[self.name] = {"status": result.get("provider_status", "unknown"), "matched": bool(result.get("matched"))}
        return data

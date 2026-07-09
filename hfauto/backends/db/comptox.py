from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.backends.db.common import load_fixture, lookup_fixture
from hfauto.backends.db.http import CachedHTTPClient


class CompToxProvider:
    name = "comptox"

    def __init__(self, local_csv: str | None = None, api_url: str | None = None, api_key_env: str = "EPA_CTX_API_KEY", api_key: str | None = None, network_enabled: bool = False, fixture_path: str | None = None, allow_network: bool | None = None, cache_path: str | None = None, cache_dir: str | None = None, **kwargs: Any):
        self.local_csv = local_csv
        self.api_url = api_url
        self.api_key = api_key or os.environ.get(api_key_env)
        self.query_param = kwargs.get("query_param", "q")
        self.api_key_header = kwargs.get("api_key_header", "x-api-key")
        self.fixture = load_fixture(fixture_path)
        if allow_network is not None:
            network_enabled = bool(allow_network)
        self.client = CachedHTTPClient(self.name, cache_path=cache_path or ".hfauto_cache/comptox.sqlite", cache_dir=cache_dir, network_enabled=network_enabled, rate_limit_per_sec=float(kwargs.get("rate_limit_per_sec", 2.0)), user_agent=kwargs.get("user_agent", "hfauto/0.7 CompToxProvider"))

    def _match_local(self, data: dict[str, Any]) -> dict[str, Any] | None:
        if not self.local_csv or not Path(self.local_csv).exists():
            return None
        df = pd.read_csv(self.local_csv)
        keys = [("inchikey", data.get("inchikey")), ("InChIKey", data.get("inchikey")), ("name", data.get("name")), ("preferredName", data.get("name")), ("CASRN", (data.get("sdf_props") or {}).get("CAS")), ("casrn", (data.get("sdf_props") or {}).get("CAS"))]
        masks = []
        for col, value in keys:
            if col in df.columns and value:
                masks.append(df[col].astype(str).str.lower() == str(value).lower())
        if not masks:
            return None
        mask = masks[0]
        for m in masks[1:]:
            mask = mask | m
        hit = df[mask].head(1)
        return None if hit.empty else hit.iloc[0].to_dict()

    def _api_lookup(self, data: dict[str, Any]) -> dict[str, Any] | None:
        if not self.api_url or not self.api_key:
            return None
        query = data.get("inchikey") or data.get("name") or data.get("canonical_smiles")
        if not query:
            return None
        resp = self.client.get_json(self.api_url, params={self.query_param: query}, cache_key=f"api:{query}")
        if not resp.ok:
            return {"provider_status": resp.status, "reason": resp.reason, "matched": False}
        return {"provider_status": "success", "matched": True, "raw": resp.value}

    @staticmethod
    def _float_from(row: dict[str, Any], keys: list[str]) -> float | None:
        for key in keys:
            if key in row and row[key] not in {None, "", "nan"}:
                try:
                    return float(row[key])
                except Exception:
                    continue
        return None

    def enrich(self, molecule_artifact):
        data = molecule_artifact.data.copy()
        public = data.setdefault("public_data", {})
        gas = public.setdefault("gas_process", {})
        fixture_hit = lookup_fixture(self.name, data, self.fixture)
        if fixture_hit is not None:
            result = dict(fixture_hit)
            result.setdefault("provider_status", "success")
            result.setdefault("source", "offline_fixture")
        else:
            row = self._match_local(data)
            if row is not None:
                result = {"matched": True, "provider_status": "success", "source": "CompTox/local CSV", "record": row, **row}
            else:
                api_result = self._api_lookup(data)
                result = {"matched": False, "provider_status": "not_configured"} if api_result is None else {"source": "EPA CTX API/configured endpoint", **api_result}
        public["comptox"] = result
        record = result.get("record") or result
        vapor_pressure = self._float_from(record, ["vapor_pressure_Pa_25C", "Vapor Pressure", "vaporPressurePa", "VP_Pa"])
        boiling_point = self._float_from(record, ["boiling_point_C", "Boiling Point", "boilingPointC", "BP_C"])
        melting_point = self._float_from(record, ["melting_point_C", "Melting Point", "meltingPointC", "MP_C"])
        if vapor_pressure is not None:
            gas["vapor_pressure_Pa_25C"] = vapor_pressure
        if boiling_point is not None:
            gas["boiling_point_C"] = boiling_point
        if melting_point is not None:
            gas["melting_point_C"] = melting_point
        feasibility = result.get("gas_process_feasibility") or "unknown"
        reasons = []
        if feasibility == "unknown" and vapor_pressure is not None:
            if vapor_pressure >= 100.0:
                feasibility = "likely_gas_process_compatible"
            elif vapor_pressure >= 1.0:
                feasibility = "borderline_low_vapor_pressure"; reasons.append("low_vapor_pressure")
            else:
                feasibility = "poor_low_vapor_pressure"; reasons.append("very_low_vapor_pressure")
        if boiling_point is not None and boiling_point > 250:
            reasons.append("high_boiling_point")
            if feasibility == "unknown":
                feasibility = "borderline_high_boiling_point"
        gas["gas_process_feasibility"] = feasibility if feasibility != "unknown" else gas.get("gas_process_feasibility", "unknown_comptox_no_data")
        gas["gas_process_reasons"] = reasons
        gas["ehs_review_flag"] = result.get("ehs_review_flag") or gas.get("ehs_review_flag", "review_if_selected")
        data.setdefault("db_provider_status", {})[self.name] = {"status": result.get("provider_status"), "matched": result.get("matched", False)}
        return data

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.backends.db.common import load_fixture, lookup_fixture


class NIOSHProvider:
    name = "niosh"

    def __init__(self, local_csv: str | None = None, fixture_path: str | None = None, **kwargs: Any):
        self.local_csv = local_csv
        self.fixture = load_fixture(fixture_path)
        self.config = kwargs

    def enrich(self, molecule_artifact):
        data = molecule_artifact.data.copy()
        public = data.setdefault("public_data", {})
        gas = public.setdefault("gas_process", {})
        fixture_hit = lookup_fixture(self.name, data, self.fixture)
        if fixture_hit is not None:
            result = dict(fixture_hit)
            result.setdefault("provider_status", "success")
            result.setdefault("source", "offline_fixture")
            public["niosh"] = result
            review = bool(result.get("manual_ehs_review_required", True))
            gas["ehs_review_flag"] = "manual_ehs_review_required" if review else "no_niosh_review_flag"
            gas["niosh_hazard_summary"] = result.get("hazard_summary")
            data.setdefault("db_provider_status", {})[self.name] = {"status": "success", "matched": True}
            return data
        result = {"matched": False, "provider_status": "not_configured", "source": "NIOSH/local curated table"}
        try:
            if self.local_csv and Path(self.local_csv).exists():
                df = pd.read_csv(self.local_csv)
                props = data.get("sdf_props") or {}
                masks = []
                for col, value in [("name", data.get("name")), ("CAS", props.get("CAS")), ("cas", props.get("cas")), ("CASRN", props.get("CASRN"))]:
                    if col in df.columns and value:
                        masks.append(df[col].astype(str).str.lower() == str(value).lower())
                if masks:
                    mask = masks[0]
                    for m in masks[1:]:
                        mask = mask | m
                    hit = df[mask].head(1)
                    if not hit.empty:
                        row = hit.iloc[0].to_dict()
                        result = {"matched": True, "provider_status": "success", "record": row, "source": "NIOSH/local curated table"}
                        review = bool(row.get("manual_ehs_review_required", True))
                        gas["ehs_review_flag"] = "manual_ehs_review_required" if review else "no_niosh_review_flag"
                        gas["niosh_idlh"] = row.get("idlh") or row.get("IDLH")
                        gas["niosh_rel"] = row.get("rel") or row.get("REL")
                        gas["niosh_pel"] = row.get("pel") or row.get("PEL")
            if not result.get("matched"):
                gas.setdefault("ehs_review_flag", "review_if_selected")
        except Exception as exc:
            result = {"matched": False, "provider_status": "failed", "reason": str(exc), "source": "NIOSH/local curated table"}
        public["niosh"] = result
        data.setdefault("db_provider_status", {})[self.name] = {"status": result.get("provider_status"), "matched": result.get("matched", False)}
        return data

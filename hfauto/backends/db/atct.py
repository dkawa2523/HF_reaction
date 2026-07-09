from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.backends.db.common import load_fixture, lookup_fixture


class ATcTProvider:
    name = "atct"

    def __init__(self, local_csv: str | None = None, fixture_path: str | None = None, **kwargs: Any):
        self.local_csv = local_csv
        self.fixture = load_fixture(fixture_path)
        self.config = kwargs

    def enrich(self, molecule_artifact):
        data = molecule_artifact.data.copy()
        fixture_hit = lookup_fixture(self.name, data, self.fixture)
        if fixture_hit is not None:
            result = dict(fixture_hit)
            result.setdefault("provider_status", "success")
            result.setdefault("source", "offline_fixture")
            data.setdefault("public_data", {})["atct"] = result
            refs = data.setdefault("public_data", {}).setdefault("reference_values", {})
            for key in ("delta_f_H_298_kj_mol", "S_298_J_molK", "Cp_298_J_molK"):
                if result.get(key) is not None:
                    refs[key] = {"value": result.get(key), "source": result.get("source", "offline_fixture")}
            data.setdefault("db_provider_status", {})[self.name] = {"status": "success", "matched": True}
            return data
        out = {"matched": False, "provider_status": "not_configured", "source": "ATcT/local curated table"}
        try:
            if self.local_csv and Path(self.local_csv).exists():
                df = pd.read_csv(self.local_csv)
                masks = []
                for col, value in [("inchikey", data.get("inchikey")), ("formula", data.get("formula")), ("name", data.get("name"))]:
                    if col in df.columns and value:
                        masks.append(df[col].astype(str).str.lower() == str(value).lower())
                if masks:
                    mask = masks[0]
                    for m in masks[1:]:
                        mask = mask | m
                    hit = df[mask].head(1)
                    if not hit.empty:
                        row = hit.iloc[0].to_dict()
                        out = {"matched": True, "provider_status": "success", "record": row, "source": "ATcT/local curated table"}
                        refs = data.setdefault("public_data", {}).setdefault("reference_values", {})
                        for key in ("delta_f_H_298_kj_mol", "S_298_J_molK", "Cp_298_J_molK"):
                            if key in row and pd.notna(row[key]):
                                refs[key] = {"value": float(row[key]), "source": "ATcT/local curated table"}
        except Exception as exc:
            out = {"matched": False, "provider_status": "failed", "reason": str(exc), "source": "ATcT/local curated table"}
        data.setdefault("public_data", {})["atct"] = out
        data.setdefault("db_provider_status", {})[self.name] = {"status": out.get("provider_status"), "matched": out.get("matched", False)}
        return data

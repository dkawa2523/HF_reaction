from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.backends.db.common import load_fixture, lookup_fixture

from hfauto.core.reference_data import load_reference_set, match_reference


class CCCBDBProvider:
    """NIST CCCBDB adapter for small-molecule method validation references."""

    name = "cccbdb"

    def __init__(self, table_path: str | None = None, reference_set_path: str | None = None, **_kwargs: Any):
        self.table_path = Path(table_path) if table_path else None
        self.reference_set = load_reference_set(reference_set_path)

    def _lookup_table(self, data: dict[str, Any]) -> dict[str, Any] | None:
        if not self.table_path or not self.table_path.exists():
            return None
        try:
            df = pd.read_csv(self.table_path)
        except Exception:
            return None
        inchikey = data.get("inchikey")
        name = str(data.get("name") or "").lower()
        if inchikey and "inchikey" in df.columns:
            hit = df[df["inchikey"] == inchikey]
        elif name and "name" in df.columns:
            hit = df[df["name"].astype(str).str.lower() == name]
        else:
            hit = pd.DataFrame()
        if hit.empty:
            return None
        return {"matched": True, "provider_status": "local_table_match", "record": hit.iloc[0].to_dict(), "table_path": str(self.table_path)}

    def enrich(self, molecule_artifact):
        data = molecule_artifact.data.copy()
        data.setdefault("public_data", {})
        result = self._lookup_table(data)
        if result is None:
            ref = match_reference(data, self.reference_set)
            if ref and "cccbdb" in (ref.get("sources") or []):
                result = {
                    "matched": True,
                    "provider_status": "builtin_reference_match",
                    "reference_id": ref.get("ref_id"),
                    "match_type": ref.get("match_type"),
                    "match_confidence": ref.get("match_confidence"),
                    "properties": {k: v for k, v in (ref.get("properties") or {}).items() if k in {"hf_stretch_cm1", "dipole_D"}},
                    "source": "bundled_hf_amine_reference_set",
                }
        if result is None:
            result = {"matched": False, "provider_status": "not_executed", "intended_use": "small-molecule geometry, frequency, dipole validation"}
        data["public_data"]["cccbdb"] = result
        data.setdefault("db_provider_status", {})[self.name] = {"status": result.get("provider_status"), "matched": bool(result.get("matched"))}
        return data

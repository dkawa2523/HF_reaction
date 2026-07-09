from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.core.reference_data import normalize_key


class CASCommonChemistryProvider:
    """Optional CAS Common Chemistry cross-check adapter.

    CAS Common Chemistry public site terms/API access should be reviewed before
    production automation.  This adapter therefore supports project-provided local
    crosswalk tables and otherwise records a non-blocking skipped status.
    """

    name = "cas_common_chemistry"

    def __init__(self, crosswalk_path: str | None = None, **_kwargs: Any):
        self.crosswalk_path = Path(crosswalk_path) if crosswalk_path else None

    def enrich(self, molecule_artifact):
        data = molecule_artifact.data.copy()
        data.setdefault("public_data", {})
        result = {"matched": False, "provider_status": "not_executed", "note": "Provide a licensed/local CAS crosswalk for automated use."}
        if self.crosswalk_path and self.crosswalk_path.exists():
            try:
                df = pd.read_csv(self.crosswalk_path)
                name_key = normalize_key(data.get("name"))
                hit = df[df.get("name", pd.Series(dtype=str)).astype(str).map(normalize_key) == name_key] if "name" in df.columns else pd.DataFrame()
                if data.get("inchikey") and "inchikey" in df.columns:
                    by_inchi = df[df["inchikey"] == data.get("inchikey")]
                    if not by_inchi.empty:
                        hit = by_inchi
                if not hit.empty:
                    result = {"matched": True, "provider_status": "local_crosswalk_match", "record": hit.iloc[0].to_dict(), "crosswalk_path": str(self.crosswalk_path)}
            except Exception as exc:
                result = {"matched": False, "provider_status": "failed", "reason": str(exc)}
        data["public_data"]["cas_common_chemistry"] = result
        data.setdefault("db_provider_status", {})[self.name] = {"status": result.get("provider_status"), "matched": bool(result.get("matched"))}
        return data

from __future__ import annotations

"""Run-level calculation cache and duplicate detection."""

from pathlib import Path
from typing import Any

import pandas as pd

from hfauto.core.hashing import fingerprint_dict, sha256_file
from hfauto.core.io import ensure_dir, write_json
from hfauto.core.schemas.manifest import Manifest


def artifact_cache_key(artifact: Any) -> str:
    paths = getattr(artifact, "paths", {}) or {}
    method = getattr(artifact, "method", None) or {}
    data = getattr(artifact, "data", {}) or {}
    path_hashes: dict[str, str] = {}
    for key, value in paths.items():
        try:
            p = Path(value)
            if p.exists() and p.is_file() and p.stat().st_size < 50_000_000:
                path_hashes[key] = sha256_file(p)
        except Exception:
            continue
    payload = {
        "artifact_type": getattr(artifact, "artifact_type", None),
        "artifact_id": getattr(artifact, "artifact_id", None),
        "species_id": data.get("species_id"),
        "reaction_id": data.get("reaction_id"),
        "state": data.get("state"),
        "hf_n": data.get("hf_n"),
        "method": method,
        "path_hashes": path_hashes,
    }
    return fingerprint_dict(payload)


def build_cache_index(manifest: Manifest, run_dir: str | Path | None = None) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for art in manifest.artifacts:
        if art.artifact_type not in {"calculation", "species", "species_preopt", "species_optimized", "thermo", "kinetics"}:
            continue
        key = artifact_cache_key(art)
        rows.append({
            "cache_key": key,
            "run_id": manifest.run_id,
            "artifact_id": art.artifact_id,
            "artifact_type": art.artifact_type,
            "status": art.status.status,
            "species_id": art.data.get("species_id"),
            "reaction_id": art.data.get("reaction_id"),
            "state": art.data.get("state"),
            "hf_n": art.data.get("hf_n"),
            "engine": (art.method or {}).get("engine"),
            "method_id": (art.method or {}).get("method_id"),
            "task": (art.method or {}).get("task"),
            "stage": (art.method or {}).get("stage"),
            "paths": art.paths,
            "run_dir": str(run_dir) if run_dir else None,
        })
    return pd.DataFrame(rows)


def write_cache_index(manifest: Manifest, out_dir: str | Path, run_dir: str | Path | None = None) -> dict[str, Path]:
    out = ensure_dir(out_dir)
    df = build_cache_index(manifest, run_dir=run_dir)
    csv_path = out / "calculation_cache_index.csv"
    json_path = out / "calculation_cache_index.json"
    df.to_csv(csv_path, index=False)
    write_json(json_path, {"schema_version": "hfauto.cache_index.v1", "n_rows": int(len(df)), "rows": df.to_dict(orient="records")})
    return {"csv": csv_path, "json": json_path}


def duplicate_cache_keys(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty or "cache_key" not in df.columns:
        return pd.DataFrame()
    counts = df.groupby("cache_key").size().reset_index(name="n")
    dup_keys = set(counts[counts["n"] > 1]["cache_key"])
    return df[df["cache_key"].isin(dup_keys)].sort_values(["cache_key", "artifact_id"])

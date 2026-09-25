from __future__ import annotations

from pathlib import Path

import pandas as pd


def build_status_snapshot(plan_csv: str | Path | None, scheduler: str = "slurm") -> pd.DataFrame:
    if plan_csv and Path(plan_csv).exists():
        try:
            plan = pd.read_csv(plan_csv)
        except pd.errors.EmptyDataError:
            plan = pd.DataFrame()
    else:
        plan = pd.DataFrame()
    rows = []
    if plan.empty:
        return pd.DataFrame(columns=["job_id", "stage", "scheduler", "scheduler_state", "hfauto_state", "note"])
    for _, row in plan.iterrows():
        rows.append(
            {
                "job_id": row.get("job_id"),
                "stage": row.get("stage"),
                "scheduler": scheduler,
                "scheduler_state": "not_submitted",
                "hfauto_state": row.get("status", "planned"),
                "target_artifact_id": row.get("target_artifact_id"),
                "command": row.get("command"),
                "note": "Offline snapshot. Use generated status command on production scheduler for live state.",
            }
        )
    return pd.DataFrame(rows)


def write_status_snapshot(plan_csv: str | Path | None, out_dir: str | Path, scheduler: str = "slurm") -> dict[str, Path]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    df = build_status_snapshot(plan_csv, scheduler)
    csv = out / "job_status_snapshot.csv"
    json = out / "job_status_snapshot.json"
    df.to_csv(csv, index=False)
    json.write_text(df.to_json(orient="records", indent=2), encoding="utf-8")
    return {"status_csv": csv, "status_json": json}

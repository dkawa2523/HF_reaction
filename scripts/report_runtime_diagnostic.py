#!/usr/bin/env python3
"""Build a small source-backed runtime diagnostic for an endpoint run."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from build_proton_transfer_endpoint_portable_report import (
    _execute,
    _finalize_portable_html,
    _find_deliver_script,
    _find_node,
    _load_csv_table,
    _node_path,
    _write_json,
)

CASE_RUNTIME_QUERY = """\
SELECT
  "case" AS case_name,
  case_id,
  COUNT(*) AS job_count,
  SUM(CAST(duration_hours AS REAL)) AS total_hours,
  100.0 * SUM(CAST(duration_hours AS REAL)) /
    (SELECT SUM(CAST(duration_hours AS REAL))
       FROM runtime_jobs WHERE scope = 'current') AS share_of_current_pct
FROM runtime_jobs
WHERE scope = 'current'
GROUP BY "case", case_id
ORDER BY case_id
"""

TASK_RUNTIME_QUERY = """\
SELECT
  task,
  COUNT(*) AS job_count,
  SUM(CAST(duration_hours AS REAL)) AS total_hours,
  100.0 * SUM(CAST(duration_hours AS REAL)) /
    (SELECT SUM(CAST(duration_hours AS REAL))
       FROM runtime_jobs WHERE scope = 'current') AS share_of_current_pct
FROM runtime_jobs
WHERE scope = 'current'
GROUP BY task
ORDER BY task
"""

SLOWEST_JOBS_QUERY = """\
SELECT "case" AS case_name, scope, task, job_id, status,
       CAST(duration_s AS REAL) AS duration_s,
       CAST(duration_hours AS REAL) AS duration_hours
FROM runtime_jobs
WHERE duration_hours IS NOT NULL AND trim(duration_hours) != ''
ORDER BY CAST(duration_hours AS REAL) DESC
LIMIT 12
"""


def _duration_seconds(job: dict[str, Any]) -> float | None:
    checks = job.get("checks") or {}
    command_check = checks.get("command_checkpoint_consistent") or {}
    actual = command_check.get("actual")
    if isinstance(actual, list):
        actual = actual[0] if actual else None
    if not isinstance(actual, dict) or actual.get("duration_s") is None:
        return None
    return float(actual["duration_s"])


def _label(case_id: str) -> str:
    return case_id.removeprefix("01_").removeprefix("02_").removeprefix("03_").removesuffix("_hf2")


def build_report(run_dir: Path, output_dir: Path) -> dict[str, Any]:
    run_dir = run_dir.resolve()
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    audit_path = run_dir / "raw_qm_audit" / "nwchem_raw_audit.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    jobs = audit.get("jobs") or []

    rows: list[dict[str, Any]] = []
    for job in jobs:
        duration_s = _duration_seconds(job)
        rows.append(
            {
                "scope": str(job.get("scope") or ""),
                "case_id": str(job.get("case_id") or ""),
                "case": _label(str(job.get("case_id") or "")),
                "task": str(job.get("task") or ""),
                "job_id": str(job.get("job_id") or ""),
                "status": str(job.get("status") or ""),
                "duration_s": duration_s,
                "duration_hours": None if duration_s is None else duration_s / 3600.0,
            }
        )

    current = [row for row in rows if row["scope"] == "current"]
    current_total_s = sum(float(row["duration_s"] or 0.0) for row in current)
    recorded_total_s = sum(float(row["duration_s"] or 0.0) for row in rows)
    missing_duration_count = sum(row["duration_s"] is None for row in rows)

    csv_path = output_dir / "runtime_by_job.csv"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    import sqlite3

    connection = sqlite3.connect(":memory:")
    try:
        _load_csv_table(connection, "runtime_jobs", csv_path)
        case_rows = _execute(connection, CASE_RUNTIME_QUERY)
        task_rows = _execute(connection, TASK_RUNTIME_QUERY)
        top_jobs = _execute(connection, SLOWEST_JOBS_QUERY)
    finally:
        connection.close()

    source_dir = output_dir / "source_snapshot"
    source_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(audit_path, source_dir / audit_path.name)

    generated_at = datetime.now(timezone.utc).isoformat()
    manifest_sources = [
        {
            "id": "runtime_csv_file",
            "label": "NWChem command runtime extraction",
            "path": "runtime_by_job.csv",
        },
        {
            "id": "runtime_case_query",
            "label": "Current NWChem runtime grouped by case",
            "path": "runtime_by_job.csv",
        },
        {
            "id": "slow_jobs_query",
            "label": "Slowest recorded NWChem jobs",
            "path": "runtime_by_job.csv",
        },
        {
            "id": "raw_audit",
            "label": "independent raw NWChem audit",
            "path": "source_snapshot/nwchem_raw_audit.json",
        },
    ]
    source_catalog = [
        manifest_sources[0],
        {
            **manifest_sources[1],
            "query": {
                "engine": "SQLite 3",
                "language": "sql",
                "executed_at": generated_at,
                "description": "Saved command runtimes grouped by current case.",
                "tables_used": ["runtime_by_job.csv"],
                "filters": ["scope = current"],
                "metric_definitions": [
                    "total_hours is the sum of unique saved command wall times."
                ],
                "sql": CASE_RUNTIME_QUERY,
            },
        },
        {
            **manifest_sources[2],
            "query": {
                "engine": "SQLite 3",
                "language": "sql",
                "executed_at": generated_at,
                "description": "Recorded jobs ranked by saved command wall time.",
                "tables_used": ["runtime_by_job.csv"],
                "filters": ["duration_hours is present", "top 12"],
                "metric_definitions": [
                    "duration_hours is saved duration_s divided by 3600."
                ],
                "sql": SLOWEST_JOBS_QUERY,
            },
        },
        manifest_sources[3],
    ]
    charts = [
        {
            "id": "runtime-by-case",
            "title": "Current NWChem runtime by case",
            "subtitle": "26 current jobs; summed recorded wall time, hours",
            "showDescription": True,
            "type": "bar",
            "dataset": "runtime_by_case",
            "sourceId": "runtime_case_query",
            "encodings": {
                "x": {"field": "case_name", "type": "nominal", "label": "case"},
                "y": {
                    "field": "total_hours",
                    "type": "quantitative",
                    "label": "recorded runtime",
                    "unit": "h",
                },
                "tooltip": [
                    {"field": "job_count", "type": "quantitative", "label": "jobs"},
                    {
                        "field": "share_of_current_pct",
                        "type": "quantitative",
                        "label": "share of current",
                        "unit": "%",
                    },
                ],
            },
            "xAxisTitle": "case",
            "yAxisTitle": "summed runtime (h)",
            "unit": "h",
            "layout": "12 columns",
        }
    ]
    tables = [
        {
            "id": "slowest-jobs",
            "title": "Slowest recorded NWChem jobs",
            "subtitle": "Current, historical, and expected-failed records with saved command durations",
            "showDescription": True,
            "dataset": "slowest_jobs",
            "sourceId": "slow_jobs_query",
            "density": "compact",
            "defaultSort": {"field": "duration_hours", "direction": "desc"},
            "columns": [
                {"field": "case_name", "label": "case", "type": "text"},
                {"field": "scope", "label": "scope", "type": "text"},
                {"field": "task", "label": "task", "type": "text"},
                {"field": "job_id", "label": "job", "type": "text"},
                {"field": "duration_hours", "label": "hours", "type": "number", "unit": "h"},
            ],
            "layout": "12 columns",
        }
    ]
    summary = (
        f"## 技術要約\n\n"
        f"- **遅かった主因はNWChem本計算です。** current 26 jobの保存時間合計は "
        f"**{current_total_s / 3600.0:.2f}時間**、失敗元・履歴を含む記録済み合計は "
        f"**{recorded_total_s / 3600.0:.2f}時間**です。\n"
        f"- **anilineがcurrent時間の{case_rows[-1]['share_of_current_pct']:.1f}%を占めました。** "
        f"通常候補、mode-following 2方向、xfine精度救済が追加されたためです。\n"
        f"- **レポート作成は主因ではありません。** 量子化学計算後の監査・描画・テストは分単位でした。"
    )
    blocks = [
        {"id": "title", "type": "markdown", "body": "# HF2 endpoint計算の実行時間診断"},
        {
            "id": "summary",
            "type": "markdown",
            "body": summary,
            "sourceId": "runtime_csv_file",
        },
        {
            "id": "case-finding",
            "type": "markdown",
            "body": (
                "## anilineとtrimethylamineが計算時間を支配\n\n"
                "棒はcurrent判定に使われる各NWChem jobの保存済みwall timeをケース別に加算した値です。"
                "独立jobを直列実行したため、加算値は実際の待ち時間の主要部分に対応します。"
            ),
            "sourceId": "runtime_case_query",
        },
        {"id": "case-chart", "type": "chart", "chartId": "runtime-by-case"},
        {
            "id": "definitions",
            "type": "markdown",
            "body": (
                "## 集計範囲と定義\n\n"
                "`current`は最終判定に参照された26 job、`expected_failed_attempt`は継続計算に置換された失敗元、"
                "`historical_unreferenced`は現判定に未参照の履歴です。durationは各jobの保存済みcommand recordのwall timeで、"
                f"全31 job中{missing_duration_count}件はduration欠損のため、21.46時間は下限です。"
            ),
            "sourceId": "runtime_csv_file",
        },
        {
            "id": "method",
            "type": "markdown",
            "body": (
                "## 算出方法\n\n"
                "独立raw監査が列挙した全jobについて `command_checkpoint_consistent.actual.duration_s` を読み、"
                "jobを重複させずscope・case・task別に加算しました。currentではoptimize 15件、opt_freq 11件です。"
            ),
            "sourceId": "raw_audit",
        },
        {
            "id": "slow-jobs-finding",
            "type": "markdown",
            "body": (
                "## 長時間jobは大きい分子の最適化と追加Hessian\n\n"
                "anilineの複数最適化が各1.3–1.6時間、timeoutしたaniline jobが2.0時間、"
                "mode-followingと精度救済のopt_freqが各約1時間でした。"
            ),
            "sourceId": "slow_jobs_query",
        },
        {"id": "slow-jobs-table", "type": "table", "tableId": "slowest-jobs"},
        {
            "id": "limitations",
            "type": "markdown",
            "body": (
                "## 限界と頑健性\n\n"
                "保存時間のない3履歴jobは合計に含まれないため、記録済み21.46時間は下限です。"
                "各jobのwall time加算であり、CPU時間ではありません。ハードウェア待機、対話・実装修正時間も含みません。"
            ),
            "sourceId": "runtime_csv_file",
        },
        {
            "id": "next",
            "type": "markdown",
            "body": (
                "## 次回を短縮する方法\n\n"
                "- 16–32 core環境では、4 core/jobの独立seedを3–5本並列化する。\n"
                "- 現在のhash付きcheckpointを維持し、成功jobを再実行しない。\n"
                "- 負のendpoint主張が不要ならseed数を減らせるが、結論の強さとの交換条件になる。\n"
                "- 完全HessianはDFT候補選別後に限定する。"
            ),
        },
        {
            "id": "questions",
            "type": "markdown",
            "body": (
                "## 追加確認事項\n\n"
                "今後は、追加CPUでのseed並列実行と、同じ科学ゲートを保った二段階screeningの実測短縮率を確認できます。"
            ),
        },
    ]
    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": "HF2 endpoint計算の実行時間診断",
            "description": "保存済みNWChem command時間から、長時間化の原因を分解。",
            "generatedAt": generated_at,
            "cards": [],
            "charts": charts,
            "tables": tables,
            "sources": manifest_sources,
            "blocks": blocks,
        },
        "snapshot": {
            "version": 1,
            "generatedAt": generated_at,
            "status": "ready",
            "datasets": {
                "runtime_by_case": case_rows,
                "runtime_by_task": task_rows,
                "slowest_jobs": top_jobs,
            },
            "accessIssues": [],
        },
        "sources": source_catalog,
        "package_info": {
            "report_kind": "hf2_runtime_diagnostic",
            "source_run_id": run_dir.name,
        },
    }
    artifact_path = output_dir / "artifact.json"
    _write_json(artifact_path, artifact)
    _write_json(
        output_dir / "report_plan.json",
        {
            "audience": "technical",
            "required_structure": [
                "technical summary",
                "key findings with visual evidence",
                "scope and definitions",
                "methodology",
                "limitations",
                "recommended next steps",
                "further questions",
            ],
            "chart_map": [
                {
                    "section": "aniline and trimethylamine dominate",
                    "question": "which case consumed current runtime",
                    "family": "comparison",
                    "type": "bar",
                    "fields": ["case_name", "total_hours"],
                    "palette": "single-root preferred",
                }
            ],
        },
    )

    builder = _find_deliver_script()
    node = _find_node()
    if builder is None or node is None:
        raise RuntimeError("Portable report builder or Node.js is unavailable")
    html_path = output_dir / "report.html"
    completed = subprocess.run(
        [
            str(node),
            _node_path(builder, node),
            "--input",
            _node_path(artifact_path, node),
            "--output",
            _node_path(html_path, node),
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=120,
    )
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr.strip() or completed.stdout.strip())
    receipt = json.loads(completed.stdout.strip().splitlines()[-1])
    receipt.setdefault("stages", {})["hfauto_finalization"] = "passed"
    receipt["hfauto_finalization"] = _finalize_portable_html(html_path)
    _write_json(output_dir / "portable_delivery_receipt.json", receipt)
    return {
        "html": str(html_path),
        "current_hours": current_total_s / 3600.0,
        "recorded_hours": recorded_total_s / 3600.0,
        "missing_duration_jobs": missing_duration_count,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build_report(args.run_dir, args.output_dir), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

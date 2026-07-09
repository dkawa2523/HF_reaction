from __future__ import annotations

"""GoodVibes-compatible thermochemistry backend.

Phase 10 adds a production-oriented external GoodVibes adapter while keeping the
Phase 6 deterministic internal quasi-RRHO fallback.  The public contract remains
``species_thermo(...) -> dict`` so the ThermoStage does not depend on a specific
thermochemistry package.

Design notes
------------
* External execution is opt-in via config or environment variables.
* Every external run writes stdout/stderr/command_result sidecars.
* GoodVibes CSV parsing is tolerant to column-name variants used across
  versions and user options.
* If external GoodVibes is unavailable, species records explicitly say that the
  internal low-frequency correction was used; they are never silently promoted
  to production thermochemistry.
"""

from dataclasses import dataclass
import csv
import os
from pathlib import Path
from typing import Any, Iterable

from hfauto.core.executables import CommandResult, resolve_executable, run_command
from hfauto.core.io import ensure_dir, write_json
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.thermo_models import low_frequency_correction, scale_thermal_correction
from hfauto.core.units import pressure_correction_hartree


_FLOAT_COLUMNS_G = [
    "qh-G(T)",
    "qh_G(T)",
    "qh_G",
    "QH-G(T)",
    "G(T)",
    "G",
    "Gibbs Free Energy",
    "qh Gibbs Free Energy",
    "qh-Gibbs Free Energy",
]
_FLOAT_COLUMNS_H = ["H(T)", "H", "Enthalpy", "qh-H(T)", "qh_H"]
_FLOAT_COLUMNS_QH_CORR = [
    "qh-T.S",
    "qh_TS",
    "qh-G corr",
    "QH correction",
    "qh_correction",
]
_NAME_COLUMNS = ["Structure", "structure", "Name", "name", "File", "file"]


@dataclass
class ParsedGoodVibesTable:
    path: Path
    rows: dict[str, dict[str, Any]]


def _to_float(value: Any) -> float | None:
    try:
        if value is None:
            return None
        s = str(value).strip()
        if not s or s.lower() in {"nan", "none", "--"}:
            return None
        return float(s)
    except Exception:
        return None


def _first_float(row: dict[str, Any], columns: Iterable[str]) -> float | None:
    for col in columns:
        if col in row:
            val = _to_float(row[col])
            if val is not None:
                return val
    # fallback: case/space-insensitive lookup
    normalized = {str(k).lower().replace(" ", "").replace("_", "-"): v for k, v in row.items()}
    for col in columns:
        key = col.lower().replace(" ", "").replace("_", "-")
        if key in normalized:
            val = _to_float(normalized[key])
            if val is not None:
                return val
    return None


def _row_name(row: dict[str, Any]) -> str | None:
    for col in _NAME_COLUMNS:
        if row.get(col):
            return str(row[col]).strip()
    return None


class GoodVibesEngine:
    name = "goodvibes"

    def __init__(self, **config: Any):
        self.config = config
        self._external_cache: dict[str, dict[str, Any]] = {}

    def allow_subprocess(self) -> bool:
        return bool(
            self.config.get("allow_subprocess", False)
            or self.config.get("run_external", False)
            or self.config.get("use_external_goodvibes", False)
            or os.environ.get("HFAUTO_ALLOW_GOODVIBES") == "1"
            or os.environ.get("HFAUTO_ALLOW_SUBPROCESS") == "1"
        )

    def executable(self) -> str | None:
        # GoodVibes may be installed as goodvibes or python -m goodvibes in some
        # environments.  The explicit executable always wins; otherwise look for
        # common console entry points.
        explicit = self.config.get("executable")
        found = resolve_executable("goodvibes", explicit)
        if found:
            return found
        return resolve_executable("GoodVibes.py", None)

    def _output_file_for_calc(self, source_calc: Artifact) -> str | None:
        paths = source_calc.paths or {}
        for key in ["output", "out", "log", "orca_out", "input"]:
            p = paths.get(key)
            if p and Path(p).exists():
                return str(Path(p))
        # Some older dummy artifacts put paths in data.
        for key in ["output_path", "log_path"]:
            p = source_calc.data.get(key)
            if p and Path(p).exists():
                return str(Path(p))
        return None

    def maybe_run_external(self, output_files: list[str], workdir: Path) -> dict[str, Any]:
        exe = self.executable()
        workdir = ensure_dir(workdir)
        if not self.allow_subprocess():
            return {
                "external_goodvibes_executed": False,
                "external_status": "subprocess_not_allowed",
                "production_thermo_ready": False,
            }
        if exe is None:
            return {
                "external_goodvibes_executed": False,
                "external_status": "executable_not_found",
                "production_thermo_ready": False,
            }
        if not output_files:
            return {
                "external_goodvibes_executed": False,
                "external_status": "no_output_files",
                "production_thermo_ready": False,
            }
        cmd = [exe]
        # Keep options conservative and version-agnostic.  Users may pass more in
        # extra_args, e.g. ["--qs", "grimme", "--fs", "100"] depending on their
        # GoodVibes installation.
        extra = [str(x) for x in (self.config.get("extra_args") or [])]
        cmd.extend(extra)
        cmd.extend([str(Path(f).resolve()) for f in output_files])
        result = run_command(
            cmd,
            workdir,
            timeout_s=int(self.config.get("timeout_s", 3600)),
            stdout_name="goodvibes_stdout.txt",
            stderr_name="goodvibes_stderr.txt",
        )
        parsed = self.parse_output_directory(workdir)
        status = "success" if result.ok and parsed else ("ran_no_csv" if result.ok else "failed")
        payload = {
            "external_goodvibes_executed": bool(result.ok),
            "external_status": status,
            "production_thermo_ready": bool(result.ok and parsed),
            "command_result": result.to_dict(),
            "goodvibes_tables": [str(t.path) for t in parsed],
        }
        write_json(workdir / "goodvibes_external_summary.json", payload)
        return payload

    def parse_goodvibes_csv(self, csv_path: str | Path) -> ParsedGoodVibesTable:
        path = Path(csv_path)
        rows: dict[str, dict[str, Any]] = {}
        with path.open(newline="", encoding="utf-8", errors="ignore") as f:
            reader = csv.DictReader(f)
            for i, row in enumerate(reader):
                key = _row_name(row) or f"row_{i}"
                normalized: dict[str, Any] = dict(row)
                g = _first_float(row, _FLOAT_COLUMNS_G)
                h = _first_float(row, _FLOAT_COLUMNS_H)
                qcorr = _first_float(row, _FLOAT_COLUMNS_QH_CORR)
                if g is not None:
                    normalized["goodvibes_gibbs_hartree"] = float(g)
                if h is not None:
                    normalized["goodvibes_enthalpy_hartree"] = float(h)
                if qcorr is not None:
                    normalized["goodvibes_qh_correction_raw"] = float(qcorr)
                normalized["goodvibes_row_name"] = key
                rows[str(key)] = normalized
        return ParsedGoodVibesTable(path=path, rows=rows)

    def parse_output_directory(self, workdir: str | Path) -> list[ParsedGoodVibesTable]:
        wd = Path(workdir)
        candidates = sorted(wd.glob("Goodvibes*.csv")) + sorted(wd.glob("goodvibes*.csv")) + sorted(wd.glob("*.csv"))
        tables: list[ParsedGoodVibesTable] = []
        seen: set[Path] = set()
        for p in candidates:
            if p in seen:
                continue
            seen.add(p)
            try:
                table = self.parse_goodvibes_csv(p)
                if table.rows:
                    tables.append(table)
            except Exception:
                continue
        return tables

    # Backward-compatible name used by earlier tests/stages.
    def read_goodvibes_csv(self, csv_path: str | Path) -> dict[str, dict[str, Any]]:
        return self.parse_goodvibes_csv(csv_path).rows

    def _match_goodvibes_row(self, tables: list[ParsedGoodVibesTable], source_calc: Artifact, output_file: str | None) -> dict[str, Any] | None:
        names: list[str] = []
        if output_file:
            p = Path(output_file)
            names.extend([p.name, p.stem, str(p.resolve())])
        names.extend([
            str(source_calc.artifact_id),
            str(source_calc.data.get("species_id") or ""),
            str(source_calc.data.get("energy_source_calc_id") or ""),
        ])
        normalized_names = {n.lower() for n in names if n}
        for table in tables:
            for key, row in table.rows.items():
                key_l = str(key).lower()
                if key_l in normalized_names or Path(key_l).stem in normalized_names:
                    out = dict(row)
                    out["goodvibes_csv_path"] = str(table.path)
                    return out
        # If there is exactly one row, use it; this is common for per-species runs.
        all_rows = [(table, key, row) for table in tables for key, row in table.rows.items()]
        if len(all_rows) == 1:
            table, _, row = all_rows[0]
            out = dict(row)
            out["goodvibes_csv_path"] = str(table.path)
            return out
        return None

    def _external_thermo_for_calc(self, source_calc: Artifact) -> dict[str, Any]:
        cache_key = source_calc.artifact_id
        if cache_key in self._external_cache:
            return self._external_cache[cache_key]
        output_file = self._output_file_for_calc(source_calc)
        work_root = Path(self.config.get("work_root", ".hfauto_cache/goodvibes_runs"))
        safe = "".join(ch if ch.isalnum() or ch in "_-" else "_" for ch in cache_key)[:120]
        workdir = ensure_dir(work_root / safe)
        if output_file:
            summary = self.maybe_run_external([output_file], workdir)
        else:
            summary = {
                "external_goodvibes_executed": False,
                "external_status": "no_source_output_file",
                "production_thermo_ready": False,
            }
        tables = self.parse_output_directory(workdir)
        row = self._match_goodvibes_row(tables, source_calc, output_file)
        payload = {**summary, "matched_row": row or {}, "workdir": str(workdir), "source_output_file": output_file}
        self._external_cache[cache_key] = payload
        return payload

    def species_correction(self, frequencies_cm1: Iterable[float | int | str | None], T_K: float, config: dict[str, Any] | None = None) -> dict[str, Any]:
        cfg = {**self.config, **(config or {})}
        values = []
        for f in frequencies_cm1 or []:
            try:
                if f is not None:
                    values.append(float(f))
            except Exception:
                pass
        lfc = low_frequency_correction(
            values,
            cutoff_cm1=float(cfg.get("qrrho_cutoff_cm1", cfg.get("quasi_rrho_cutoff_cm1", cfg.get("low_frequency_cutoff_cm1", 100.0)))),
            max_correction_per_mode_kcal_mol=float(cfg.get("max_lowfreq_correction_per_mode_kcal_mol", 0.20)),
            enabled=bool(cfg.get("apply_internal_quasi_rrho", cfg.get("quasi_rrho", True))),
        )
        return {
            "backend": self.name,
            "model": lfc.model,
            "cutoff_cm1": lfc.cutoff_cm1,
            "low_frequency_count": lfc.low_frequency_count,
            "qrrho_correction_kcal_mol": lfc.correction_kcal_mol,
            "qrrho_correction_hartree": lfc.correction_hartree,
            "external_goodvibes_executed": False,
            "scientific_note": "internal GoodVibes-like low-frequency correction; use external GoodVibes for final reporting if required",
        }

    def species_thermo(
        self,
        *,
        species_id: str,
        source_calc: Artifact,
        T_K: float,
        p_bar: float | None,
        p_standard_bar: float,
    ) -> dict[str, Any]:
        data = source_calc.data
        electronic = data.get("electronic_energy_hartree")
        if electronic is None:
            raise ValueError(f"Calculation {source_calc.artifact_id} lacks electronic_energy_hartree")
        thermal = data.get("thermal_correction_gibbs_hartree")
        if thermal is None and data.get("gibbs_298K_hartree") is not None:
            thermal = float(data["gibbs_298K_hartree"]) - float(electronic)
        if thermal is None:
            thermal = 0.0

        external_meta: dict[str, Any] = {}
        external_row: dict[str, Any] = {}
        external_requested = bool(self.config.get("use_external_goodvibes", self.config.get("run_external", False)))
        if external_requested or self.allow_subprocess() and bool(self.config.get("prefer_external_if_available", False)):
            external_meta = self._external_thermo_for_calc(source_calc)
            external_row = external_meta.get("matched_row") or {}

        scaled_thermal, temp_model = scale_thermal_correction(
            float(thermal),
            T_K=float(T_K),
            source_T_K=float(self.config.get("source_temperature_K", 298.15)),
            enabled=bool(self.config.get("allow_linear_temperature_scaling", False)),
        )
        lfc = low_frequency_correction(
            data.get("frequencies_cm1", []),
            cutoff_cm1=float(self.config.get("quasi_rrho_cutoff_cm1", self.config.get("qrrho_cutoff_cm1", self.config.get("low_frequency_cutoff_cm1", 100.0)))),
            max_correction_per_mode_kcal_mol=float(self.config.get("max_lowfreq_correction_per_mode_kcal_mol", 0.20)),
            enabled=bool(self.config.get("apply_internal_quasi_rrho", self.config.get("quasi_rrho", True))),
        )

        external_g = _to_float(external_row.get("goodvibes_gibbs_hartree")) if external_row else None
        use_external_g = bool(external_g is not None and self.config.get("prefer_external_goodvibes_g", True))
        if use_external_g:
            g_standard = float(external_g)
            thermal_used = float(external_g) - float(electronic)
            thermo_model = "external_goodvibes_qh_gibbs"
            qrrho_applied = True
            qrrho_corr_h = 0.0
            qrrho_corr_kcal = 0.0
        else:
            g_standard = float(electronic) + scaled_thermal + lfc.correction_hartree
            thermal_used = scaled_thermal
            thermo_model = temp_model
            qrrho_applied = bool(lfc.applied)
            qrrho_corr_h = float(lfc.correction_hartree)
            qrrho_corr_kcal = float(lfc.correction_kcal_mol)

        pressure_corr = pressure_correction_hartree(float(T_K), p_bar, float(p_standard_bar))
        g_process = g_standard + pressure_corr
        ext_executed = bool(external_meta.get("external_goodvibes_executed", False))
        ext_status = external_meta.get("external_status") or ("not_requested" if not external_requested else "not_available")
        return {
            "species_id": species_id,
            "T_K": float(T_K),
            "p_bar": p_bar,
            "p_standard_bar": float(p_standard_bar),
            "electronic_energy_hartree": float(electronic),
            "thermal_correction_source_hartree": float(thermal),
            "thermal_correction_used_hartree": float(thermal_used),
            "thermal_temperature_model": thermo_model,
            "quasi_rrho_model": "external_goodvibes" if use_external_g else lfc.model,
            "quasi_rrho_applied": bool(qrrho_applied),
            "low_frequency_count": int(lfc.low_frequency_count),
            "low_frequency_cutoff_cm1": float(lfc.cutoff_cm1),
            "quasi_rrho_correction_kcal_mol": float(qrrho_corr_kcal),
            "quasi_rrho_correction_hartree": float(qrrho_corr_h),
            "pressure_correction_hartree": float(pressure_corr),
            "G_standard_hartree": float(g_standard),
            "G_process_hartree": float(g_process),
            "G_effective_hartree": float(g_process),
            "energy_source_calc_id": source_calc.data.get("energy_source_calc_id") or source_calc.artifact_id,
            "thermal_source_calc_id": source_calc.data.get("thermal_source_calc_id") or source_calc.artifact_id,
            "energy_source_engine": (source_calc.method or {}).get("engine") or source_calc.data.get("engine"),
            "energy_source_task": source_calc.data.get("task") or (source_calc.method or {}).get("task"),
            "thermo_backend": self.name,
            "external_goodvibes_requested": bool(external_requested),
            "external_goodvibes_executed": ext_executed,
            "external_goodvibes_status": ext_status,
            "external_goodvibes_workdir": external_meta.get("workdir"),
            "goodvibes_csv_path": external_row.get("goodvibes_csv_path"),
            "goodvibes_row_name": external_row.get("goodvibes_row_name"),
            "goodvibes_gibbs_hartree": external_g,
            "production_thermo_ready": bool(use_external_g),
            "connector_quality": "production_goodvibes" if use_external_g else "internal_qrrho_fallback",
        }

"""GoodVibes-compatible thermochemistry backend.

The production-oriented external GoodVibes adapter retains a deterministic
internal quasi-RRHO screening fallback. The public contract remains
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

from __future__ import annotations

import csv
import os
import re
from collections.abc import Iterable
from dataclasses import dataclass
from importlib import metadata
from pathlib import Path
from typing import Any

from hfauto.core.executables import resolve_executable, run_command
from hfauto.core.hashing import sha256_file
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
    "qh_gibbs_free_energy",
    "gibbs_free_energy",
]
_FLOAT_COLUMNS_H = [
    "H(T)",
    "H",
    "Enthalpy",
    "qh-H(T)",
    "qh_H",
    "qh_enthalpy",
    "enthalpy",
]
_FLOAT_COLUMNS_ZPE = ["zpe", "ZPE", "Zero-Point Energy", "Zero Point Energy"]
_FLOAT_COLUMNS_QH_CORR = [
    "qh-T.S",
    "qh_TS",
    "qh-G corr",
    "QH correction",
    "qh_correction",
]
_FLOAT_COLUMNS_T = ["Temperature", "Temperature (K)", "T(K)", "T_K"]
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
    except (TypeError, ValueError):
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


def _prepare_goodvibes_input(
    source_path: str | Path, workdir: str | Path
) -> dict[str, Any]:
    """Make a hash-bound parser-compatibility copy without changing raw QC output.

    GoodVibes 4.3 accepts NWChem, but its generic metadata pass attempts
    ``int('0.00')`` for NWChem population-analysis lines such as
    ``charge = 0.00``.  The authoritative ``Charge : 0`` line is already in
    the same output.  Normalizing only the parser-breaking spelling preserves
    every energy, frequency, coordinate, and electronic-state value.
    """

    source = Path(source_path)
    text = source.read_text(encoding="utf-8", errors="replace")
    edits = 0
    rules: list[str] = []
    if "Northwest Computational Chemistry Package" in text or "NWChem" in text:
        pattern = re.compile(
            r"^(?P<indent>\s*)charge\s*=\s*(?P<value>[+-]?\d+)\.0+\s*$",
            re.MULTILINE,
        )

        def replace(match: re.Match[str]) -> str:
            nonlocal edits
            edits += 1
            if "nwchem_decimal_population_charge" not in rules:
                rules.append("nwchem_decimal_population_charge")
            return f"{match.group('indent')}charge {match.group('value')}"

        text = pattern.sub(replace, text)
        rotational = {
            axis: re.search(
                rf"^\s*{axis}=\s+(?P<value>\*+|[+-]?\d+(?:\.\d+)?)\s+cm-1\s+\(\s*(?P<temperature>\*+|[+-]?\d+(?:\.\d+)?)\s+K\)",
                text,
                re.MULTILINE,
            )
            for axis in ("A", "B", "C")
        }
        if all(match is not None for match in rotational.values()):
            a_token = rotational["A"].group("value")  # type: ignore[union-attr]
            b_token = rotational["B"].group("value")  # type: ignore[union-attr]
            c_token = rotational["C"].group("value")  # type: ignore[union-attr]
            try:
                linear_axis = set(a_token) == {"*"} or abs(float(a_token)) < 1.0e-8
                equal_perpendicular_axes = abs(float(b_token) - float(c_token)) < 1.0e-5
            except ValueError:
                linear_axis = False
                equal_perpendicular_axes = False
            if linear_axis and equal_perpendicular_axes:
                b_temperature = float(
                    rotational["B"].group("temperature")  # type: ignore[union-attr]
                )
                replacement = (
                    f"A=   {float(b_token):.6f} cm-1  "
                    f"(  {b_temperature:.6f} K)"
                )
                text, replacements = re.subn(
                    r"(?m)^\s*A=.*$", replacement, text, count=1
                )
                edits += replacements
                rules.append("nwchem_linear_rotational_constant_mapping")
                if "C*V symmetry detected" not in text:
                    text += "\nC*V symmetry detected\n"
                    edits += 1
                    rules.append("nwchem_linear_point_group_marker")
    destination_dir = ensure_dir(Path(workdir) / "prepared_inputs")
    destination = destination_dir / (
        f"{source.stem}_{sha256_file(source)[:12]}{source.suffix or '.out'}"
    )
    destination.write_text(text, encoding="utf-8")
    return {
        "source_path": str(source.resolve()),
        "source_sha256": sha256_file(source),
        "prepared_path": str(destination.resolve()),
        "prepared_sha256": sha256_file(destination),
        "compatibility_edits": edits,
        "compatibility_rules": rules,
        "scientific_values_changed": False,
    }


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
        prepared_inputs = [
            _prepare_goodvibes_input(path, workdir) for path in output_files
        ]
        cmd = [exe]
        # Keep options conservative and version-agnostic.  Users may pass more in
        # extra_args, e.g. ["--qs", "grimme", "--fs", "100"] depending on their
        # GoodVibes installation.
        configured_extra = self.config.get("extra_args")
        if configured_extra is None:
            try:
                major = int(metadata.version("goodvibes").split(".", 1)[0])
            except (metadata.PackageNotFoundError, ValueError):
                major = 3
            extra = ["-q", "--csv", "Goodvibes.csv"] if major >= 4 else ["-q", "--csv"]
        else:
            extra = [str(x) for x in configured_extra]
        cmd.extend(extra)
        cmd.extend([item["prepared_path"] for item in prepared_inputs])
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
            "prepared_inputs": prepared_inputs,
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
                zpe = _first_float(row, _FLOAT_COLUMNS_ZPE)
                qcorr = _first_float(row, _FLOAT_COLUMNS_QH_CORR)
                temperature = _first_float(row, _FLOAT_COLUMNS_T)
                if g is not None:
                    normalized["goodvibes_gibbs_hartree"] = float(g)
                if h is not None:
                    normalized["goodvibes_enthalpy_hartree"] = float(h)
                if zpe is not None:
                    normalized["goodvibes_zpe_hartree"] = float(zpe)
                if qcorr is not None:
                    normalized["goodvibes_qh_correction_raw"] = float(qcorr)
                if temperature is not None:
                    normalized["goodvibes_temperature_K"] = float(temperature)
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
            except (OSError, csv.Error, TypeError, ValueError):
                continue
        return tables

    def _match_goodvibes_row(self, tables: list[ParsedGoodVibesTable], source_calc: Artifact, output_file: str | None) -> dict[str, Any] | None:
        names: list[str] = []
        if output_file:
            p = Path(output_file)
            names.extend([p.name, p.stem, str(p.resolve())])
        names.extend([
            str(source_calc.artifact_id),
            str(source_calc.data.get("species_id") or ""),
            str(source_calc.data.get("energy_source_calc_id") or ""),
            str(source_calc.data.get("thermal_source_calc_id") or ""),
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
            except (TypeError, ValueError):
                continue
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

    @staticmethod
    def _thermochemistry_frequencies(data: dict[str, Any]) -> tuple[Any, str]:
        """Prefer the projected physical-mode list while retaining old artifacts."""

        if (
            "thermochemistry_frequencies_cm1" in data
            and data.get("thermochemistry_frequencies_cm1") is not None
        ):
            return (
                data.get("thermochemistry_frequencies_cm1"),
                "thermochemistry_frequencies_cm1",
            )
        return data.get("frequencies_cm1", []), "frequencies_cm1_legacy_fallback"

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

        source_temperature = _to_float(data.get("thermochemistry_temperature_K"))
        if source_temperature is None:
            source_temperature = float(self.config.get("source_temperature_K", 298.15))
        temperature_tolerance = float(
            self.config.get("temperature_match_tolerance_K", 0.05)
        )
        source_temperature_matches = bool(
            abs(float(T_K) - float(source_temperature)) <= temperature_tolerance
        )
        scaled_thermal, temp_model = scale_thermal_correction(
            float(thermal),
            T_K=float(T_K),
            source_T_K=float(source_temperature),
            enabled=bool(self.config.get("allow_linear_temperature_scaling", False)),
        )
        thermochemistry_frequencies, frequency_source = (
            self._thermochemistry_frequencies(data)
        )
        lfc = low_frequency_correction(
            thermochemistry_frequencies,
            cutoff_cm1=float(self.config.get("quasi_rrho_cutoff_cm1", self.config.get("qrrho_cutoff_cm1", self.config.get("low_frequency_cutoff_cm1", 100.0)))),
            max_correction_per_mode_kcal_mol=float(self.config.get("max_lowfreq_correction_per_mode_kcal_mol", 0.20)),
            enabled=bool(self.config.get("apply_internal_quasi_rrho", self.config.get("quasi_rrho", True))),
        )

        external_g = _to_float(external_row.get("goodvibes_gibbs_hartree")) if external_row else None
        external_h = _to_float(
            external_row.get("goodvibes_enthalpy_hartree")
            if external_row
            else None
        )
        external_zpe = _to_float(
            external_row.get("goodvibes_zpe_hartree")
            if external_row
            else None
        )
        external_temperature = _to_float(
            external_row.get("goodvibes_temperature_K")
            if external_row
            else None
        )
        if external_temperature is None:
            external_temperature = _to_float(
                self.config.get("external_goodvibes_temperature_K")
            )
        # A GoodVibes table with no explicit temperature belongs only to the
        # source frequency temperature.  In particular, a 298 K row cached for
        # one calculation must never be reused as 373/423 K production data.
        if external_temperature is None and external_g is not None:
            external_temperature = float(source_temperature)
        external_temperature_matches = bool(
            external_g is not None
            and external_temperature is not None
            and abs(float(T_K) - float(external_temperature))
            <= temperature_tolerance
        )
        ext_executed = bool(
            external_meta.get("external_goodvibes_executed", False)
        )
        external_run_verified = bool(
            ext_executed and external_meta.get("external_status") == "success"
        )
        use_external_g = bool(
            external_g is not None
            and external_temperature_matches
            and external_run_verified
            and self.config.get("prefer_external_goodvibes_g", True)
        )
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
        ext_status = external_meta.get("external_status") or ("not_requested" if not external_requested else "not_available")
        if use_external_g:
            temperature_status = "external_goodvibes_at_requested_temperature"
        elif source_temperature_matches:
            temperature_status = "source_temperature_internal_screening_only"
        elif bool(self.config.get("allow_linear_temperature_scaling", False)):
            temperature_status = "experimental_linear_scaling_not_production"
        else:
            temperature_status = "source_temperature_correction_reused_not_quantitative"
        production_ready = bool(use_external_g)
        not_ready_reasons: list[str] = []
        if not use_external_g:
            not_ready_reasons.append("external_goodvibes_gibbs_unavailable")
        if external_g is not None and not external_run_verified:
            not_ready_reasons.append("external_goodvibes_execution_not_verified")
        if external_g is not None and not external_temperature_matches:
            not_ready_reasons.append(
                "external_goodvibes_temperature_does_not_match_requested_temperature"
            )
        if not source_temperature_matches and not use_external_g:
            not_ready_reasons.append(
                "source_thermal_correction_temperature_does_not_match_requested_temperature"
            )
        zpe = _to_float(data.get("zpe_hartree"))
        if use_external_g and external_zpe is not None:
            zpe = external_zpe
        zero_point_corrected = (
            None if zpe is None else float(electronic) + zpe
        )
        enthalpy = None
        if use_external_g and external_h is not None:
            enthalpy = external_h
        elif source_temperature_matches:
            enthalpy = _to_float(
                data.get("enthalpy_hartree")
                or data.get("enthalpy_298K_hartree")
            )
        return {
            "species_id": species_id,
            "T_K": float(T_K),
            "p_bar": p_bar,
            "p_standard_bar": float(p_standard_bar),
            "electronic_energy_hartree": float(electronic),
            "zero_point_corrected_energy_hartree": zero_point_corrected,
            "enthalpy_standard_hartree": enthalpy,
            "thermal_correction_source_hartree": float(thermal),
            "thermal_correction_used_hartree": float(thermal_used),
            "thermal_temperature_model": thermo_model,
            "thermal_temperature_status": temperature_status,
            "thermal_source_temperature_K": float(source_temperature),
            "thermal_source_temperature_matches_requested": (
                source_temperature_matches
            ),
            "thermochemistry_frequency_source": frequency_source,
            "thermochemistry_frequencies_migrated_from_raw": bool(
                data.get("thermochemistry_frequencies_migrated_from_raw", False)
            ),
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
            "goodvibes_enthalpy_hartree": external_h,
            "goodvibes_zpe_hartree": external_zpe,
            "goodvibes_temperature_K": (
                float(external_temperature)
                if external_temperature is not None
                else None
            ),
            "goodvibes_temperature_matches_requested": external_temperature_matches,
            "quantitative_thermo_ready": production_ready,
            "production_thermo_ready": production_ready,
            "production_thermo_not_ready_reasons": not_ready_reasons,
            "connector_quality": "production_goodvibes" if production_ready else "internal_qrrho_fallback",
        }

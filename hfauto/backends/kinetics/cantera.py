from __future__ import annotations

"""Cantera-oriented exports and process-screening reactor models.

Phase 10 keeps the safe pseudo-mechanism draft while adding optional production
validation/execution hooks.  If Cantera is installed and explicitly enabled, the
adapter attempts to load the generated or user-supplied mechanism and run a very
small batch-reactor smoke test.  Failures are recorded as data, never hidden.
"""

import csv
import importlib.util
import math
import os
from pathlib import Path
from typing import Any

from hfauto.core.io import ensure_dir, write_json


def _safe_name(value: Any) -> str:
    s = str(value).replace("-", "_").replace("/", "_").replace(".", "_")
    return "".join(ch if ch.isalnum() or ch == "_" else "_" for ch in s)


def _to_float(value: Any, default: float = 0.0) -> float:
    try:
        if value is None or value == "":
            return float(default)
        return float(value)
    except Exception:
        return float(default)


class CanteraEngine:
    name = "cantera"

    def __init__(self, **config: Any):
        self.config = config

    def allow_external(self) -> bool:
        return bool(
            self.config.get("allow_external", False)
            or self.config.get("run_external", False)
            or os.environ.get("HFAUTO_ALLOW_CANTERA") == "1"
        )

    def cantera_available(self) -> bool:
        return importlib.util.find_spec("cantera") is not None

    def write_mechanism_draft(
        self,
        kinetics_records: list[dict[str, Any]],
        arrhenius_records: list[dict[str, Any]],
        out_path: str | Path,
        config: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        path = Path(out_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        arr_by_rxn = {r.get("reaction_id"): r for r in arrhenius_records}
        reaction_ids = sorted({str(r.get("reaction_id")) for r in kinetics_records if r.get("reaction_id")})
        species_names: list[str] = []
        lines: list[str] = [
            "# hfauto Phase 10 Cantera draft mechanism",
            "# Pseudo-species preserve reaction IDs until Arkane/NASA thermo is available.",
            "# This file is valid as a smoke-test mechanism only; review before production use.",
            "units:",
            "  length: cm",
            "  time: s",
            "  quantity: mol",
            "  activation-energy: kcal/mol",
            "",
            "phases:",
            "- name: gas",
            "  thermo: ideal-gas",
            "  elements: [X]",
            "  species: [",
        ]
        for rid in reaction_ids:
            rc = f"RC_{_safe_name(rid)}"
            ip = f"IP_{_safe_name(rid)}"
            species_names.extend([rc, ip])
        lines.append("    " + ", ".join(species_names) + "]")
        lines.extend([
            "  kinetics: gas",
            "  reactions: all",
            "  state: {T: 300.0, P: 1 atm}",
            "",
            "species:",
        ])
        for s in species_names:
            lines.extend([
                f"- name: {s}",
                "  composition: {X: 1}",
                "  thermo:",
                "    model: constant-cp",
                "    h0: 0.0 kcal/mol",
                "    s0: 0.0 cal/mol/K",
                "    cp0: 0.0 cal/mol/K",
            ])
        lines.append("\nreactions:")
        for rid in reaction_ids:
            rc = f"RC_{_safe_name(rid)}"
            ip = f"IP_{_safe_name(rid)}"
            arr = arr_by_rxn.get(rid, {})
            A = max(_to_float(arr.get("A_s-1"), 1.0e0), 0.0)
            Ea = _to_float(arr.get("Ea_kcal_mol"), 0.0)
            b = _to_float(arr.get("b"), 0.0)
            lines.extend([
                f"- equation: {rc} => {ip}",
                f"  rate-constant: {{A: {A:.8e}, b: {b:.6g}, Ea: {Ea:.8g}}}",
                f"  note: hfauto reaction_id={rid}; pseudo first-order RC-to-IP TST step",
            ])
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return {
            "format": "cantera_yaml_draft",
            "path": str(path),
            "n_reactions": len(reaction_ids),
            "n_species": len(species_names),
            "pseudo_species": True,
            "production_note": "Replace pseudo constant-Cp species with Arkane NASA polynomials for production Cantera simulations.",
        }

    def export(self, kinetics_records: list[dict[str, Any]], out_dir: str | Path) -> dict[str, str]:
        out = ensure_dir(out_dir)
        native_path = out / "hfauto_kinetics.yaml"
        cantera_path = out / "cantera_mechanism.yaml"
        json_path = out / "hfauto_kinetics.json"
        reactor_path = out / "reactor_screening.csv"
        reactor_template_path = out / "cantera_reactor_template.py"
        validation_path = out / "cantera_validation.json"
        production_reactor_path = out / "cantera_reactor_results.csv"

        write_json(json_path, {"kinetics_records": kinetics_records, "schema": "hfauto.kinetics.v1"})
        lines = ["schema: hfauto.kinetics.v1", "reactions:"]
        for rec in kinetics_records:
            lines.extend([
                f"  - reaction_id: {rec.get('reaction_id')}",
                f"    mol_id: {rec.get('mol_id')}",
                f"    site_id: {rec.get('site_id')}",
                f"    hf_n: {rec.get('hf_n')}",
                f"    T_K: {rec.get('T_K')}",
                f"    delta_G_act_kcal_mol: {rec.get('delta_G_act_kcal_mol')}",
                f"    k_TST_s-1: {rec.get('k_TST_s-1')}",
                f"    transmission_coefficient: {rec.get('transmission_coefficient')}",
                f"    k_corrected_s-1: {rec.get('k_corrected_s-1')}",
                f"    kinetics_quality: {rec.get('kinetics_quality')}",
            ])
        native_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

        arr = [{"reaction_id": r.get("reaction_id"), "A_s-1": r.get("k_corrected_s-1"), "b": 0.0, "Ea_kcal_mol": 0.0} for r in kinetics_records]
        draft_meta = self.write_mechanism_draft(kinetics_records, arr, cantera_path)

        residence = self.config.get("residence_times_s", [0.001, 0.01, 0.1, 1.0])
        rows = self.simulate_first_order_batch(kinetics_records, residence)
        self._write_rows(reactor_path, rows, ["reaction_id", "mol_id", "site_id", "hf_n", "T_K", "residence_time_s", "k_corrected_s-1", "first_order_conversion_RC_to_IP", "reactor_model"])

        reactor_template_path.write_text(
            "# generated_by=hfauto Phase 10\n"
            "# Requires cantera and production thermo/species definitions for real reactor use.\n"
            "import cantera as ct\n"
            "gas = ct.Solution('cantera_mechanism.yaml')\n"
            "print(gas.n_species, gas.n_reactions)\n",
            encoding="utf-8",
        )

        validation = self.validate_mechanism(cantera_path)
        write_json(validation_path, {**validation, "draft_meta": draft_meta})
        reactor_meta = self.run_cantera_reactor_smoke_test(cantera_path, production_reactor_path)
        # If no real Cantera run happened, write a stable empty CSV.
        if not production_reactor_path.exists():
            production_reactor_path.write_text("status,reason\nnot_run,cantera_not_available_or_not_allowed\n", encoding="utf-8")
        return {
            "hfauto_yaml": str(native_path),
            "cantera_yaml": str(cantera_path),
            "json": str(json_path),
            "reactor_csv": str(reactor_path),
            "reactor_template_py": str(reactor_template_path),
            "cantera_validation_json": str(validation_path),
            "cantera_reactor_results_csv": str(production_reactor_path),
        }

    def _write_rows(self, path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
        with path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                writer.writerow({k: row.get(k) for k in fieldnames})

    def validate_mechanism(self, mechanism_path: str | Path) -> dict[str, Any]:
        if not self.allow_external():
            return {"cantera_available": self.cantera_available(), "cantera_validation_status": "not_requested", "production_cantera_ready": False}
        if not self.cantera_available():
            return {"cantera_available": False, "cantera_validation_status": "module_not_installed", "production_cantera_ready": False}
        try:
            import cantera as ct  # type: ignore
            gas = ct.Solution(str(mechanism_path))
            return {
                "cantera_available": True,
                "cantera_validation_status": "success",
                "production_cantera_ready": bool(not self.config.get("pseudo_species", True)),
                "n_species": int(gas.n_species),
                "n_reactions": int(gas.n_reactions),
                "cantera_version": getattr(ct, "__version__", None),
            }
        except Exception as exc:
            return {"cantera_available": True, "cantera_validation_status": "failed", "production_cantera_ready": False, "reason": str(exc)}

    def run_cantera_reactor_smoke_test(self, mechanism_path: str | Path, out_csv: str | Path) -> dict[str, Any]:
        if not self.allow_external() or not self.config.get("run_reactor", False):
            return {"cantera_reactor_executed": False, "reactor_status": "not_requested"}
        if not self.cantera_available():
            return {"cantera_reactor_executed": False, "reactor_status": "module_not_installed"}
        try:
            import cantera as ct  # type: ignore
            gas = ct.Solution(str(mechanism_path))
            T = float(self.config.get("reactor_temperature_K", 300.0))
            P = float(self.config.get("reactor_pressure_Pa", 101325.0))
            if gas.n_species:
                X = {gas.species_names[0]: 1.0}
                gas.TPX = T, P, X
            reactor = ct.IdealGasReactor(gas)
            sim = ct.ReactorNet([reactor])
            times = [float(x) for x in self.config.get("residence_times_s", [1e-6, 1e-4, 1e-2])]
            rows: list[dict[str, Any]] = []
            for t in times:
                sim.advance(max(t, 0.0))
                row = {"time_s": t, "temperature_K": reactor.T, "pressure_Pa": reactor.thermo.P, "status": "success"}
                for s in gas.species_names[:10]:
                    row[f"X_{s}"] = reactor.thermo[s].X[0]
                rows.append(row)
            fieldnames = sorted({k for r in rows for k in r})
            self._write_rows(Path(out_csv), rows, fieldnames)
            return {"cantera_reactor_executed": True, "reactor_status": "success", "n_rows": len(rows)}
        except Exception as exc:
            Path(out_csv).write_text(f"status,reason\nfailed,{str(exc)!r}\n", encoding="utf-8")
            return {"cantera_reactor_executed": False, "reactor_status": "failed", "reason": str(exc)}

    def simulate_first_order_batch(
        self,
        kinetics_records: list[dict[str, Any]],
        residence_times_s: list[float],
    ) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        for rec in kinetics_records:
            k = rec.get("k_corrected_s-1")
            if k is None:
                continue
            kf = max(float(k), 0.0)
            for tau in residence_times_s:
                tau_f = max(float(tau), 0.0)
                conversion = 1.0 - math.exp(-kf * tau_f) if kf * tau_f < 700 else 1.0
                rows.append(
                    {
                        "reaction_id": rec.get("reaction_id"),
                        "mol_id": rec.get("mol_id"),
                        "site_id": rec.get("site_id"),
                        "hf_n": rec.get("hf_n"),
                        "T_K": rec.get("T_K"),
                        "residence_time_s": tau_f,
                        "k_corrected_s-1": kf,
                        "first_order_conversion_RC_to_IP": conversion,
                        "reactor_model": "pseudo_first_order_batch_RC_to_IP",
                    }
                )
        return rows

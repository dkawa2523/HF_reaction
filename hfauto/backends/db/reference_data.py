from __future__ import annotations

from typing import Any


def _norm(text: str | None) -> str:
    return "".join(ch.lower() for ch in str(text or "") if ch.isalnum())


# Compact offline reference fixture for CI/offline review.  Production runs should
# enrich these records from NIST/ATcT/CCCBDB directly; the fixture is versioned and
# marked as such in every provider result.
LOCAL_REFERENCE_VERSION = "local_reference_v1"

LOCAL_COMPOUND_REFERENCES: dict[str, dict[str, Any]] = {
    "ammonia": {
        "aliases": ["ammonia", "azane", "NH3"],
        "nist_webbook": {"matched": True, "proton_affinity_kj_mol": 853.6, "gas_basicity_kj_mol": 819.0, "has_ir": True, "source_note": "offline fixture; verify against NIST for production"},
        "cccbdb": {"matched": True, "dipole_D": 1.47, "reference_kind": "small_molecule_benchmark_fixture"},
        "comptox": {"matched": True, "gas_process_feasibility": "manual_review_required", "ehs_review_flag": "manual_ehs_review_required"},
        "niosh": {"matched": True, "manual_ehs_review_required": True, "hazard_tags": ["toxic_or_irritant_gas", "corrosive_or_alkaline"]},
    },
    "trimethylamine": {
        "aliases": ["trimethylamine", "N,N-dimethylmethanamine", "TMA"],
        "nist_webbook": {"matched": True, "proton_affinity_kj_mol": 948.9, "gas_basicity_kj_mol": 918.0, "has_ir": True, "source_note": "offline fixture; verify against NIST for production"},
        "cccbdb": {"matched": True, "reference_kind": "small_molecule_benchmark_fixture"},
        "comptox": {"matched": True, "gas_process_feasibility": "candidate_gas_process_review", "ehs_review_flag": "manual_ehs_review_required"},
        "niosh": {"matched": True, "manual_ehs_review_required": True, "hazard_tags": ["flammable_gas", "amine_odor", "corrosive_or_irritant"]},
    },
    "pyridine": {
        "aliases": ["pyridine", "azabenzene"],
        "nist_webbook": {"matched": True, "proton_affinity_kj_mol": 930.0, "gas_basicity_kj_mol": 898.0, "has_ir": True, "source_note": "offline fixture; verify against NIST for production"},
        "cccbdb": {"matched": True, "reference_kind": "small_molecule_benchmark_fixture"},
        "comptox": {"matched": True, "gas_process_feasibility": "liquid_vapor_review", "ehs_review_flag": "manual_ehs_review_required"},
        "niosh": {"matched": True, "manual_ehs_review_required": True, "hazard_tags": ["flammable", "toxic_or_irritant"]},
    },
    "aniline": {
        "aliases": ["aniline", "benzenamine"],
        "nist_webbook": {"matched": True, "proton_affinity_kj_mol": 882.5, "gas_basicity_kj_mol": 851.0, "has_ir": True, "source_note": "offline fixture; verify against NIST for production"},
        "cccbdb": {"matched": False, "reference_kind": "not_in_small_fixture"},
        "comptox": {"matched": True, "gas_process_feasibility": "low_volatility_or_manual_review", "ehs_review_flag": "manual_ehs_review_required"},
        "niosh": {"matched": True, "manual_ehs_review_required": True, "hazard_tags": ["toxic", "manual_review"]},
    },
    "water": {
        "aliases": ["water", "oxidane", "H2O"],
        "nist_webbook": {"matched": True, "proton_affinity_kj_mol": 691.0, "gas_basicity_kj_mol": 660.0, "has_ir": True, "source_note": "offline fixture; verify against NIST for production"},
        "cccbdb": {"matched": True, "dipole_D": 1.85, "reference_kind": "small_molecule_benchmark_fixture"},
    },
    "hydrogenfluoride": {
        "aliases": ["hydrogen fluoride", "HF", "hydrogenfluoride"],
        "nist_webbook": {"matched": True, "has_ir": True, "source_note": "offline fixture for HF reference species"},
        "atct": {"matched": True, "species": "HF", "source_note": "offline fixture; verify against ATcT for production"},
        "niosh": {"matched": True, "manual_ehs_review_required": True, "hazard_tags": ["highly_toxic", "corrosive", "HF"]},
    },
}

_ALIAS_INDEX: dict[str, str] = {}
for key, rec in LOCAL_COMPOUND_REFERENCES.items():
    _ALIAS_INDEX[_norm(key)] = key
    for alias in rec.get("aliases", []):
        _ALIAS_INDEX[_norm(alias)] = key


def lookup_local_reference(name: str | None = None, smiles: str | None = None, inchikey: str | None = None) -> dict[str, Any] | None:
    # Name-based lookup is intentionally preferred for the offline fixture because
    # example SDF names are stable and no network/InChIKey is required.
    for candidate in [name, smiles, inchikey]:
        key = _ALIAS_INDEX.get(_norm(candidate))
        if key:
            out = dict(LOCAL_COMPOUND_REFERENCES[key])
            out["reference_key"] = key
            out["local_reference_version"] = LOCAL_REFERENCE_VERSION
            return out
    return None


def default_reference_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for key, rec in LOCAL_COMPOUND_REFERENCES.items():
        nist = rec.get("nist_webbook", {})
        rows.append({
            "reference_key": key,
            "aliases": ";".join(rec.get("aliases", [])),
            "proton_affinity_kj_mol": nist.get("proton_affinity_kj_mol"),
            "gas_basicity_kj_mol": nist.get("gas_basicity_kj_mol"),
            "has_ir": nist.get("has_ir"),
            "reference_status": "offline_fixture_for_ci; production_runs_should_refresh_from_public_db",
        })
    return rows

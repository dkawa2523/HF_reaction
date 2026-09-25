from __future__ import annotations

from typing import Any


def identity_confidence_score(identity: dict[str, Any] | None) -> float:
    identity = identity or {}
    if identity.get("identity_conflict"):
        return 0.0
    table = {"high": 1.0, "medium": 0.7, "structure_only": 0.45, "unknown": 0.2, None: 0.2}
    return float(table.get(identity.get("identity_confidence"), 0.3))


def reference_availability(public_data: dict[str, Any] | None) -> dict[str, bool]:
    public_data = public_data or {}
    refs = public_data.get("reference_values", {}) or {}
    thermo = public_data.get("thermochemistry_reference", {}) or {}
    bench = public_data.get("benchmark_reference", {}) or {}
    return {
        "has_pa_gb_reference": bool(refs.get("proton_affinity_kj_mol") or refs.get("gas_basicity_kj_mol")),
        "has_thermochemistry_reference": bool(
            refs.get("delta_f_H_298_kj_mol") or thermo.get("delta_f_H_298_kj_mol") or thermo.get("delta_f_G_298_kj_mol")
        ),
        "has_vibrational_reference": bool(refs.get("hf_stretch_cm1") or bench.get("hf_stretch_cm1") or public_data.get("nist_webbook", {}).get("has_ir")),
        "has_process_property_reference": bool(
            public_data.get("gas_process", {}).get("vapor_pressure_Pa_25C")
            or public_data.get("gas_process", {}).get("boiling_point_C")
        ),
    }


def public_db_support_score(public_data: dict[str, Any] | None) -> float:
    public_data = public_data or {}
    refs = reference_availability(public_data)
    score = 0.0
    if public_data.get("pubchem", {}).get("matched"):
        score += 0.25
    if refs["has_pa_gb_reference"]:
        score += 0.20
    if refs["has_thermochemistry_reference"]:
        score += 0.20
    if refs["has_vibrational_reference"]:
        score += 0.15
    if refs["has_process_property_reference"]:
        score += 0.15
    if public_data.get("niosh", {}).get("matched") or public_data.get("comptox", {}).get("matched"):
        score += 0.05
    return max(0.0, min(1.0, score))


def gas_process_flags(public_data: dict[str, Any] | None) -> dict[str, Any]:
    public_data = public_data or {}
    gas = public_data.get("gas_process", {}) or {}
    ehs = gas.get("ehs_review_flag")
    manual = ehs in {True, "manual_review_required", "yes", "unknown_comptox_not_executed"} or bool(
        public_data.get("niosh", {}).get("manual_ehs_review_required")
    )
    return {
        "gas_process_feasibility": gas.get("gas_process_feasibility", "unknown"),
        "ehs_review_flag": ehs or "unknown",
        "manual_ehs_review_required": bool(manual),
        "vapor_pressure_Pa_25C": gas.get("vapor_pressure_Pa_25C"),
        "boiling_point_C": gas.get("boiling_point_C"),
    }


# Demo/offline seeds for workflow tests only.  Production workflows should use
# live/cached providers or a reviewed internal reference snapshot.
STATIC_NIST_PA_GB = {
    "ammonia": {"proton_affinity_kj_mol": 853.6, "gas_basicity_kj_mol": 819.0},
    "trimethylamine": {"proton_affinity_kj_mol": 948.9, "gas_basicity_kj_mol": 918.0},
    "pyridine": {"proton_affinity_kj_mol": 930.0, "gas_basicity_kj_mol": 898.0},
    "dimethyl_ether": {"proton_affinity_kj_mol": 792.0, "gas_basicity_kj_mol": 760.0},
}

STATIC_ATCT = {
    "hydrogen_fluoride": {"delta_f_H_298_kj_mol": -273.3},
    "ammonia": {"delta_f_H_298_kj_mol": -45.9},
}

STATIC_CCCBDB = {
    "hydrogen_fluoride": {"vibrational_frequencies_cm1": [4138.0], "dipole_D": 1.826},
    "ammonia": {"vibrational_frequencies_cm1": [3337.0], "dipole_D": 1.47},
}

STATIC_PUBCHEM = {
    "ammonia": {
        "cid": 222,
        "properties": {
            "MolecularFormula": "H3N",
            "MolecularWeight": 17.031,
            "CanonicalSMILES": "N",
            "IsomericSMILES": "N",
            "InChIKey": "QGZKDVFQNNGYKY-UHFFFAOYSA-N",
            "IUPACName": "azane",
            "XLogP": -0.7,
            "TPSA": 1.0,
            "HBondDonorCount": 1,
            "HBondAcceptorCount": 1,
            "HeavyAtomCount": 1,
            "RotatableBondCount": 0,
        },
    },
    "trimethylamine": {
        "cid": 1146,
        "properties": {
            "MolecularFormula": "C3H9N",
            "MolecularWeight": 59.112,
            "CanonicalSMILES": "CN(C)C",
            "IsomericSMILES": "CN(C)C",
            "InChIKey": "GETQZCLCWQTVFV-UHFFFAOYSA-N",
            "IUPACName": "N,N-dimethylmethanamine",
            "XLogP": 0.2,
            "TPSA": 3.2,
            "HBondDonorCount": 0,
            "HBondAcceptorCount": 1,
            "HeavyAtomCount": 4,
            "RotatableBondCount": 0,
        },
    },
    "pyridine": {
        "cid": 1049,
        "properties": {
            "MolecularFormula": "C5H5N",
            "MolecularWeight": 79.102,
            "CanonicalSMILES": "C1=CC=NC=C1",
            "IsomericSMILES": "C1=CC=NC=C1",
            "InChIKey": "JUJWROOIHBZHMG-UHFFFAOYSA-N",
            "IUPACName": "pyridine",
            "XLogP": 0.7,
            "TPSA": 12.9,
            "HBondDonorCount": 0,
            "HBondAcceptorCount": 1,
            "HeavyAtomCount": 6,
            "RotatableBondCount": 0,
        },
    },
    "dimethyl_ether": {
        "cid": 8254,
        "properties": {
            "MolecularFormula": "C2H6O",
            "MolecularWeight": 46.069,
            "CanonicalSMILES": "COC",
            "IsomericSMILES": "COC",
            "InChIKey": "LCGLNKUTAGEVQW-UHFFFAOYSA-N",
            "IUPACName": "methoxymethane",
            "XLogP": 0.1,
            "TPSA": 9.2,
            "HBondDonorCount": 0,
            "HBondAcceptorCount": 1,
            "HeavyAtomCount": 3,
            "RotatableBondCount": 0,
        },
    },
}


def _norm_static_key(value: Any) -> str:
    import re

    return re.sub(r"[^a-z0-9]+", "", str(value or "").lower())


def first_static_match(data: dict[str, Any], table: dict[str, dict[str, Any]]) -> tuple[str | None, dict[str, Any] | None]:
    keys = [_norm_static_key(data.get("name")), _norm_static_key(data.get("canonical_smiles")), _norm_static_key(data.get("inchikey"))]
    for raw_key, rec in table.items():
        candidates = [_norm_static_key(raw_key)]
        props = rec.get("properties", {}) if isinstance(rec, dict) else {}
        candidates.extend([_norm_static_key(props.get("CanonicalSMILES")), _norm_static_key(props.get("InChIKey"))])
        if any(k and k in candidates for k in keys):
            return raw_key, rec
    return None, None


def process_penalty_from_public_data(public_data: dict[str, Any] | None) -> float:
    """Return a small, reviewable process/EHS penalty from public DB flags."""
    public_data = public_data or {}
    flags = gas_process_flags(public_data)
    feasibility = str(flags.get("gas_process_feasibility") or "").lower()
    penalty = 0.0
    if "low_volatility" in feasibility or "not_suitable" in feasibility:
        penalty += 0.50
    elif "process_condition" in feasibility or "moderate" in feasibility or "review" in feasibility:
        penalty += 0.15
    if flags.get("manual_ehs_review_required"):
        penalty += 0.35
    if public_data.get("niosh", {}).get("manual_ehs_review_required"):
        penalty += 0.20
    return min(1.5, penalty)

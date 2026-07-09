from __future__ import annotations

import re
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Any, Iterable

import yaml


def normalize_key(value: Any) -> str:
    """Normalize names/aliases for permissive public-reference matching."""
    if value is None:
        return ""
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


def load_package_yaml(relative_name: str) -> dict[str, Any]:
    """Load a YAML file shipped under ``hfauto/data``.

    The helper also works from an editable source checkout.  Keeping reference
    data access behind this small function makes it easy to replace bundled
    reference sets by project-curated tables later.
    """
    try:
        text = resources.files("hfauto.data").joinpath(relative_name).read_text(encoding="utf-8")
        return yaml.safe_load(text) or {}
    except Exception:
        path = Path(__file__).resolve().parents[1] / "data" / relative_name
        return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def load_reference_set(path: str | Path | None = None) -> dict[str, Any]:
    if path:
        return yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    return load_package_yaml("reference_hf_amine.yaml")


def load_niosh_flags(path: str | Path | None = None) -> dict[str, Any]:
    if path:
        return yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    return load_package_yaml("niosh_flags.yaml")


def iter_reference_entries(reference_set: dict[str, Any]) -> Iterable[dict[str, Any]]:
    for entry in reference_set.get("references", []) or []:
        yield entry


def match_reference(data: dict[str, Any], reference_set: dict[str, Any] | None = None) -> dict[str, Any] | None:
    """Return the best bundled/reference-set match for a molecule record.

    Match priority: exact InChIKey, normalized name/alias, formula.  Formula-only
    matches are marked by the caller as lower confidence because isomers can share
    formulae.
    """
    refs = reference_set or load_reference_set()
    inchikey = data.get("inchikey") or (data.get("identity") or {}).get("inchikey")
    name = data.get("name") or data.get("sdf_props", {}).get("Name")
    formula = data.get("formula")
    name_key = normalize_key(name)
    formula_key = normalize_key(formula)
    for entry in iter_reference_entries(refs):
        if inchikey and entry.get("inchikey") == inchikey:
            return {**entry, "match_type": "inchikey", "match_confidence": "high"}
    if name_key:
        for entry in iter_reference_entries(refs):
            aliases = [entry.get("name"), *(entry.get("aliases") or [])]
            if name_key in {normalize_key(a) for a in aliases if a}:
                return {**entry, "match_type": "name_alias", "match_confidence": "medium"}
    if formula_key:
        for entry in iter_reference_entries(refs):
            if normalize_key(entry.get("formula")) == formula_key:
                return {**entry, "match_type": "formula", "match_confidence": "low"}
    return None


def match_flag_entry(data: dict[str, Any], flag_set: dict[str, Any] | None = None) -> dict[str, Any] | None:
    flags = flag_set or load_niosh_flags()
    inchikey = data.get("inchikey")
    name_key = normalize_key(data.get("name"))
    for entry in flags.get("entries", []) or []:
        if inchikey and entry.get("inchikey") == inchikey:
            return {**entry, "match_type": "inchikey"}
    for entry in flags.get("entries", []) or []:
        aliases = [entry.get("name"), *(entry.get("aliases") or [])]
        if name_key and name_key in {normalize_key(a) for a in aliases if a}:
            return {**entry, "match_type": "name_alias"}
    return None


@dataclass(frozen=True)
class ReferenceComparison:
    property_name: str
    reference_value: float
    computed_value: float
    error: float
    abs_error: float
    units: str | None = None


def compare_property(property_name: str, reference_value: Any, computed_value: Any, units: str | None = None) -> ReferenceComparison | None:
    try:
        ref = float(reference_value)
        comp = float(computed_value)
    except Exception:
        return None
    err = comp - ref
    return ReferenceComparison(property_name, ref, comp, err, abs(err), units)

from __future__ import annotations

import re
from hashlib import sha256


def slug(text: str) -> str:
    text = re.sub(r"[^A-Za-z0-9_\-]+", "_", str(text).strip())
    return text.strip("_") or "unnamed"


def species_id(*parts: object) -> str:
    return "spc_" + "_".join(slug(str(p)) for p in parts if p is not None and str(p) != "")


def species_artifact_id(species_id: str) -> str:
    """The artifact id of a SpeciesRecord, in every stage that emits or names one."""
    return f"species_{species_id}"


def reaction_id(*parts: object) -> str:
    return "rxn_" + "_".join(slug(str(p)) for p in parts if p is not None and str(p) != "")


def path_token(value: object) -> str:
    """A readable, stable filesystem component of at most 64 characters."""
    token = slug(str(value))
    if len(token) <= 64:
        return token
    return f"{token[:53]}_{sha256(token.encode('utf-8')).hexdigest()[:10]}"

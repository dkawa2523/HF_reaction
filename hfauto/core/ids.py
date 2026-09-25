from __future__ import annotations

import re
from hashlib import sha256


def slug(text: str) -> str:
    text = re.sub(r"[^A-Za-z0-9_\-]+", "_", str(text).strip())
    return text.strip("_") or "unnamed"


def mol_id_from_index(idx: int) -> str:
    return f"mol{idx:05d}"


def site_id(mol_id: str, site_type: str, atom_index: int) -> str:
    return f"{slug(mol_id)}_{slug(site_type)}_{int(atom_index)}"


def conformer_id(mol_id: str, conf_index: int) -> str:
    return f"{slug(mol_id)}_conf{int(conf_index):04d}"


def species_id(*parts: object) -> str:
    return "spc_" + "_".join(slug(str(p)) for p in parts if p is not None and str(p) != "")


def reaction_id(*parts: object) -> str:
    return "rxn_" + "_".join(slug(str(p)) for p in parts if p is not None and str(p) != "")


def path_token(value: object, max_length: int = 64) -> str:
    """Return a readable, stable filesystem component with bounded length."""

    token = slug(str(value))
    if len(token) <= int(max_length):
        return token
    digest = sha256(token.encode("utf-8")).hexdigest()[:10]
    prefix_length = max(1, int(max_length) - len(digest) - 1)
    return f"{token[:prefix_length]}_{digest}"

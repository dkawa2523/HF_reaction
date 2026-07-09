from __future__ import annotations

import re


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

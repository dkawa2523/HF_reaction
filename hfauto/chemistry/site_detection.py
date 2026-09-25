from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class SiteRule:
    site_type: str
    smarts: str
    priority: int
    base_confidence: float
    role: str = "active"  # active | low_priority | excluded
    reason: str | None = None


SITE_RULES = [
    SiteRule("aniline_like", "[NX3;!$([N+]);$([N][c]);!$(NC=O);!$(NS(=O)=O);!$(N=O)]", 1, 0.84),
    SiteRule("aliphatic_amine", "[NX3;!$([N+]);!$(NC=O);!$(NS(=O)=O);!$(N=O);!$([nH])]", 1, 0.95),
    SiteRule("pyridine_like", "[n;H0;!$([n+]);!$([nH])]", 1, 0.92),
    SiteRule("imine", "[NX2;!$([N+])]=[CX3]", 2, 0.82),
    SiteRule("amidine_or_guanidine_N", "[NX3,NX2][CX3](=[NX2,NX3])[NX3,NX2]", 1, 0.88),
    SiteRule("sulfide", "[SX2;!$([S+]);!$(S=O)]", 3, 0.70),
    SiteRule("ether_oxygen", "[OX2;H0;!$(O=C)]", 4, 0.58, role="low_priority"),
    SiteRule("carbonyl_oxygen", "[OX1]=[CX3]", 5, 0.45, role="low_priority"),
]

EXCLUSION_RULES = [
    SiteRule("amide_N_excluded", "[NX3;$(NC=O)]", 99, 0.05, role="excluded", reason="amide_N_low_basicity"),
    SiteRule("sulfonamide_N_excluded", "[NX3;$(NS(=O)=O)]", 99, 0.05, role="excluded", reason="sulfonamide_N_low_basicity"),
    SiteRule("nitro_N_excluded", "[N+](=O)[O-]", 99, 0.05, role="excluded", reason="nitro_N_not_HF_acceptor"),
    SiteRule("pyrrole_N_excluded", "[nH]", 99, 0.10, role="excluded", reason="pyrrole_N_lone_pair_aromatic"),
    SiteRule("quaternary_N_excluded", "[N+]", 99, 0.05, role="excluded", reason="quaternary_N_no_lone_pair"),
]


def _atom_environment(atom) -> dict[str, Any]:
    hyb = str(atom.GetHybridization()).replace("HybridizationType.", "")
    degree = int(atom.GetDegree())
    heavy_degree = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() > 1)
    # Very cheap steric proxy: more heavy neighbors and higher degree around site.
    neighbor_heavy_sum = sum(sum(1 for nn in n.GetNeighbors() if nn.GetAtomicNum() > 1) for n in atom.GetNeighbors())
    steric_proxy = min(1.0, (heavy_degree + 0.25 * neighbor_heavy_sum) / 6.0)
    return {
        "atomic_symbol": atom.GetSymbol(),
        "formal_charge": int(atom.GetFormalCharge()),
        "is_aromatic": bool(atom.GetIsAromatic()),
        "degree": degree,
        "heavy_degree": heavy_degree,
        "hybridization": hyb,
        "steric_proxy": round(float(steric_proxy), 3),
    }


def detect_basic_sites_from_smiles(smiles: str, include_low_priority: bool = True, include_excluded: bool = True) -> list[dict]:
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    excluded_by_atom: dict[int, SiteRule] = {}
    for rule in EXCLUSION_RULES:
        patt = Chem.MolFromSmarts(rule.smarts)
        if patt is None:
            continue
        for match in mol.GetSubstructMatches(patt):
            excluded_by_atom.setdefault(int(match[0]), rule)

    candidates_by_atom: dict[int, dict[str, Any]] = {}
    for rule in SITE_RULES:
        if rule.role == "low_priority" and not include_low_priority:
            continue
        patt = Chem.MolFromSmarts(rule.smarts)
        if patt is None:
            continue
        for match in mol.GetSubstructMatches(patt):
            atom_idx = int(match[0])
            atom = mol.GetAtomWithIdx(atom_idx)
            env = _atom_environment(atom)
            excluded_rule = excluded_by_atom.get(atom_idx)
            excluded = excluded_rule is not None
            rec = {
                "atom_index": atom_idx,
                "site_type": rule.site_type,
                "site_smarts": rule.smarts,
                "priority": 99 if excluded else rule.priority,
                "site_confidence": 0.05 if excluded else rule.base_confidence,
                "local_environment": env,
                "excluded": excluded,
                "exclude_reason": excluded_rule.reason if excluded_rule else None,
                "site_role": "excluded" if excluded else rule.role,
                "basicity_proxy": {
                    "rule_priority": rule.priority,
                    "steric_penalty": env["steric_proxy"],
                    "public_pa_gb_support": False,
                },
            }
            prev = candidates_by_atom.get(atom_idx)
            # Keep highest-priority active rule per atom; exclusions override active candidates.
            if prev is None or rec["priority"] < prev["priority"]:
                candidates_by_atom[atom_idx] = rec

    if include_excluded:
        for atom_idx, rule in excluded_by_atom.items():
            if atom_idx in candidates_by_atom:
                continue
            atom = mol.GetAtomWithIdx(atom_idx)
            candidates_by_atom[atom_idx] = {
                "atom_index": atom_idx,
                "site_type": rule.site_type,
                "site_smarts": rule.smarts,
                "priority": rule.priority,
                "site_confidence": rule.base_confidence,
                "local_environment": _atom_environment(atom),
                "excluded": True,
                "exclude_reason": rule.reason,
                "site_role": "excluded",
                "basicity_proxy": {"rule_priority": rule.priority, "public_pa_gb_support": False},
            }

    return sorted(candidates_by_atom.values(), key=lambda r: (r["excluded"], r["priority"], r["atom_index"]))

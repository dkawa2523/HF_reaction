
from __future__ import annotations

from pathlib import Path
from typing import Any

from hfauto.backends.db.common import load_fixture, lookup_fixture
from hfauto.backends.db.http import CachedHTTPClient


class PubChemProvider:
    """PubChem PUG-REST provider with fixture/cache-first behavior."""
    name = "pubchem"
    PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

    def __init__(self, cache_path: str | None = None, cache_dir: str | None = None, fixture_path: str | None = None, rate_limit_per_sec: float = 4.0, fetch_3d: bool = False, fetch_synonyms: bool = True, network_enabled: bool = False, allow_network: bool | None = None, timeout_s: int = 30, **_: Any):
        if allow_network is not None:
            network_enabled = bool(allow_network)
        self.fixture = load_fixture(fixture_path)
        self.cache_dir = Path(cache_dir or ".hfauto_cache/pubchem")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.fetch_3d = bool(fetch_3d)
        self.fetch_synonyms = bool(fetch_synonyms)
        self.client = CachedHTTPClient(self.name, cache_path or self.cache_dir / "pubchem.sqlite", network_enabled=network_enabled, rate_limit_per_sec=rate_limit_per_sec, timeout_s=timeout_s)

    def _lookup_cids(self, data: dict[str, Any]) -> tuple[list[int], dict[str, Any]]:
        attempts: list[dict[str, Any]] = []
        pairs = []
        if data.get("inchikey"):
            pairs.append(("inchikey", str(data["inchikey"])))
        if data.get("canonical_smiles"):
            pairs.append(("smiles", str(data["canonical_smiles"])))
        if data.get("name"):
            pairs.append(("name", str(data["name"])))
        for kind, value in pairs:
            url = f"{self.PUG}/compound/{kind}/{value}/cids/JSON"
            r = self.client.get_json(url)
            attempts.append({"kind": kind, "ok": r.ok, "reason": r.reason, "from_cache": r.from_cache, "url": r.url})
            if r.ok and isinstance(r.json_data, dict):
                cids = r.json_data.get("IdentifierList", {}).get("CID", []) or []
                if cids:
                    return [int(c) for c in cids], {"attempts": attempts}
        return [], {"attempts": attempts}

    def _properties(self, cid: int) -> dict[str, Any]:
        props = ["MolecularFormula", "MolecularWeight", "CanonicalSMILES", "IsomericSMILES", "InChI", "InChIKey", "IUPACName", "XLogP", "TPSA", "HBondDonorCount", "HBondAcceptorCount", "HeavyAtomCount", "RotatableBondCount", "ExactMass", "MonoisotopicMass", "Complexity"]
        r = self.client.get_json(f"{self.PUG}/compound/cid/{cid}/property/{','.join(props)}/JSON")
        if not r.ok or not isinstance(r.json_data, dict):
            return {"provider_status": "failed", "reason": r.reason, "url": r.url}
        return (r.json_data.get("PropertyTable", {}).get("Properties", [{}]) or [{}])[0]

    def _synonyms(self, cid: int) -> list[str]:
        if not self.fetch_synonyms:
            return []
        r = self.client.get_json(f"{self.PUG}/compound/cid/{cid}/synonyms/JSON")
        if not r.ok or not isinstance(r.json_data, dict):
            return []
        infos = r.json_data.get("InformationList", {}).get("Information", []) or []
        return list((infos[0].get("Synonym", []) if infos else []) or [])[:50]

    def _sdf3d(self, cid: int) -> str | None:
        if not self.fetch_3d:
            return None
        out = self.cache_dir / "sdf3d" / f"CID_{cid}_3d.sdf"
        out.parent.mkdir(parents=True, exist_ok=True)
        if out.exists() and out.stat().st_size > 0:
            return str(out)
        r = self.client.get_text(f"{self.PUG}/compound/cid/{cid}/record/SDF", params={"record_type": "3d"})
        if r.ok and r.text and "$$$$" in r.text:
            out.write_text(r.text, encoding="utf-8")
            return str(out)
        return None

    def _result_from_fixture(self, data: dict[str, Any]) -> dict[str, Any] | None:
        hit = lookup_fixture(self.name, data, self.fixture)
        if not hit:
            return None
        result = dict(hit)
        result.setdefault("matched", True)
        result.setdefault("provider_status", "success")
        result.setdefault("source", "offline_fixture")
        result.setdefault("from_fixture", True)
        return result

    def _result_from_network(self, data: dict[str, Any]) -> dict[str, Any]:
        cids, lookup_status = self._lookup_cids(data)
        if not cids:
            return {"matched": False, "provider_status": "skipped_or_no_match", "lookup_status": lookup_status, "reason": "no_cid_or_network_disabled", "source": "PubChem PUG-REST"}
        cid = cids[0]
        return {"matched": True, "provider_status": "success", "cid": cid, "cids": cids[:20], "properties": self._properties(cid), "synonyms": self._synonyms(cid), "sdf_3d_path": self._sdf3d(cid), "lookup_status": lookup_status, "source": "PubChem PUG-REST"}

    def enrich(self, molecule_artifact):
        data = molecule_artifact.data.copy()
        data.setdefault("identity", {})
        data.setdefault("public_data", {})
        try:
            result = self._result_from_fixture(data) or self._result_from_network(data)
        except Exception as exc:
            result = {"matched": False, "provider_status": "failed", "reason": str(exc), "source": "PubChem PUG-REST"}
        data["public_data"]["pubchem"] = result
        if result.get("matched"):
            props = result.get("properties", {}) or {}
            pubchem_inchikey = props.get("InChIKey")
            identity_conflict = bool(pubchem_inchikey and data.get("inchikey") and pubchem_inchikey != data.get("inchikey"))
            data["identity"].update({"pubchem_cid": result.get("cid"), "identity_confidence": "high" if not identity_conflict and pubchem_inchikey else "medium", "identity_conflict": identity_conflict, "cas_rn": data["identity"].get("cas_rn") or result.get("cas_rn")})
            data["public_data"].setdefault("basic_props", {}).update({"molecular_weight": props.get("MolecularWeight"), "exact_mass": props.get("ExactMass"), "xlogp": props.get("XLogP"), "tpsa": props.get("TPSA"), "hbd": props.get("HBondDonorCount"), "hba": props.get("HBondAcceptorCount"), "rotatable_bonds": props.get("RotatableBondCount"), "complexity": props.get("Complexity")})
            if result.get("sdf_3d_path"):
                data["public_data"].setdefault("structure_seeds", {})["pubchem_3d_sdf"] = result.get("sdf_3d_path")
        else:
            data["identity"].setdefault("identity_confidence", "structure_only")
            data["identity"].setdefault("identity_conflict", False)
        data.setdefault("db_provider_status", {})[self.name] = {"status": "success" if result.get("matched") else result.get("provider_status", "miss"), "matched": bool(result.get("matched")), "cid": result.get("cid"), "reason": result.get("reason"), "source": result.get("source")}
        return data

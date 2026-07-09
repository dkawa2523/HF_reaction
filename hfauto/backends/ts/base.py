from __future__ import annotations

from pathlib import Path
from typing import Any

from hfauto.chemistry.geometry_qc import geometry_qc_from_xyz
from hfauto.chemistry.hf_builder import endpoints_atom_order_match
from hfauto.chemistry.reaction_path_qc import make_midpoint_ts_xyz, reaction_coordinate_value
from hfauto.core.schemas.artifact import Artifact


class TSSearchResult(dict):
    """Dict-compatible TS backend result with attribute access."""

    def __init__(self, artifacts: list[Artifact], record: Any = None, success: bool | None = None, ts_species_id: str | None = None, ts_calc_id: str | None = None):
        if success is None:
            success = not any(a.status.status == "failed" for a in artifacts)
        super().__init__(artifacts=artifacts, record=record, success=success, ts_species_id=ts_species_id, ts_calc_id=ts_calc_id)
        self.artifacts = artifacts
        self.record = record
        self.success = bool(success)
        self.ts_species_id = ts_species_id
        self.ts_calc_id = ts_calc_id
        ts_species = next((a for a in artifacts if a.artifact_type == "species" and a.data.get("state") == "transition_state"), None)
        calc = next((a for a in artifacts if a.artifact_type == "calculation"), None)
        if ts_species is not None:
            self["ts_species"] = ts_species
            if self.ts_species_id is None:
                self.ts_species_id = ts_species.artifact_id
                self["ts_species_id"] = ts_species.artifact_id
        if calc is not None:
            self["calculation"] = calc
            if self.ts_calc_id is None:
                self.ts_calc_id = calc.artifact_id
                self["ts_calc_id"] = calc.artifact_id

    def __getattr__(self, name: str) -> Any:
        try:
            return self[name]
        except KeyError as exc:
            raise AttributeError(name) from exc


class IRCResult(dict):
    def __init__(self, artifacts: list[Artifact], record: Any = None, success: bool | None = None):
        if success is None:
            success = not any(a.status.status == "failed" for a in artifacts)
        super().__init__(artifacts=artifacts, record=record, success=success)
        self.artifacts = artifacts
        self.record = record
        self.success = bool(success)
        self.artifact = next((a for a in artifacts if a.artifact_type == "irc"), artifacts[0] if artifacts else None)

    def __getattr__(self, name: str) -> Any:
        if self.artifact is not None and hasattr(self.artifact, name):
            return getattr(self.artifact, name)
        try:
            return self[name]
        except KeyError as exc:
            raise AttributeError(name) from exc


# Backward-compatible name used by earlier Phase 5 drafts.
TSBackendResult = TSSearchResult


def canonical_species_id(species: Artifact) -> str:
    return str(species.data.get("species_id") or species.data.get("source_species_id") or species.artifact_id)


def midpoint_ts_xyz(reactant: Artifact, product: Artifact, out_xyz: str | Path) -> Path:
    return make_midpoint_ts_xyz(reactant.data.get("xyz_path") or reactant.paths.get("xyz"), product.data.get("xyz_path") or product.paths.get("xyz"), out_xyz)


def make_ts_species_artifact(reaction: Artifact, reactant: Artifact, product: Artifact, ts_xyz: str | Path, ts_id: str, source: str, extra_qc: dict[str, Any] | None = None) -> Artifact:
    data = {
        **reactant.data,
        "species_id": ts_id,
        "source_species_id": ts_id,
        "state": "transition_state",
        "xyz_path": str(ts_xyz),
        "reactant_species_id": canonical_species_id(reactant),
        "product_species_id": canonical_species_id(product),
        "reaction_id": reaction.data.get("reaction_id", reaction.artifact_id),
        "reaction_coordinate": reaction.data.get("reaction_coordinate") or reactant.data.get("reaction_coordinate"),
        "ts_builder": source,
    }
    geom_qc = geometry_qc_from_xyz(data, ts_xyz)
    progress = reaction_coordinate_progress_score(reactant, product, ts_xyz, data.get("reaction_coordinate") or {})
    return Artifact(
        artifact_id=ts_id,
        artifact_type="species",
        parents=[reaction.artifact_id, reactant.artifact_id, product.artifact_id],
        paths={"xyz": str(ts_xyz), "reactant_xyz": reactant.data.get("xyz_path", "") or reactant.paths.get("xyz", ""), "product_xyz": product.data.get("xyz_path", "") or product.paths.get("xyz", "")},
        data=data,
        method={"engine": source, "task": "ts_guess"},
        qc={"ts_guess_source": source, "geometry_qc": geom_qc, **progress, **(extra_qc or {})},
    )


def endpoints_compatible(reactant: Artifact, product: Artifact) -> tuple[bool, dict[str, Any]]:
    rc_xyz = reactant.data.get("xyz_path") or reactant.paths.get("xyz")
    ip_xyz = product.data.get("xyz_path") or product.paths.get("xyz")
    if not rc_xyz or not ip_xyz:
        return False, {"reason": "missing_xyz"}
    try:
        atom_order_ok = endpoints_atom_order_match(rc_xyz, ip_xyz)
    except Exception as exc:
        return False, {"reason": str(exc), "atom_order_ok": False}
    map_ok = reactant.data.get("atom_order_key") == product.data.get("atom_order_key")
    return bool(atom_order_ok and map_ok), {"atom_order_ok": atom_order_ok, "atom_order_key_ok": map_ok}


def reaction_coordinate_progress_score(reactant: Artifact, product: Artifact, ts_xyz: str | Path, reaction_coordinate: dict[str, Any]) -> dict[str, Any]:
    rc_data = {**reactant.data, "reaction_coordinate": reaction_coordinate}
    ip_data = {**product.data, "reaction_coordinate": reaction_coordinate}
    ts_data = {**reactant.data, "state": "transition_state", "reaction_coordinate": reaction_coordinate}
    q_rc = reaction_coordinate_value(rc_data, reactant.data.get("xyz_path") or reactant.paths.get("xyz"))
    q_ip = reaction_coordinate_value(ip_data, product.data.get("xyz_path") or product.paths.get("xyz"))
    q_ts = reaction_coordinate_value(ts_data, ts_xyz)
    if q_rc is None or q_ip is None or q_ts is None or abs(q_ip - q_rc) < 1e-12:
        return {"reaction_coordinate_progress_score": 0.0, "reaction_coordinate_between_endpoints": False, "q_reactant_A": q_rc, "q_product_A": q_ip, "q_ts_A": q_ts}
    lo, hi = sorted([q_rc, q_ip])
    between = lo <= q_ts <= hi
    progress = (q_ts - q_rc) / (q_ip - q_rc)
    score = max(0.0, 1.0 - 2.0 * abs(float(progress) - 0.5))
    return {"reaction_coordinate_progress_score": float(score), "reaction_coordinate_between_endpoints": bool(between), "reaction_coordinate_progress_fraction": float(progress), "q_reactant_A": q_rc, "q_product_A": q_ip, "q_ts_A": q_ts}

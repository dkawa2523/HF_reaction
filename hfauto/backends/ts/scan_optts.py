from __future__ import annotations

from pathlib import Path
from typing import Any

from hfauto.backends.ts.dummy import DummyTSEngine
from hfauto.core.schemas.artifact import Artifact


class ScanOptTSEngine:
    """Constrained scan -> OptTS fallback adapter skeleton.

    The category is implemented as an independent backend so production code can
    later add xTB/r2SCAN scan grids without changing TSSearchStage.  In Phase 5
    it either delegates to explicit dummy fallback or returns a recoverable
    failure Artifact.
    """

    name = "scan_optts"

    def __init__(self, **kwargs: Any):
        self.config = kwargs

    def search_ts(self, reaction: Artifact, reactant: Artifact, product: Artifact, method: dict[str, Any], workdir: str | Path) -> dict[str, Any]:
        if method.get("fallback_to_dummy", self.config.get("fallback_to_dummy", False)):
            return DummyTSEngine().search_ts(reaction, reactant, product, {**self.config, **method, "method_id": method.get("method_id", "scan_optts_dummy")}, Path(workdir) / "dummy_fallback")
        fail = Artifact.failure(
            "ts_failed_" + str(reaction.data.get("reaction_id") or reaction.artifact_id),
            "ts_result",
            "scan_optts production runner is not implemented in Phase 5",
            category="backend_not_implemented",
            parents=[reaction.artifact_id, reactant.artifact_id, product.artifact_id],
            recommended_fallback="orca_nebts_or_enable_dummy_fallback",
        )
        return {"artifacts": [fail], "record": fail.model_dump(), "success": False}

    def run_irc(self, reaction: Artifact, ts_species: Artifact, reactant: Artifact, product: Artifact, method: dict[str, Any], workdir: str | Path) -> dict[str, Any]:
        return DummyTSEngine().run_irc(reaction, ts_species, reactant, product, {**self.config, **method}, workdir)

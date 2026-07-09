from __future__ import annotations

from pathlib import Path
from typing import Any

from hfauto.backends.ts.dummy import DummyTSEngine
from hfauto.core.schemas.artifact import Artifact


class PysisyphusEngine:
    """pysisyphus NEB/GSM/IRC backend skeleton.

    The interface is live and stage-compatible; production integration can write
    YAML inputs and parse pysisyphus trajectories here without touching stage code.
    """

    name = "pysisyphus"

    def __init__(self, **kwargs: Any):
        self.config = kwargs

    def search_ts(self, reaction: Artifact, reactant: Artifact, product: Artifact, method: dict[str, Any], workdir: str | Path) -> dict[str, Any]:
        if method.get("fallback_to_dummy", self.config.get("fallback_to_dummy", False)):
            return DummyTSEngine().search_ts(reaction, reactant, product, {**self.config, **method, "method_id": method.get("method_id", "pysisyphus_dummy")}, Path(workdir) / "dummy_fallback")
        fail = Artifact.failure(
            "ts_failed_" + str(reaction.data.get("reaction_id") or reaction.artifact_id),
            "ts_result",
            "pysisyphus production runner is not implemented in Phase 5",
            category="backend_not_implemented",
            parents=[reaction.artifact_id, reactant.artifact_id, product.artifact_id],
            recommended_fallback="orca_nebts_or_enable_dummy_fallback",
        )
        return {"artifacts": [fail], "record": fail.model_dump(), "success": False}

    def run_irc(self, reaction: Artifact, ts_species: Artifact, reactant: Artifact, product: Artifact, method: dict[str, Any], workdir: str | Path) -> dict[str, Any]:
        return DummyTSEngine().run_irc(reaction, ts_species, reactant, product, {**self.config, **method}, workdir)

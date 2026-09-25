"""report stage (design §4.1 #8, §8.2): ranking, discovery coverage and method panel tables.

Reads the whole view; the HTML page is added in Wave 5 (reporting.html).
"""

from __future__ import annotations

from typing import ClassVar, Literal

from hfauto.core.evidence import Evidence
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.records import (
    ArtifactType,
    DiscoveryRecord,
    MinimumRecord,
    ReactionRecord,
    ReactionThermo,
    ReportRecord,
)
from hfauto.reporting import summary
from hfauto.stages.spec import StageConfig, StageRuntime, StageSpec


class ReportConfig(StageConfig):
    T_K: float | None = None  # None: the first of conditions.temperatures_K
    standard_state: Literal["1atm", "1bar", "1M"] | None = None  # None: the first state
    metric: Literal["dG_act", "dG_rxn"] = "dG_act"


class ReportStage:
    spec: ClassVar[StageSpec] = StageSpec(
        name="report", config=ReportConfig, consumes=(), produces=(ArtifactType.REPORT,)
    )

    def run(self, inputs: Manifest, config: ReportConfig, rt: StageRuntime) -> list[Artifact]:
        reactions = inputs.records(ArtifactType.REACTION, ReactionRecord)
        minima = {m.minimum_id: m for m in inputs.records(ArtifactType.MINIMUM, MinimumRecord)}
        calculations = {
            a.artifact_id: a.payload
            for a in inputs.of(ArtifactType.CALCULATION)
            if isinstance(a.payload, Evidence)
        }
        panel, disagreements = summary.method_panel(calculations, reactions, minima=minima)
        rows = summary.rank_rows(
            reactions,
            inputs.records(ArtifactType.REACTION_THERMO, ReactionThermo),
            summary.participant_notes(reactions, minima, disagreements),
            config.T_K if config.T_K is not None else rt.conditions.temperatures_K[0],
            config.standard_state or rt.conditions.standard_states[0],
            metric=config.metric,
            policy=rt.policy,
        )
        counts = summary.coverage(
            inputs.records(ArtifactType.DISCOVERY, DiscoveryRecord),
            [a for a in inputs.artifacts if a.status == "failed"],
        )
        paths = summary.write_tables(rt.stage_dir, rows, counts, panel)
        record = ReportRecord(
            rows=tuple(rows), tables={name: rt.file_ref(path) for name, path in paths.items()}
        )
        return [
            Artifact(artifact_id=f"report_{rt.stage_id}", type=ArtifactType.REPORT, payload=record)
        ]

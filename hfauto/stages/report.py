"""report stage (design §4.1 #8, §8.2): ranking, discovery coverage and method panel tables,
and the static HTML page (reporting.html) of the whole view.
"""

from __future__ import annotations

from typing import ClassVar, Literal

from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.records import (
    ArtifactType,
    DiscoveryRecord,
    MinimumRecord,
    ReactionRecord,
    ReactionThermo,
    ReportRecord,
)
from hfauto.reporting import html, summary
from hfauto.stages.spec import StageConfig, StageRuntime, StageSpec


class ReportConfig(StageConfig):
    """None: the first (T, state) of the ReactionThermo records in the view."""

    T_K: float | None = None
    standard_state: Literal["1atm", "1bar", "1M"] | None = None


class ReportStage:
    spec: ClassVar[StageSpec] = StageSpec(
        name="report", config=ReportConfig, consumes=(), produces=(ArtifactType.REPORT,)
    )

    def run(self, inputs: Manifest, config: ReportConfig, rt: StageRuntime) -> list[Artifact]:
        reactions = inputs.records(ArtifactType.REACTION, ReactionRecord)
        minima = {m.minimum_id: m for m in inputs.records(ArtifactType.MINIMUM, MinimumRecord)}
        panel = summary.method_panel(inputs.of(ArtifactType.CALCULATION), reactions,
                                     minima=minima)
        thermo = inputs.records(ArtifactType.REACTION_THERMO, ReactionThermo)
        T, state = config.T_K, config.standard_state
        if thermo:
            T = thermo[0].T_K if T is None else T
            state = state or thermo[0].standard_state
        rows = summary.rank_rows(reactions, thermo, T, state)
        counts = summary.coverage(
            inputs.records(ArtifactType.DISCOVERY, DiscoveryRecord),
            [a for a in inputs.artifacts if a.status == "failed"],
        )
        paths = summary.write_tables(rt.stage_dir, rows, counts, panel)
        tables = {name: rt.file_ref(path) for name, path in paths.items()}
        artifact_id = f"report_{rt.stage_id}"
        # The page shows this stage's ranks, so it renders the view plus the table-only record.
        tables_only = Artifact(
            artifact_id=artifact_id,
            type=ArtifactType.REPORT,
            payload=ReportRecord(rows=tuple(rows), tables=tables),
        )
        view = inputs.model_copy(update={"artifacts": [*inputs.artifacts, tables_only]})
        page = html.render(view, rt.stage_dir, load_xyz=rt.load_xyz)
        record = ReportRecord(rows=tuple(rows), tables={**tables, "report.html": rt.file_ref(page)})
        return [Artifact(artifact_id=artifact_id, type=ArtifactType.REPORT, payload=record)]

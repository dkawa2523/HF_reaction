from __future__ import annotations

from hfauto.stages.build_hf import BuildHFStage
from hfauto.stages.calibrate import CalibrateStage
from hfauto.stages.conformers import ConformersStage
from hfauto.stages.connector_audit import ConnectorAuditStage
from hfauto.stages.descriptors import DescriptorsStage
from hfauto.stages.detect_sites import DetectSitesStage
from hfauto.stages.dft_minima import DFTMinimaStage
from hfauto.stages.enrich import EnrichStage
from hfauto.stages.ingest import IngestStage
from hfauto.stages.irc import IRCStage
from hfauto.stages.ops import OpsStage
from hfauto.stages.hpc_plan import HPCPlanStage
from hfauto.stages.kinetics import KineticsStage
from hfauto.stages.preopt import PreoptStage
from hfauto.stages.rank import RankStage
from hfauto.stages.sp import SinglePointStage
from hfauto.stages.thermo import ThermoStage
from hfauto.stages.ts_search import TSSearchStage
from hfauto.stages.viz import VizStage

STAGES = {
    "ingest": IngestStage,
    "enrich": EnrichStage,
    "detect-sites": DetectSitesStage,
    "detect_sites": DetectSitesStage,
    "conformers": ConformersStage,
    "build-hf": BuildHFStage,
    "build_hf": BuildHFStage,
    "preopt": PreoptStage,
    "dft-minima": DFTMinimaStage,
    "dft_minima": DFTMinimaStage,
    "ts-search": TSSearchStage,
    "ts_search": TSSearchStage,
    "irc": IRCStage,
    "sp": SinglePointStage,
    "thermo": ThermoStage,
    "descriptors": DescriptorsStage,
    "kinetics": KineticsStage,
    "calibrate": CalibrateStage,
    "connector-audit": ConnectorAuditStage,
    "connector_audit": ConnectorAuditStage,
    "production-audit": ConnectorAuditStage,
    "production_audit": ConnectorAuditStage,
    "method-calibration": CalibrateStage,
    "method_calibration": CalibrateStage,
    "method-validation": CalibrateStage,
    "method_validation": CalibrateStage,
    "rank": RankStage,
    "viz": VizStage,
    "hpc-plan": HPCPlanStage,
    "hpc_plan": HPCPlanStage,
    "ops-plan": HPCPlanStage,
    "ops_plan": HPCPlanStage,
    "ops": OpsStage,
    "operations": OpsStage,
    "hpc-ops": OpsStage,
    "hpc_ops": OpsStage,
}


def get_stage(name: str):
    if name not in STAGES:
        raise KeyError(f"Unknown stage: {name}")
    return STAGES[name]()

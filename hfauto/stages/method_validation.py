from __future__ import annotations

from hfauto.stages.calibrate import CalibrateStage


class MethodValidationStage(CalibrateStage):
    """Backward-compatible alias for Phase 7 CalibrateStage."""

    name = "method-validation"

from __future__ import annotations

import pytest

from hfauto.core.qc import ts_qc


def test_ts_mode_projection_threshold_authorizes_irc_without_claiming_connectivity() -> None:
    result = ts_qc(
        1,
        -1131.57,
        0.547,
        mode_overlap_threshold=0.50,
        imaginary_frequency_cutoff_cm1=-50.0,
    )

    assert result["ts_validated_by_frequency"] is True
    assert result["ts_target_mode_projection_gate_passed"] is True
    assert "irc_validated" not in result


def test_ts_mode_projection_threshold_is_explicit_and_bounded() -> None:
    assert ts_qc(1, -500.0, 0.49)["ts_validated_by_frequency"] is False
    with pytest.raises(ValueError, match="mode_overlap_threshold"):
        ts_qc(1, -500.0, 1.0, mode_overlap_threshold=0.0)

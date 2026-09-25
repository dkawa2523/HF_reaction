"""Backend-neutral result for one bounded reaction discovery attempt."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass(frozen=True)
class DiscoveryResult:
    driver: str
    success: bool
    product_xyz_path: str | None = None
    reactant_endpoint_xyz_path: str | None = None
    ts_xyz_path: str | None = None
    trajectory_path: str | None = None
    electronic_energy_hartree: float | None = None
    low_level_ts_validated: bool = False
    low_level_irc_connected: bool = False
    imaginary_mode_count: int | None = None
    biased_energy_used_as_barrier: bool = False
    paths: dict[str, str] = field(default_factory=dict)
    data: dict[str, Any] = field(default_factory=dict)
    failure_reason: str | None = None


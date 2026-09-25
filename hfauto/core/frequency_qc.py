"""Shared numerical-noise policy for harmonic frequency validation."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from typing import Any

IMAGINARY_FREQUENCY_CUTOFF_SETTING = "imaginary_frequency_cutoff_cm1"
# NWChem minimum parsing has historically treated frequencies below -1 cm-1
# as physical imaginary modes.  Keep that established convention as the shared
# default so minimum and TS gates cannot silently disagree.
DEFAULT_IMAGINARY_FREQUENCY_CUTOFF_CM1 = -1.0


def resolve_imaginary_frequency_cutoff(value: Any = None) -> float:
    """Return a finite, strictly negative imaginary-frequency cutoff."""

    cutoff = (
        DEFAULT_IMAGINARY_FREQUENCY_CUTOFF_CM1 if value is None else float(value)
    )
    if not math.isfinite(cutoff) or cutoff >= 0.0:
        raise ValueError(
            f"{IMAGINARY_FREQUENCY_CUTOFF_SETTING} must be finite and < 0"
        )
    return cutoff


def imaginary_frequency_cutoff_from_method(
    method: Mapping[str, Any] | None,
) -> float:
    """Resolve the one public cutoff setting from a method dictionary."""

    return resolve_imaginary_frequency_cutoff(
        (method or {}).get(IMAGINARY_FREQUENCY_CUTOFF_SETTING)
    )


def is_significant_imaginary_frequency(
    frequency_cm1: Any,
    cutoff_cm1: Any = DEFAULT_IMAGINARY_FREQUENCY_CUTOFF_CM1,
) -> bool:
    """Whether a finite frequency lies below the numerical-noise boundary."""

    cutoff = resolve_imaginary_frequency_cutoff(cutoff_cm1)
    try:
        frequency = float(frequency_cm1)
    except (TypeError, ValueError):
        return False
    return math.isfinite(frequency) and frequency < cutoff


def significant_imaginary_frequencies(
    frequencies_cm1: Iterable[Any],
    cutoff_cm1: Any = DEFAULT_IMAGINARY_FREQUENCY_CUTOFF_CM1,
) -> list[float]:
    """Return finite frequencies below the shared numerical-noise cutoff."""

    cutoff = resolve_imaginary_frequency_cutoff(cutoff_cm1)
    significant: list[float] = []
    for value in frequencies_cm1:
        try:
            frequency = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(frequency) and frequency < cutoff:
            significant.append(frequency)
    return significant

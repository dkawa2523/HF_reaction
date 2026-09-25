"""Evidence gates for stationary minima on one electronic-structure surface."""

from __future__ import annotations

from hfauto.core.schemas.artifact import Artifact


def is_accepted_optimized_minimum(species: Artifact) -> bool:
    """Return whether an optimized species carries a complete minimum claim."""

    return bool(
        species.artifact_type == "species_optimized"
        and species.status.status == "success"
        and species.qc.get("is_minimum") is True
        and species.qc.get("minimum_accepted") is True
        and species.qc.get("scf_converged") is True
        and species.qc.get("geometry_converged") is True
        and species.qc.get("geometry_sane") is True
        and species.qc.get("fallback_dummy") is False
    )


def has_frequency_validated_minimum_evidence(
    species: Artifact,
    calculation: Artifact | None,
    *,
    require_real_qm: bool = True,
) -> bool:
    """Require matching opt/freq evidence before registry admission."""

    if not is_accepted_optimized_minimum(species) or calculation is None:
        return False
    n_imag = calculation.data.get(
        "n_imag", calculation.qc.get("n_imag")
    )
    frequency_complete = calculation.data.get(
        "frequency_count_complete",
        calculation.qc.get("frequency_count_complete"),
    )
    real_qm_ok = (
        calculation.qc.get("real_qm_executed") is True
        if require_real_qm
        else calculation.qc.get("fallback_dummy") is not True
    )
    return bool(
        calculation.status.status == "success"
        and calculation.qc.get("is_minimum") is True
        and calculation.qc.get("minimum_accepted") is True
        and n_imag == 0
        and frequency_complete is True
        and calculation.qc.get("fallback_dummy") is False
        and real_qm_ok
    )

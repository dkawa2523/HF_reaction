from __future__ import annotations

import hashlib
import json

import pytest
from pydantic import ValidationError

from hfauto.core.constants import BOHR_TO_ANGSTROM, CM1_TO_HARTREE
from hfauto.core.evidence import Evidence, Failure, FailureKind, FileRef, Geometry, Level

LEVEL = Level(program="NWChem", version="7.2.3", method="PBE0", basis="def2-SVPD",
              dispersion="d3bj", charge=0, multiplicity=1, grid="fine", scf_tol=1e-7)


def test_level_keys():
    assert (LEVEL.program, LEVEL.method, LEVEL.basis) == ("nwchem", "pbe0", "def2-svpd")
    canonical = json.dumps(LEVEL.model_dump(mode="json"), sort_keys=True, separators=(",", ":"))
    assert LEVEL.full_key() == hashlib.sha256(canonical.encode()).hexdigest()[:16]
    numerics = LEVEL.model_copy(update={"grid": "xfine", "scf_tol": 1e-8})
    assert numerics.full_key() != LEVEL.full_key()


def test_models_are_frozen_and_forbid_extra_fields():
    with pytest.raises(ValidationError):
        LEVEL.method = "b3lyp"  # type: ignore[misc]
    with pytest.raises(ValidationError):
        FileRef(path="a", sha256="b", size=1)  # type: ignore[call-arg]


@pytest.mark.parametrize(("freqs", "modes", "valid"), [
    ((-900.0, 1000.0, 3000.0), 1, True), ((1000.0, 3000.0), 0, False),  # 3N - 6 = 3 for HCN
    ((-900.0, 1000.0, 3000.0), 0, False), (None, 0, False)])
def test_a_freq_evidence_holds_3n_minus_external_modes_and_its_imaginary_ones(freqs, modes, valid):
    ref = FileRef(path="a", sha256="0")
    geo = Geometry(file=ref, fingerprint="f", symbols=("H", "C", "N"))
    make = lambda task: Evidence(
        engine="nwchem", task=task, level=LEVEL, start=geo, final=geo, energy_hartree=-93.0,
        frequencies_cm1=freqs, n_external=6, imaginary_modes=((0.1,) * 9,) * modes, output=ref,
        job_key="k")
    stored = make("sp").model_dump(mode="json")  # only a freq is checked
    del stored["gradient"]  # a JobStore record written before Evidence.gradient
    assert Evidence.model_validate(stored).gradient is None
    if valid:
        make("freq")
    else:
        with pytest.raises(ValidationError):
            make("freq")


def test_failure_kinds_and_constants():
    assert len(FailureKind) == 10
    failure = Failure(kind="gate_rejected", reason="no_collision_free_seed")
    assert failure.kind is FailureKind.GATE_REJECTED and failure.energy_hartree is None
    assert CM1_TO_HARTREE * 219474.6313705 == pytest.approx(1.0, rel=1e-9)
    assert BOHR_TO_ANGSTROM == pytest.approx(0.529177210903)

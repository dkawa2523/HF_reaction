from __future__ import annotations

import hashlib
import json

import pytest
from pydantic import ValidationError

from hfauto.core.constants import BOHR_TO_ANGSTROM, CM1_TO_HARTREE
from hfauto.core.evidence import Failure, FailureKind, FileRef, Level

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


def test_failure_kinds_and_constants():
    assert len(FailureKind) == 11
    failure = Failure(kind="gate_rejected", reason="no_collision_free_seed")
    assert failure.kind is FailureKind.GATE_REJECTED
    assert CM1_TO_HARTREE * 219474.6313705 == pytest.approx(1.0, rel=1e-9)
    assert BOHR_TO_ANGSTROM == pytest.approx(0.529177210903)

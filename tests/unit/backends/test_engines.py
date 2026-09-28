"""Engine registry: the final table, KeyError for unknown keys, test overrides (§6.2)."""

import fakes
import pytest

from hfauto.backends import engines
from hfauto.backends.protocols import Capability as C
from hfauto.backends.protocols import Requirements
from hfauto.core.method import EngineSite


def create(capability, name):
    return engines.create(capability, name, jobs=None, site=EngineSite(version="0"))


def test_table_holds_exactly_the_seven_final_entries() -> None:
    nw = "hfauto.backends.nwchem.engine"
    assert engines._TABLE == {
        C.QM: {"nwchem": f"{nw}:NWChemEngine", "xtb": "hfauto.backends.xtb:XTBEngine"},
        C.PATH: {"nwchem_string": f"{nw}:NWChemString",
                 "pysis_neb": "hfauto.backends.pysis.engine:PysisNEB"},
        C.SADDLE: {"nwchem_saddle": f"{nw}:NWChemSaddle"},
        C.CONFORMERS: {"crest": "hfauto.backends.crest:CRESTEngine"},
        C.DISCOVERY: {"readuct": "hfauto.backends.readuct.engine:ReaDuctEngine"},
    }


def test_unknown_keys_raise_key_error() -> None:
    with pytest.raises(KeyError):
        create(C.SADDLE, "xtb")  # a registered name, but not for this capability
    with pytest.raises(KeyError):
        engines.requirements_for({C.DISCOVERY: ["internal"]})
    with pytest.raises(KeyError), engines.override(C.QM, "dummy", fakes.FakeQM):
        pass


def test_overrides_nest_restore_and_check_the_protocol(tmp_run, override_engine) -> None:
    first, second = (fakes.FakeQM(tmp_run, fakes.harmonic()) for _ in range(2))
    with engines.override(C.QM, "xtb", lambda **_: first):
        with engines.override(C.QM, "xtb", lambda **_: second):
            assert create(C.QM, "xtb") is second
        assert create(C.QM, "xtb") is first
        assert engines.requirements_for({C.QM: ["xtb"]}) == {"xtb": Requirements()}
    assert not engines._OVERRIDES
    override_engine(C.QM, "nwchem", fakes.FakePath(tmp_run, fakes.harmonic()))
    with pytest.raises(TypeError, match="QMEngine"):
        create(C.QM, "nwchem")
    override_engine(C.DISCOVERY, "readuct", fakes.FakeDiscovery(["scripted"]))
    discovery = create(C.DISCOVERY, "readuct")
    assert discovery.explore(None, None, None, None) == "scripted"
    assert discovery.calls == [(None, None, None, None)]
    with pytest.raises(AssertionError, match="script exhausted"):
        discovery.explore(None, None, None, None)

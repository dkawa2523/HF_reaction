"""Engine registry: the final table, KeyError for unknown keys, the Protocol check (§6.2);
the NWChem command line."""

import fakes
import pytest

from hfauto.backends import engines
from hfauto.backends.protocols import Capability as C
from hfauto.core.method import EngineSite, ExecutionSpec
from hfauto.execution.jobs import Task


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


def test_requirements_come_from_the_classes() -> None:
    found = engines.requirements_for({C.QM: ["xtb"], C.SADDLE: ["nwchem_saddle"]})
    assert found["xtb"].version_command == ("xtb", "--version")
    assert found["nwchem_saddle"].executables == ("nwchem",)


def test_create_checks_the_capability_protocol(tmp_run, monkeypatch) -> None:
    path = fakes.FakePath(tmp_run, fakes.harmonic())
    monkeypatch.setattr(engines, "_load", lambda target: lambda **_: path)
    with pytest.raises(TypeError, match="QMEngine"):
        create(C.QM, "nwchem")
    assert create(C.PATH, "nwchem_string") is path


def test_nwchem_with_fewer_ranks_than_the_site_runs_unbound() -> None:
    exe = {"nwchem": "hfauto-no-nwchem", "mpirun": "hfauto-no-mpirun"}
    site = EngineSite(version="7.2.3", executables=exe, execution=ExecutionSpec(ranks=4))
    qm = engines.create(C.QM, "nwchem", jobs=None, site=site)

    def argv(ranks):
        return qm._argv(Task("nwchem", "7.2.3", "energy", {}, ExecutionSpec(ranks=ranks)))

    assert argv(4) == ("hfauto-no-mpirun", "-np", "4", "hfauto-no-nwchem", "job.nw")
    assert argv(2) == ("hfauto-no-mpirun", "-np", "2", "--bind-to", "none",
                       "hfauto-no-nwchem", "job.nw")

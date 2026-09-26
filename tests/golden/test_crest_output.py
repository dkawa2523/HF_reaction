"""G18 / G19: CREST ensembles, missing energies and topology stops (§10.2, BUG-08, CH-18)."""

import shutil

import pytest

from hfauto.backends.crest import TOPOLOGY_STOP, parse_outputs
from hfauto.chemistry.xyz import read_xyz
from hfauto.core.constants import HARTREE_TO_KCAL_MOL
from hfauto.core.evidence import FailureKind
from hfauto.execution.jobstore import JobStore

pytestmark = pytest.mark.golden
TMA_HF2 = ["C", "N", "C", "C", *["H"] * 10, "F", "H", "F"]


def _parse(tmp_path, stdout, symbols):
    return parse_outputs(tmp_path, stdout, symbols, JobStore(tmp_path / "jobs").file_ref)


def test_g18_conformers_energies_and_a_missing_energy(tmp_path, golden):
    shutil.copy(golden.path("crest/G18/crest_conformers.xyz"), tmp_path)
    stdout = golden.text("crest/G18/crest.stdout")
    ens = _parse(tmp_path, stdout, TMA_HF2)
    energies = [e for _, e in ens.members]
    relative = [float(r.split()[1]) for r in golden.text("crest/G18/crest.energies").splitlines()]
    assert [(e - energies[0]) * HARTREE_TO_KCAL_MOL for e in energies] == pytest.approx(
        relative, abs=0.01) and energies[0] == pytest.approx(-24.40971451)
    assert (ens.version, ens.topology_removed, ens.topology_stops) == ("3.0.2", 0, ())
    lines = (tmp_path / "crest_conformers.xyz").read_text(encoding="utf-8").split("\n")
    lines[20] = ""  # comment line of the second conformer (BUG-08: never 0 kcal/mol)
    (tmp_path / "crest_conformers.xyz").write_text("\n".join(lines), encoding="utf-8")
    assert _parse(tmp_path, stdout, TMA_HF2).members[1][1] is None


def test_g19_topology_change_is_a_topology_stop(tmp_path, golden):
    stdout = golden.text("crest/G19/crest.stdout")
    missing = _parse(tmp_path, stdout, ["N", "H"])  # a stop without crestopt.log
    assert TOPOLOGY_STOP in stdout and missing.kind is FailureKind.INCOMPLETE_OUTPUT
    frame = "2\n Etot= {}\nN 0 0 0\nH 0 0 {}\n"  # crestopt.log: the last frame is the stop
    (tmp_path / "crestopt.log").write_text(frame.format(-1.0, 1.1) + frame.format(-1.2, 1.02))
    ens = _parse(tmp_path, stdout, ["N", "H"])
    assert ens.members == () and len(ens.topology_stops) == 1
    assert read_xyz(tmp_path / ens.topology_stops[0].file.path).coords[1, 2] == pytest.approx(1.02)

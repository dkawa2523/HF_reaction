"""G09, G11, G12: real pysisyphus IRC branches judged by gates.connection (CH-03, CH-04)."""

import re
from itertools import accumulate

import pytest

from hfauto.chemistry.gates import connection
from hfauto.chemistry.identity import assign
from hfauto.chemistry.xyz import read_xyz
from hfauto.core.evidence import Evidence, FileRef, Geometry, Level

pytestmark = pytest.mark.golden
LEVEL = Level(program="nwchem", version="7.2.3", method="pbe0", charge=0, multiplicity=1)
GEO = Geometry(file=FileRef(path="x.xyz", sha256="0" * 64), fingerprint="f", symbols=("H",))


def _ev(task, e):
    return Evidence(engine="nwchem", task=task, level=LEVEL, start=GEO, final=GEO, output=GEO.file,
                    energy_hartree=e[-1], trajectory_energies_hartree=tuple(e), job_key=task)


def branches(text):
    """TS freq and both branches: the displaced start (E_TS − actual lowering) + each step's dE."""
    e_ts = float(re.search(r"TS:\s+(-?[\d.]+) hartree", text).group(1))
    sides = []
    for part in text.split("# IRC - ")[1:]:
        lowering = float(re.search(r"Actual energy lowering:\s+(-?[\d.]+) au", part).group(1))
        steps = re.findall(r"(?m)^\s+\d+\s+-?[\d.]+\s+(-?[\d.]+)\s+[\d.]+\s+[\d.]+\s*$", part)
        sides.append(_ev("opt", list(accumulate(map(float, steps), initial=e_ts - lowering))))
    return _ev("freq", [e_ts]), (sides[0], sides[1])


@pytest.mark.parametrize("gid,why", [("G09", "side0:no_initial_descent"),
                                     ("G11", "side1:no_descent")])
def test_irc_branches_that_climb_fail_the_connection_gate(golden, gid, why) -> None:
    ts, sides = branches(golden.text(f"pysisyphus/{gid}/pysis_irc.out"))
    gate, label = connection(ts, sides, ("a", "b"), frozenset({"a", "b"}), degenerate=False)
    assert label == "failed" and why in gate.reasons


def test_g12_hono_branches_fail_and_their_assignment_is_unique(golden) -> None:
    ts, sides = branches(golden.text("pysisyphus/G12/pysis_irc.out"))
    gate, _ = connection(ts, sides, ("trans", "cis"), frozenset({"trans", "cis"}), degenerate=False)
    assert {"side0:no_initial_descent", "side1:no_initial_descent"} <= set(gate.reasons)
    candidates = {n: (read_xyz(golden.path(f"nwchem/G03/hono_{n}_final.xyz")).coords, float(
        re.findall(r"Total DFT energy =\s+(\S+)", golden.text(f"nwchem/G03/hono_{n}.out"))[-1]))
        for n in ("trans", "cis")}
    for branch, side, want in (("forward", sides[0], "trans"), ("backward", sides[1], "cis")):
        x = read_xyz(golden.path(f"pysisyphus/G12/{branch}_final.xyz"))  # 0.64 Å from the other
        assert assign(x.symbols, x.coords, side.energy_hartree, candidates) == want  # CH-04

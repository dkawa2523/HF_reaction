"""Gates on real outputs. G09, G11, G12: pysisyphus IRC branches judged by gates.connection
(CH-03, CH-04): a side counts by where it ends, not by the path there (U6-P5). G25, G26, G28:
NWChem opt and freq pairs judged by gates.is_minimum (U3-P7)."""

import re
from itertools import accumulate

import pytest

from hfauto.backends.nwchem import output as nw
from hfauto.chemistry.gates import Gate, connection, is_minimum
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


@pytest.mark.parametrize("gid,reasons", [
    ("G09", ("side0:no_descent", "side1:no_descent")),  # both branches end above the TS
    ("G11", ()),  # the backward branch rises at step 0 but ends below E_TS - drop
])
def test_a_branch_is_judged_by_where_it_ends(golden, gid, reasons) -> None:
    ts, sides = branches(golden.text(f"pysisyphus/{gid}/pysis_irc.out"))
    gate, label = connection(ts, sides, ("a", "b"), frozenset({"a", "b"}), degenerate=False)
    assert gate.reasons == reasons and label == ("failed" if reasons else "elementary")


def test_g12_hono_branches_start_above_the_ts_and_connect_trans_and_cis(golden) -> None:
    ts, sides = branches(golden.text("pysisyphus/G12/pysis_irc.out"))  # displaced 1-2 mEh up
    assert all(side.trajectory_energies_hartree[0] > ts.energy_hartree for side in sides)
    assert connection(ts, sides, ("trans", "cis"), frozenset({"trans", "cis"}),
                      degenerate=False) == (Gate(True), "elementary")
    candidates = {n: (read_xyz(golden.path(f"nwchem/G03/hono_{n}_final.xyz")).coords, float(
        re.findall(r"Total DFT energy =\s+(\S+)", golden.text(f"nwchem/G03/hono_{n}.out"))[-1]))
        for n in ("trans", "cis")}
    for branch, side, want in (("forward", sides[0], "trans"), ("backward", sides[1], "cis")):
        x = read_xyz(golden.path(f"pysisyphus/G12/{branch}_final.xyz"))  # 0.64 Å from the other
        assert assign(x.symbols, x.coords, side.energy_hartree, candidates) == want  # CH-04


@pytest.mark.parametrize(("stem", "n_external", "gap"), [
    ("G25/oh", 5, 9.6e-7), ("G26/fhf", 5, 3.0e-8), ("G28/cl", 3, 0.0)])
def test_freq_reproduces_its_opt_energy_in_the_same_state(golden, stem, n_external, gap):
    """OH. (UKS doublet), FHF- and Cl-: NWChem's freq, which tightens the grid and screening,
    stays within the SCF noise floor of its opt; the opt shifted by 1e-4 Eh does not."""
    calcs = {}
    for task in ("opt", "freq"):
        text = golden.text(f"nwchem/{stem}_{task}.out")
        symbols, _ = nw.geometry_block(text, -1)
        n_modes = 3 * len(symbols) - n_external
        geo = Geometry(file=GEO.file, fingerprint="final", symbols=symbols)
        calcs[task] = Evidence(
            engine="nwchem", task=task, level=nw.observe_level(text), start=geo, final=geo,
            output=GEO.file, energy_hartree=nw.total_energy(text), job_key=task,
            frequencies_cm1=(1000.0,) * n_modes if task == "freq" else None,  # not under test
            n_external=n_external if task == "freq" else None)
    freq, opt = calcs["freq"], calcs["opt"]
    assert abs(freq.energy_hartree - opt.energy_hartree) == pytest.approx(gap, abs=1e-8)
    assert is_minimum(freq, opt=opt) == Gate(True)
    shifted = opt.model_copy(update={"energy_hartree": opt.energy_hartree - 1e-4})
    assert is_minimum(freq, opt=shifted).reasons == ("state_mismatch",)

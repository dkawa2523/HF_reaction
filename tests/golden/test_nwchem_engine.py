"""NWChem engines through JobRunner with a stub executable that replays recorded outputs:
G07 (HCN TS: level lines, vibrational block and .hess), G01 (two blocks), G13 (string) and
G29 (OH. ROHF-CCSD(T)). A driver job prints a gradient block at its last frame and, with
STUB_MAXITER in its environment, stops at maxiter; a deck with a zcoord block and STUB_AUTOZ
fails in autoz (G33's message)."""

import json
import shutil
import sys

import numpy as np
import pytest

from hfauto.backends.nwchem.engine import NWChemEngine, NWChemSaddle, NWChemString
from hfauto.backends.nwchem.input import hess_text
from hfauto.backends.nwchem.output import classify_failure, geometry_block, read_hess
from hfauto.chemistry.vibrations import KAPPA_MIN, external_basis, shape_hessian
from hfauto.chemistry.xyz import XYZ, Molecule, read_xyz, write_xyz_trajectory
from hfauto.core.constants import BOHR_TO_ANGSTROM
from hfauto.core.evidence import Evidence, Failure, FailureKind, FileRef, PathProfile
from hfauto.core.hashing import sha256_file, sha256_text
from hfauto.core.method import EngineSite, MethodSpec
from hfauto.execution.jobs import JobRunner, Task
from hfauto.execution.jobstore import JobStore

pytestmark = pytest.mark.golden
FINE = MethodSpec(id="fine", kind="dft", functional="pbe0", basis="def2-svpd",
                  dispersion="d3zero", grid="fine", scf_energy_tol=1e-7)  # G07's level
XFINE = FINE.model_copy(update={"grid": "xfine", "scf_energy_tol": 1e-8})  # G01 / G13
STUB = r'''
import os, re, shutil, sys
from pathlib import Path
here, deck = Path(__file__).parent, Path(sys.argv[-1]).read_text()
block = re.search(r"^geometry[^\n]*\n(.*?)^end", deck, re.M | re.S).group(1)
atoms = [row for row in block.splitlines() if len(row.split()) == 4]  # not the zcoord lines
task = deck.rsplit("task ", 1)[1].strip()
frame = lambda rows: f"{len(rows)}\n geometry\n" + "\n".join(rows) + "\n"
if os.environ.get("STUB_AUTOZ") and "zcoord" in block:
    sys.exit(print(" !! There are insufficient internal variables: expected    11 got    12"))
if task == "dft string":
    Path("job.string_final.xyz").write_text(frame(atoms) * 11)
    sys.exit(print((here / "G13.out").read_text()))
if task == "tce energy":
    sys.exit(print((here / "G29.out").read_text()))
if len(atoms) == 4:
    sys.exit(print((here / "G01.out").read_text()))
g07 = (here / "G07.out").read_text()
echo = "".join(f"{i:5d} {r.split()[0]} 1.0 {' '.join(r.split()[1:])}\n"
               for i, r in enumerate(atoms, 1))
print("".join(g07.splitlines(True)[:2]) + ' Geometry "geometry" -> ""\n'
      + " Output coordinates in angstroms\n" + echo + "".join(g07.splitlines(True)[13:50])
      + " Total DFT energy =  -93.166994283613\n")
if task == "dft frequencies":
    shutil.copyfile(here / "G07.hess", "job.hess")
    print(g07[g07.index("  Vibrational analysis"):g07.index(" Task  times")])
if task in ("dft optimize", "dft saddle"):
    restarted = "vectors input" in deck
    moved = [f"{r.split()[0]} " + " ".join(f"{float(v) + 0.01:.8f}" for v in r.split()[1:])
             for r in atoms]
    frames = [atoms] * (3 if restarted else 1) + ([] if restarted else [moved])
    for n, rows in enumerate(frames):
        Path(f"final-{n:03d}.xyz").write_text(frame(rows))
    Path("job.movecs").write_text("v"), Path("job.drv.hess").write_text("h")
    print(" DFT ENERGY GRADIENTS\n" + "".join(
        f"{i:4d} {r.split()[0]} " + " ".join(f"{float(v) / 0.529177210903:.6f}" for v in
                                             r.split()[1:]) + " 0.000010 0.000000 0.000000\n"
        for i, r in enumerate(frames[-1], 1)))
    print("@    0    -93.16699428\n@    1    -93.16699428\n")
    if os.environ.get("STUB_MAXITER"):
        sys.exit(print(" Failed to converge in maximum number of steps"))
    if not restarted:
        sys.exit(124)
    print("      Optimization converged\n")
print(" Total times  cpu:        1.0s")
'''


@pytest.fixture
def nwchem(golden, tmp_path):
    """(JobRunner, EngineSite) whose executables['nwchem'] runs STUB with sys.executable."""
    stub = tmp_path / "stub"
    stub.mkdir()
    for name, rel in (("G01.out", "G01/nwchem.out"), ("G13.out", "G13/nwchem_string.out"),
                      ("G07.out", "G07/nwchem.out"), ("G07.hess", "G07/hfauto_job.hess"),
                      ("G29.out", "G29/oh_ccsdt.out")):
        shutil.copyfile(golden.path(f"nwchem/{rel}"), stub / name)
    (stub / "nwchem.py").write_text(STUB)
    if sys.platform == "win32":
        exe = stub / "nwchem.cmd"
        exe.write_text(f'@"{sys.executable}" "{stub / "nwchem.py"}" %*\n')
    else:
        exe = stub / "nwchem"
        exe.write_text(f'#!/bin/sh\nexec "{sys.executable}" "{stub / "nwchem.py"}" "$@"\n')
        exe.chmod(0o755)
    site = EngineSite(version="7.2.3", executables={"nwchem": str(exe)})
    return JobRunner(JobStore(tmp_path / "run" / "jobs"), cores=1), site


def _hcn_ts(golden) -> Molecule:
    return Molecule(read_xyz(golden.path("nwchem/G07/final.xyz")), 0, 1)


def test_frequencies_cache_and_input_checks(nwchem, golden):
    jobs, site = nwchem
    engine, mol = NWChemEngine(jobs=jobs, site=site), _hcn_ts(golden)
    ev = engine.frequencies(mol, FINE)
    assert isinstance(ev, Evidence) and ev.task == "freq" and ev.n_external == 6
    assert ev.frequencies_cm1[0] == pytest.approx(-1131.6, abs=1.0)
    assert len(ev.imaginary_modes) == 1
    assert np.load(jobs.store.run_dir / ev.hessian.path).shape == (9, 9)
    assert engine.frequencies(mol, FINE) == ev and jobs.stats().hits == 1
    far = Molecule(XYZ(mol.xyz.symbols, mol.xyz.coords + [0.6, 0, 0]), 0, 1)
    failure = engine.optimize(far, FINE, init_hessian=ev)  # beyond 0.5 Å per atom
    assert failure.kind is FailureKind.INPUT_INVALID
    assert failure.reason == "hessian_geometry_mismatch"
    triplet = engine.frequencies(Molecule(mol.xyz, 0, 3), FINE)  # ported K case: state conflict
    assert triplet.kind is FailureKind.METHOD_MISMATCH and "multiplicity" in triplet.reason
    grid = engine.frequencies(mol, XFINE)
    assert grid.kind is FailureKind.METHOD_MISMATCH and "grid" in grid.reason


def test_a_freq_starts_from_the_converged_vectors_of_its_opt(nwchem, golden):
    """The freq at an opt's final structure reads the opt's job.movecs under its own key: it
    stays on the opt's electronic state (from scratch, the UKS OH···CH4 complex found the other
    OH π component)."""
    jobs, site = nwchem
    engine = NWChemEngine(jobs=jobs, site=site)
    opt = engine.optimize(_hcn_ts(golden), FINE)
    final = Molecule(read_xyz(jobs.store.resolve(opt.final.file)), 0, 1)
    guided, plain = engine.frequencies(final, FINE, scf_guess=opt), engine.frequencies(final, FINE)
    assert isinstance(guided, Evidence) and guided.job_key != plain.job_key
    decks = [(jobs.store.attempt_dir(ev.job_key, 0) / "job.nw").read_text()
             for ev in (guided, plain)]
    assert "vectors input job.movecs" in decks[0] and "vectors input" not in decks[1]


def test_a_first_order_saddle_hessian_starts_a_minimization_as_its_positive_definite_model(
        nwchem, golden):
    """K1: the TS Hessian of G07 (one mode below -saddle_cm1) is written as its positive-definite
    model under a new job key, and so is the same Hessian with no saddle mode (a minimum's or an
    R6 seed's, X1); with a second saddle mode (a higher-order saddle) it is written as it is
    under the old key."""
    jobs, site = nwchem
    engine, mol = NWChemEngine(jobs=jobs, site=site), _hcn_ts(golden)
    ts = engine.frequencies(mol, FINE)
    side = Molecule(XYZ(mol.xyz.symbols, mol.xyz.coords + [0.1, 0, 0]), 0, 1)
    raw = np.load(jobs.store.resolve(ts.hessian))
    higher = ts.model_copy(update={"frequencies_cm1": (-120.0, *ts.frequencies_cm1)})
    none = ts.model_copy(update={"frequencies_cm1": tuple(abs(f) for f in ts.frequencies_cm1)})
    model = shape_hessian(raw, side.xyz.coords)
    decks = []
    for freq, written, reshaped in ((ts, model, True), (none, model, True), (higher, raw, False)):
        opt = engine.optimize(side, FINE, init_hessian=freq)
        assert isinstance(opt, Evidence) and opt.task == "opt"
        payload = {"method": FINE.signature(), "molecule": side.fingerprint(),
                   "hessian": freq.hessian.sha256}
        old = jobs.store.key(Task(engine="nwchem", version_pin="7.2.3", kind="optimize",
                                  key_payload=payload, execution=site.execution))
        assert (opt.job_key != old) is reshaped
        attempt = jobs.store.attempt_dir(opt.job_key, 0)
        assert (attempt / "job.hess").read_text() == hess_text(written)
        decks.append((attempt / "job.nw").read_text())
    assert decks[0] == decks[1] == decks[2] and "trust 0.3\n  inhess 2" in decks[0]


def test_a_scan_point_holds_its_bond_and_starts_from_the_previous_vectors(nwchem, golden):
    """G2-P4: fixed_bond is a frozen zcoord bond (1-based atoms) and scf_guess the previous
    point's job.movecs, both in the job key only when given (every other key is unchanged); the
    start must carry the bond (NWChem moves an input to its zcoord value). An autoz failure
    continues in Cartesian coordinates with a spring restraint (r0 in bohr)."""
    jobs, site = nwchem
    engine, mol = NWChemEngine(jobs=jobs, site=site), _hcn_ts(golden)
    r = round(float(np.linalg.norm(mol.xyz.coords[2] - mol.xyz.coords[0])), 4)
    plain = engine.optimize(mol, FINE)
    point = engine.optimize(mol, FINE, fixed_bond=(2, 0, r), scf_guess=plain)
    assert isinstance(point, Evidence) and point.task == "opt"
    old = {"method": FINE.signature(), "molecule": mol.fingerprint(), "hessian": None}
    keys = [jobs.store.key(Task(engine="nwchem", version_pin="7.2.3", kind="optimize",
                                key_payload=p, execution=site.execution))
            for p in (old, {**old, "fixed_bond": [0, 2, r], "scf_guess": plain.job_key})]
    assert keys == [plain.job_key, point.job_key]
    first = jobs.store.attempt_dir(point.job_key, 0)
    deck = (first / "job.nw").read_text()
    assert f"\n  zcoord\n    bond 1 3 {r:.4f} rc constant\n  end\nend\ncharge 0" in deck
    assert "vectors input job.movecs" in deck and (first / "job.movecs").is_file()
    off = engine.optimize(mol, FINE, fixed_bond=(0, 2, r + 0.001))
    assert (off.kind, off.reason) == (FailureKind.INPUT_INVALID, "fixed_bond_mismatch")
    autoz = site.model_copy(update={"execution": site.execution.model_copy(
        update={"env": {"STUB_AUTOZ": "1"}})})
    held = NWChemEngine(jobs=jobs, site=autoz).optimize(mol, FINE, fixed_bond=(0, 2, r))
    assert isinstance(held, Evidence)
    decks = [(jobs.store.attempt_dir(held.job_key, i) / "job.nw").read_text() for i in range(3)]
    assert "zcoord" in decks[0] and "constraints" not in decks[0]
    spring = f"\n\nconstraints\n  spring bond 1 3 20.0 {r / BOHR_TO_ANGSTROM:.6f}\nend\n\ndft\n"
    for later in decks[1:]:  # the Cartesian start, then its timeout continuation
        assert "noautosym noautoz" in later and spring in later and "zcoord" not in later


def test_double_hybrids_are_rejected_and_autoz_falls_back_to_cartesians(nwchem, golden, tmp_path):
    """U0-P4: without dftmp2 NWChem drops the PT2 part silently. U9-P2: Cartesian coordinates
    are the autoz failure's continuation, not a site setting."""
    jobs, site = nwchem
    engine, mol = NWChemEngine(jobs=jobs, site=site), _hcn_ts(golden)
    for xc in ("B2PLYP", "b2gpplyp", "dsd-pbep86", "PWPB95"):
        failure = engine.energy(mol, FINE.model_copy(update={"functional": xc}))
        assert (failure.kind, failure.reason) == (FailureKind.INPUT_INVALID,
                                                  "unsupported_method:fine")
    assert jobs.stats().misses == 0  # no job was run
    task = Task(engine="nwchem", version_pin="7.2.3", kind="saddle", key_payload={},
                execution=site.execution, inputs={"mol": mol, "start": mol, "method": FINE})
    retry = engine.continuation(task, tmp_path, Failure(kind=FailureKind.INPUT_INVALID,
                                                        reason="autoz"))
    assert retry is not None and retry.execution == task.execution
    engine.prepare(retry, tmp_path)
    assert "noautosym noautoz" in (tmp_path / "job.nw").read_text()


def test_an_scf_failure_continues_with_cgmin_for_any_dft_and_old_vectors_for_wft(nwchem,
                                                                                tmp_path):
    """G2-P6: one plain SCF after cgmin prints <S2> (input._task; the P3c H3 bead, G34), so
    the rescue serves open-shell DFT too; WFT restarts from its old vectors."""
    jobs, site = nwchem
    engine, scf = NWChemEngine(jobs=jobs, site=site), Failure(kind=FailureKind.SCF_NOT_CONVERGED,
                                                               reason="scf")
    h3 = XYZ(["H", "H", "H"], np.array([[0.0, 0.0, 0.0], [0.93, 0.0, 0.0], [1.86, 0.0, 0.0]]))
    (tmp_path / "job.movecs").write_text("v")
    ccsd_t = MethodSpec(id="ccsd-t", kind="wft", wft_method="ccsd(t)", basis="def2-tzvpd")

    def retry(method: MethodSpec, charge: int, multiplicity: int) -> Task | None:
        mol = Molecule(h3, charge, multiplicity)
        task = Task(engine="nwchem", version_pin="7.2.3", kind="energy", key_payload={},
                    execution=site.execution, inputs={"mol": mol, "start": mol, "method": method})
        return engine.continuation(task, tmp_path, scf)

    closed, doublet, wft = retry(FINE, 1, 1), retry(FINE, 0, 2), retry(ccsd_t, 0, 2)
    assert closed.inputs["scf_rescue"] and doublet.inputs["scf_rescue"]  # H3+ and P3c's H3
    assert not wft.inputs["scf_rescue"] and wft.inputs["restart"] == [tmp_path / "job.movecs"]


def test_open_shell_ccsd_t_is_a_rohf_tce_job_and_all_electron_iodine_is_rejected(nwchem,
                                                                                 golden):
    jobs, site = nwchem
    engine = NWChemEngine(jobs=jobs, site=site)
    ccsd_t = MethodSpec(id="ccsd-t", kind="wft", wft_method="ccsd(t)", basis="def2-tzvpd")
    symbols, coords = geometry_block(golden.text("nwchem/G29/oh_ccsdt.out"))
    ev = engine.energy(Molecule(XYZ(list(symbols), coords), 0, 2), ccsd_t)
    assert isinstance(ev, Evidence) and ev.task == "sp" and ev.s2 is None
    assert (ev.level.method, ev.level.multiplicity) == ("ccsd(t)", 2)
    assert ev.energy_hartree == pytest.approx(-75.640597871999560, abs=1e-9)
    deck = (jobs.store.attempt_dir(ev.job_key, 0) / "job.nw").read_text()
    assert "  rohf\n  nopen 1\n" in deck and deck.rstrip().endswith("task tce energy")
    hi = Molecule(XYZ(["H", "I"], np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.61]])), 0, 1)
    failure = engine.energy(hi, ccsd_t.model_copy(update={"basis": "cc-pVDZ"}))
    assert (failure.kind, failure.reason) == (FailureKind.INPUT_INVALID,
                                              "no_ecp_for_basis:cc-pVDZ:I")
    assert (jobs.stats().hits, jobs.stats().misses) == (0, 1)
    keyed = engine.energy(hi, ccsd_t)  # the stub's DFT output fails it; only its key counts

    def key(mol: Molecule, **extra: int) -> str:
        payload = {"method": ccsd_t.signature(), "molecule": mol.fingerprint(), **extra}
        return jobs.store.key(Task(engine="nwchem", version_pin="7.2.3", kind="energy",
                                   key_payload=payload, execution=site.execution))

    oh = Molecule(XYZ(list(symbols), coords), 0, 2)
    assert ev.job_key == key(oh)  # G6-P3: an ECP atom's key names its freeze 4 (I 4s4p)
    assert keyed.job_key == key(hi, frozen_core=4) and "\n  freeze 4\n" in (
        jobs.store.attempt_dir(keyed.job_key, 0) / "job.nw").read_text()


def test_g01_two_vibrational_blocks_are_not_a_freq_evidence(nwchem, golden):
    symbols, coords = geometry_block(golden.text("nwchem/G01/nwchem.out"))
    mol = Molecule(XYZ(list(symbols), coords), 0, 1)
    failure = NWChemEngine(jobs=nwchem[0], site=nwchem[1]).frequencies(mol, XFINE)
    assert (failure.kind, failure.reason) == (FailureKind.INCOMPLETE_OUTPUT, "frequency_blocks:2")


def test_timeout_continues_from_the_latest_frame(nwchem, golden):
    jobs, site = nwchem
    mol = _hcn_ts(golden)
    ev = NWChemEngine(jobs=jobs, site=site).optimize(mol, FINE)
    assert isinstance(ev, Evidence) and ev.task == "opt"
    second = jobs.store.attempt_dir(ev.job_key, 1)
    assert (second / "job.drv.hess").is_file() and (second / "job.movecs").is_file()
    assert "vectors input job.movecs" in (second / "job.nw").read_text()
    final = read_xyz(jobs.store.run_dir / ev.final.file.path)
    assert np.allclose(final.coords, mol.xyz.coords + 0.01, atol=1e-7)
    assert ev.gradient == pytest.approx([1e-5, 0.0, 0.0] * 3)  # at the final frame
    assert read_xyz(jobs.store.run_dir / ev.start.file.path).coords == pytest.approx(mol.xyz.coords)


def test_saddle_shapes_the_hessian_along_the_mode_and_always_follows_mode_1(nwchem, golden):
    jobs, site = nwchem
    mol = _hcn_ts(golden)
    freq = NWChemEngine(jobs=jobs, site=site).frequencies(mol, FINE)
    saddle, mode = NWChemSaddle(jobs=jobs, site=site), np.eye(9)[8]  # not the TS mode
    ts = saddle.refine(mol, FINE, hessian=freq, mode=mode)
    assert isinstance(ts, Evidence) and ts.task == "saddle"
    first, resumed = (jobs.store.attempt_dir(ts.job_key, i) for i in (0, 1))
    assert "inhess 2\n  moddir 1" in (first / "job.nw").read_text()
    raw = np.load(jobs.store.run_dir / freq.hessian.path)
    shaped = read_hess(first / "job.hess", 3)
    assert np.allclose(shaped, shape_hessian(raw, mol.xyz.coords, mode), rtol=1e-8, atol=1e-12)
    values, vectors = np.linalg.eigh(shaped)  # G1-P1: -κ along the mode's internal part only
    rigid = external_basis(["H"] * 3, mol.xyz.coords)  # unweighted rigid motions
    internal = mode - rigid @ (rigid.T @ mode)
    assert np.count_nonzero(values < -1e-8) == 1 and values[0] <= -KAPPA_MIN
    assert abs(vectors[:, 0] @ internal) / np.linalg.norm(internal) == pytest.approx(1.0)
    payload = {"method": FINE.signature(), "molecule": mol.fingerprint(),
               "hessian": freq.hessian.sha256, "mode": sha256_text(json.dumps(mode.tolist()))}
    keys = [jobs.store.key(Task(engine="nwchem_saddle", version_pin="7.2.3", kind="saddle",
                                key_payload=p, execution=site.execution))
            for p in (payload, {**payload, "hessian_model": "negative_along_mode"})]
    assert keys[1] == ts.job_key != keys[0]  # the model is in the key: no older model reused
    deck = (resumed / "job.nw").read_text()  # timeout: the latest frame and driver Hessian
    assert "  moddir 1" in deck and "inhess" not in deck and not (resumed / "job.hess").exists()
    other = saddle.refine(mol, FINE, hessian=freq, mode=freq.imaginary_modes[0])
    assert other.job_key != ts.job_key  # the mode is part of the job key


def test_a_saddle_at_maxiter_fails_with_its_last_frame_and_takes_a_hessian_nearby(nwchem,
                                                                                   golden):
    """U6-P6: no continuation with the stalled driver Hessian; U6-P3: a freq within 0.5 Å."""
    jobs, site = nwchem
    mol = _hcn_ts(golden)
    freq = NWChemEngine(jobs=jobs, site=site).frequencies(mol, FINE)
    capped = site.model_copy(update={"execution": site.execution.model_copy(
        update={"env": {"STUB_MAXITER": "1"}})})
    near, far = (Molecule(XYZ(mol.xyz.symbols, mol.xyz.coords + [d, 0, 0]), 0, 1)
                 for d in (0.3, 0.6))
    failure = NWChemSaddle(jobs=jobs, site=capped).refine(near, FINE, hessian=freq,
                                                          mode=freq.imaginary_modes[0])
    assert (failure.kind, failure.reason) == (FailureKind.GEOMETRY_MAXITER, "maxiter")
    assert failure.energy_hartree == pytest.approx(-93.16699428)  # its last '@' step line
    last = read_xyz(jobs.store.run_dir / failure.final.file.path)
    assert np.allclose(last.coords, near.xyz.coords + 0.01, atol=1e-7)
    assert not jobs.store.attempt_dir(failure.job_key, 1).exists()  # one attempt only
    mismatch = NWChemSaddle(jobs=jobs, site=site).refine(far, FINE, hessian=freq,
                                                         mode=freq.imaginary_modes[0])
    assert (mismatch.kind, mismatch.reason) == (FailureKind.INPUT_INVALID,
                                                "hessian_geometry_mismatch")


def test_an_opt_autoz_failure_continues_from_its_latest_frame_in_cartesians(nwchem, golden,
                                                                           tmp_path):
    """G4-P4: VAL7 s10's QRC side stopped in autoz at frame 21 (gmax 4e-5), and its attempt_01
    restarted from the start (19 steps, 41.6 s). It continues from frame 21 in Cartesian
    coordinates with its vectors only (the driver Hessian is in internal coordinates): 1 step,
    3.3 s, the same minimum within 0.001 Å and 2e-8 Eh (r9 S1a probe). Before the first frame,
    and for a saddle, an autoz failure restarts from the same start."""
    jobs, site = nwchem
    engine, text = NWChemEngine(jobs=jobs, site=site), golden.text("nwchem/G33/opt_autoz.out")
    failure = classify_failure(text, returncode=236, timed_out=False)
    symbols, coords = geometry_block(text, 0)
    start, first = Molecule(XYZ(list(symbols), coords), 0, 1), tmp_path / "attempt_00"
    first.mkdir()
    for name in ("job.movecs", "job.drv.hess"):
        (first / name).write_text("x")

    def retry(kind: str) -> Task:
        task = Task(engine="nwchem", version_pin="7.2.3", kind=kind, key_payload={},
                    execution=site.execution, inputs={"mol": start, "start": start, "method": FINE,
                                                      "hessian": tmp_path / "side.npy"})
        return engine.continuation(task, first, failure)

    early = retry("optimize")
    assert early.inputs["cartesian"] and early.inputs["mol"] is start and "restart" not in early.inputs
    shutil.copyfile(golden.path("nwchem/G33/final-021.xyz"), first / "final-021.xyz")
    late, saddle = retry("optimize"), retry("saddle")
    assert saddle.inputs["cartesian"] and saddle.inputs["mol"] is start
    assert late.inputs["cartesian"] and late.inputs["restart"] == [first / "job.movecs"]
    assert late.inputs["hessian"] is None and late.inputs["start"] is start
    assert np.allclose(late.inputs["mol"].xyz.coords, read_xyz(first / "final-021.xyz").coords)
    second = tmp_path / "attempt_01"
    second.mkdir()
    engine.prepare(late, second)
    deck = (second / "job.nw").read_text()
    assert "noautosym noautoz" in deck and "vectors input job.movecs" in deck
    assert "inhess" not in deck and not (second / "job.drv.hess").exists()


def test_g13_string_profile_from_its_initial_path(nwchem, golden):
    text = golden.text("nwchem/G13/nwchem_string.out")
    (s0, c0), (s1, c1) = geometry_block(text, 0), geometry_block(text, 1)
    start, end = Molecule(XYZ(list(s0), c0), 0, 1), Molecule(XYZ(list(s1), c1), 0, 1)
    jobs, site = nwchem
    initial = write_xyz_trajectory([start.xyz, end.xyz], jobs.store.run_dir / "initial.xyz")
    ref = FileRef(path="initial.xyz", sha256=sha256_file(initial))
    profile = NWChemString(jobs=jobs, site=site).find_path(start, end, XFINE, images=11,
                                                           initial_path=ref)
    assert isinstance(profile, PathProfile) and len(profile.energies_hartree) == 11
    assert profile.ts is None
    attempt = (jobs.store.run_dir / profile.images.path).parent
    assert "xyz_path initial_path.xyz" in (attempt / "job.nw").read_text()
    assert sha256_file(attempt / "initial_path.xyz") == ref.sha256

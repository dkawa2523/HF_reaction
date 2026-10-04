"""sp targets (X5): the points thermo.reaction_points reads, per reaction outcome; each from its
freq's SCF (X3), CCSD(T) only on spin-clean freqs; a failure on an auxiliary point only is no
artifact."""

from pathlib import Path

import numpy as np
import pytest
from fakes import FakeQM, double_well, fake_level, write_geometry

from hfauto.backends.protocols import Capability
from hfauto.chemistry.topology import state_label
from hfauto.chemistry.xyz import geometry_fingerprint
from hfauto.core import records as R
from hfauto.core.evidence import Evidence, Failure, FailureKind
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec
from hfauto.core.system import CompositionInput, SpeciesInput, SystemConfig
from hfauto.drivers.minimum import calc_id
from hfauto.stages.single_point import SinglePointConfig, SinglePointStage

T, Out = R.ArtifactType, R.CaseOutcome
DFT, BIG = (MethodSpec(id=b, kind="dft", functional="xfake", basis=b) for b in ("svp", "tzvp"))
# minimum -> (symbols, composition_id, state_label, tier). The double well HNO is composition
# c = NH + O of the system (its reactant's fragments); c2 = NH + F (Hill FHN) is another one. r2
# is a second conformer of the reactant state, x another state of HNO.
MONOMER = {k: (s, (1 + i) * np.eye(3)[: len(s)]) for i, (k, s) in enumerate(
    (("nh", "NH"), ("o", "O"), ("f", "F")))}
LABEL = {k: state_label(list(s), x) for k, (s, x) in MONOMER.items()}
MINIMA = {"r": ("NHO", "HNO", "r", "dft"), "r2": ("NHO", "HNO", "r", "dft"),
          "p": ("NHO", "HNO", "p", "dft"), "x": ("NHO", "HNO", "x", "dft"),
          "r_screen": ("NHO", "HNO", "r", "screen"), "nh": ("NH", "nh", LABEL["nh"], "dft"),
          "o": ("O", "o", LABEL["o"], "dft"), "f": ("F", "f", LABEL["f"], "dft")}
REFERENCE = {"m_r", "m_r2", "m_nh", "m_o"}  # the reactant state and its separated NH + O


class SpQM(FakeQM):
    """Single points on any molecule (the monomers are off the double well) at 0 Eh; the key
    names scf_guess, as NWChem's does for DFT. ``fail``: geometries whose single point fails."""

    fail: frozenset[str] = frozenset()

    def energy(self, mol, method, *, scf_guess=None):
        self.guesses.append(scf_guess)
        key = self._key("sp", mol.fingerprint(), method.signature(), scf_guess and scf_guess.job_key)
        if geometry_fingerprint(mol.xyz.symbols, mol.xyz.coords) in self.fail:
            return Failure(kind=FailureKind.SCF_NOT_CONVERGED, reason="scf")
        start = self._start(mol, key)
        return Evidence(engine=self.name, task="sp", level=fake_level(method, mol), start=start,
                        final=start, energy_hartree=0.0, output=start.file, job_key=key)


def _inputs(tmp_run: Path, outcome: Out | None) -> tuple[Manifest, dict[str, Evidence], str]:
    """One reaction with ``outcome`` (degenerate: r -> r; barrierless: no saddle claim); m_p and
    m_twin (another record of the product basin) share one freq calculation."""
    qm = FakeQM(tmp_run, pes := double_well())
    point = {"r": "reactant", "p": "product"}
    freq = {k: qm.frequencies(pes.molecule(point[k]), DFT) for k in point}
    ts = qm.frequencies(pes.molecule("ts"), DFT)
    for i, (k, (symbols, *_)) in enumerate(MINIMA.items()):
        if k not in freq:
            x = MONOMER[k][1] if k in MONOMER else (1 + i) * np.eye(3)[: len(symbols)]
            geom = write_geometry(tmp_run, f"{k}.xyz", list(symbols), x)
            external = {1: 3, 2: 5}.get(len(symbols), 6)  # an atom, a diatomic, the rest
            freq[k] = freq["r"].model_copy(update={
                "start": geom, "final": geom, "job_key": k, "n_external": external,
                "frequencies_cm1": (1000.0,) * (3 * len(symbols) - external)})
    arts = [Artifact(artifact_id=calc_id(e), type=T.CALCULATION, payload=e)
            for e in (*freq.values(), ts)]
    arts += [Artifact(artifact_id=f"m_{k}", type=T.MINIMUM, payload=R.MinimumRecord(
        minimum_id=f"m_{k}", basin_id=k, composition_id=comp, species_id=k, tier=tier,
        level_key="x", opt_calc="o", freq_calc=calc_id(freq[k]), energy_hartree=0.0,
        state_label=label)) for k, (_, comp, label, tier) in MINIMA.items()]
    twin = next(a for a in arts if a.artifact_id == "m_p")
    arts.append(twin.model_copy(update={"artifact_id": "m_twin", "payload": twin.payload.model_copy(
        update={"minimum_id": "m_twin"})}))
    arts += [Artifact(artifact_id=k, type=T.SPECIES, payload=R.SpeciesRecord(
        species_id=k, composition_id=k, charge=0, multiplicity=1, geometry=freq[k].final,
        source="input", state_label=LABEL[k])) for k in MONOMER]
    term = (R.StoichTerm(composition_id="HNO", coefficient=1),)
    saddle = R.SaddleClaim(saddle_calc="s", freq_calc=calc_id(ts), imag_cm1=-900.0,
                           energy_hartree=0.0)
    rx = R.ReactionRecord(
        reaction_id="rx", reactants=term, products=term, endpoints=("a", "b"), source="declared",
        minima=("m_r", "m_r") if outcome is Out.DEGENERATE else ("m_r", "m_p"),
        outcome=outcome, saddle=None if outcome is Out.BARRIERLESS else saddle)
    arts.append(Artifact(artifact_id="rx", type=T.REACTION, payload=rx))
    return (Manifest(run_id="r", stage_id="v", created_at="t", artifacts=arts),
            {f"m_{k}": e for k, e in freq.items()} | {"m_twin": freq["p"], calc_id(ts): ts},
            calc_id(ts))


def _run(fake_runtime, tmp_run, inputs, methods, fail=frozenset()):
    system = SystemConfig(
        system_id="s", species=[SpeciesInput(id=k, xyz=Path(f"{k}.xyz"), multiplicity=1)
                                for k in MONOMER],
        compositions=[CompositionInput(id="c", components={"nh": 1, "o": 1}),
                      CompositionInput(id="c2", components={"nh": 1, "f": 1})])
    qm = SpQM(tmp_run, double_well())
    qm.guesses, qm.fail = [], fail
    rt = fake_runtime(system, {(Capability.QM, "nwchem"): qm}, methods=methods)
    config = SinglePointConfig(engine="nwchem", methods=list(methods))
    return SinglePointStage().run(inputs, config, rt), qm.guesses


@pytest.mark.parametrize(("outcome", "expected"), [
    (Out.ELEMENTARY_STEP, REFERENCE | {"m_p", "m_twin", "TS"}),
    (Out.REASSIGNED, REFERENCE | {"m_p", "m_twin", "TS"}),
    (Out.DEGENERATE, REFERENCE | {"TS"}),
    (Out.BARRIERLESS, REFERENCE | {"m_p", "m_twin"}),
    *((o, set()) for o in (Out.MULTI_STEP, Out.SAME_BASIN, Out.OUT_OF_WINDOW, Out.UNRESOLVED,
                           Out.BLOCKED, None)),
])
def test_sp_computes_only_the_points_the_reactions_read(fake_runtime, tmp_run, outcome, expected):
    inputs, freq, ts = _inputs(tmp_run, outcome)
    out, guesses = _run(fake_runtime, tmp_run, inputs, {"tzvp": BIG})
    assert {"TS" if p == ts else p for a in out for p in a.parents} == expected
    for a in out:  # each on the final geometry of its subject's freq calculation
        assert {freq[p].final.fingerprint for p in a.parents} == {a.payload.start.fingerprint}
    assert all(guess is not None for guess in guesses)  # X3: each from its freq's SCF
    shared = [set(a.parents) for a in out if {"m_p", "m_twin"} & set(a.parents)]
    assert shared in ([], [{"m_p", "m_twin"}])  # one freq calculation, one sp


def test_a_failure_on_an_auxiliary_point_only_is_no_artifact(fake_runtime, tmp_run):
    """A barrierless step reads its own ends (dG_rxn); its separated NH + O give only dG_assoc:
    their failed single points are left out, its product's stays a failed artifact."""
    inputs, freq, _ = _inputs(tmp_run, Out.BARRIERLESS)
    fail = frozenset(freq[m].final.fingerprint for m in ("m_nh", "m_o", "m_p"))
    out, _ = _run(fake_runtime, tmp_run, inputs, {"tzvp": BIG}, fail)
    failed = {p for a in out if a.status == "failed" for p in a.parents}
    assert failed == {"m_p", "m_twin"} and {p for a in out for p in a.parents} == {
        "m_r", "m_r2", "m_p", "m_twin"}
    connected, _ = _run(fake_runtime, tmp_run, _inputs(tmp_run, Out.ELEMENTARY_STEP)[0],
                        {"tzvp": BIG}, fail)  # with a TS, R_sep is a G_ref candidate
    assert {"m_nh", "m_o"} <= {p for a in connected if a.status == "failed" for p in a.parents}


def test_ccsd_t_runs_only_where_the_freq_passes_spin_ok(fake_runtime, tmp_run):
    """X3-4: a broken-symmetry point (the product's freq <S2> 1.0 in a singlet) gets the DFT
    layer, from its freq's SCF, but no CCSD(T)."""
    inputs, freq, _ = _inputs(tmp_run, Out.ELEMENTARY_STEP)
    hot = freq["m_p"].model_copy(update={"s2": 1.0})
    inputs = inputs.model_copy(update={"artifacts": [
        a.model_copy(update={"payload": hot}) if a.payload == freq["m_p"] else a
        for a in inputs.artifacts]})
    ccsd_t = MethodSpec(id="ccsd-t", kind="wft", wft_method="ccsd(t)", basis="tzvp")
    out, guesses = _run(fake_runtime, tmp_run, inputs, {"tzvp": BIG, "ccsd-t": ccsd_t})
    by_level = {(a.payload.level.method, p) for a in out for p in a.parents}
    assert ("xfake", "m_p") in by_level and ("ccsd(t)", "m_p") not in by_level
    assert ("ccsd(t)", "m_r") in by_level and hot in guesses

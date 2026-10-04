"""sp targets (review U8-P4): only the points the ranking reads, per reaction outcome."""

from pathlib import Path

import numpy as np
import pytest
from fakes import FakeQM, double_well, fake_level, write_geometry

from hfauto.backends.protocols import Capability
from hfauto.core import records as R
from hfauto.core.evidence import Evidence
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec
from hfauto.core.system import CompositionInput, SpeciesInput, SystemConfig
from hfauto.drivers.minimum import calc_id
from hfauto.stages.single_point import SinglePointConfig, SinglePointStage

T, Out = R.ArtifactType, R.CaseOutcome
DFT, BIG = (MethodSpec(id=b, kind="dft", functional="xfake", basis=b) for b in ("svp", "tzvp"))
# minimum -> (symbols, composition_id, state_label, tier). The double well HNO is composition
# c = NH + O of the system; c2 = NH + F (Hill FHN) is another one. r2 is a second conformer of
# the reactant state, x another state of HNO.
MINIMA = {"r": ("NHO", "HNO", "r", "dft"), "r2": ("NHO", "HNO", "r", "dft"),
          "p": ("NHO", "HNO", "p", "dft"), "x": ("NHO", "HNO", "x", "dft"),
          "r_screen": ("NHO", "HNO", "r", "screen"), "nh": ("NH", "nh", "nh", "dft"),
          "o": ("O", "o", "o", "dft"), "f": ("F", "f", "f", "dft")}
REFERENCE = {"m_r", "m_r2", "m_nh", "m_o"}  # the reactant state and the monomers of c


class SpQM(FakeQM):
    """Single points on any molecule (the monomers are off the double well) at 0 Eh."""

    def energy(self, mol, method, *, scf_guess=None):
        key = self._key("sp", mol.fingerprint(), method.signature())
        start = self._start(mol, key)
        return Evidence(engine=self.name, task="sp", level=fake_level(method, mol), start=start,
                        final=start, energy_hartree=0.0, output=start.file, job_key=key)


def _inputs(tmp_run: Path, outcome: Out | None) -> tuple[Manifest, dict[str, Evidence], str]:
    """One reaction with ``outcome`` (degenerate: r -> r; barrierless: no saddle claim)."""
    qm = FakeQM(tmp_run, pes := double_well())
    point = {"r": "reactant", "p": "product"}
    freq = {k: qm.frequencies(pes.molecule(point[k]), DFT) for k in point}
    ts = qm.frequencies(pes.molecule("ts"), DFT)
    for i, (k, (symbols, *_)) in enumerate(MINIMA.items()):
        if k not in freq:
            geom = write_geometry(tmp_run, f"{k}.xyz", list(symbols),
                                  (1 + i) * np.eye(3)[: len(symbols)])
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
    arts += [Artifact(artifact_id=k, type=T.SPECIES, payload=R.SpeciesRecord(
        species_id=k, composition_id=k, charge=0, multiplicity=1, geometry=freq[k].final,
        source="input", state_label=k)) for k in ("nh", "o", "f")]
    term = (R.StoichTerm(composition_id="HNO", coefficient=1),)
    saddle = R.SaddleClaim(saddle_calc="s", freq_calc=calc_id(ts), imag_cm1=-900.0,
                           energy_hartree=0.0)
    rx = R.ReactionRecord(
        reaction_id="rx", reactants=term, products=term, endpoints=("a", "b"), source="declared",
        minima=("m_r", "m_r") if outcome is Out.DEGENERATE else ("m_r", "m_p"),
        outcome=outcome, saddle=None if outcome is Out.BARRIERLESS else saddle)
    arts.append(Artifact(artifact_id="rx", type=T.REACTION, payload=rx))
    return (Manifest(run_id="r", stage_id="v", created_at="t", artifacts=arts),
            {f"m_{k}": e for k, e in freq.items()} | {calc_id(ts): ts}, calc_id(ts))


@pytest.mark.parametrize(("outcome", "expected"), [
    (Out.ELEMENTARY_STEP, REFERENCE | {"m_p", "TS"}),
    (Out.REASSIGNED, REFERENCE | {"m_p", "TS"}),
    (Out.DEGENERATE, REFERENCE | {"TS"}),
    (Out.BARRIERLESS, REFERENCE | {"m_p"}),
    *((o, set()) for o in (Out.MULTI_STEP, Out.SAME_BASIN, Out.OUT_OF_WINDOW, Out.UNRESOLVED,
                           Out.BLOCKED, None)),
])
def test_sp_computes_only_the_points_the_ranking_reads(fake_runtime, tmp_run, outcome, expected):
    inputs, freq, ts = _inputs(tmp_run, outcome)
    system = SystemConfig(
        system_id="s", species=[SpeciesInput(id=k, xyz=Path(f"{k}.xyz"), multiplicity=1)
                                for k in ("nh", "o", "f")],
        compositions=[CompositionInput(id="c", components={"nh": 1, "o": 1}),
                      CompositionInput(id="c2", components={"nh": 1, "f": 1})])
    rt = fake_runtime(system, {(Capability.QM, "nwchem"): SpQM(tmp_run, double_well())},
                      methods={"tzvp": BIG})
    out = SinglePointStage().run(inputs, SinglePointConfig(engine="nwchem", methods=["tzvp"]), rt)
    assert {"TS" if p == ts else p for a in out for p in a.parents} == expected
    for a in out:  # each on the final geometry of its subject's freq calculation
        assert {freq[p].final.fingerprint for p in a.parents} == {a.payload.start.fingerprint}

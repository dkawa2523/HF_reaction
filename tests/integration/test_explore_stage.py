"""explore stage with FakeDiscovery (§8.2): AFIR only after NT2 found no maximum, the kcal/mol
window, failed attempts, one species per product basin (degenerate products kept apart from the
source), relaxation discovery (legacy relaxation K case), two lowest sources per state, the
screen-optimized start structure and the source's charge and multiplicity (HCN⁺•: dissociations
are trials), and parallel units recorded as serial ones."""

import dataclasses

import numpy as np
from fakes import NH3_HF, NH3_HF_EXCHANGED, NH3_HF_SYMBOLS, FakeDiscovery, write_geometry

from hfauto.backends.protocols import NO_NT2_MAXIMUM, Capability, DiscoveryResult
from hfauto.chemistry.topology import state_label
from hfauto.chemistry.trials import perturb_linear
from hfauto.chemistry.xyz import XYZ
from hfauto.core.evidence import Evidence, Failure, FailureKind, Level
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec
from hfauto.core.records import ArtifactType, MinimumRecord, SpeciesRecord
from hfauto.core.system import SpeciesInput, SystemConfig
from hfauto.stages.explore import ExploreConfig, ExploreStage, _Basin, _basin_of

SYMBOLS, SPECIES, MINIMUM = ["H", "C", "N"], ArtifactType.SPECIES, ArtifactType.MINIMUM
HCN = np.array([[0, 0, -1.066], [0, 0, 0], [0, 0, 1.156]])
HNC = np.array([[0, 0, 2.17], [0, 0, 0], [0, 0, 1.17]])
LABEL = state_label(SYMBOLS, HCN)
LEVEL = Level(program="fake", version="0", method="gfn2", charge=1, multiplicity=2)
CATION = "CHN_q1_m2"
CH, CN = ((0, 1),), ((1, 2),)  # the two dissociation drives of HCN⁺•


def _optimized(i):  # the screen-optimized HCN of m{i}: the input displaced by 0.02 (i + 1) Å
    return HCN + 0.02 * (i + 1)


def _inputs(root):
    """Three HCN⁺• screen minima of one state (each with its opt calculation); m0 also holds a
    seed that relaxed from HNC. m3 is an Ar atom: never a source."""
    geo, ar = (write_geometry(root, f"in/{n}.xyz", s, x) for n, s, x in (
        ("hcn", SYMBOLS, HCN), ("ar", ["Ar"], np.zeros((1, 3)))))
    out = [Artifact(artifact_id=sid, type=SPECIES, payload=SpeciesRecord(
        species_id=sid, composition_id=comp, charge=q, multiplicity=q + 1, geometry=g,
        source="conformer", state_label=state_label(g.symbols, x)))
        for sid, comp, q, g, x in (
            ("seed", CATION, 1, geo, HNC), ("s0", CATION, 1, geo, HCN), ("s1", CATION, 1, geo, HCN),
            ("s2", CATION, 1, geo, HCN), ("s3", "Ar_q0_m1", 0, ar, np.zeros((1, 3))))]
    for i in range(4):
        final = ar if i == 3 else write_geometry(root, f"opt/m{i}.xyz", SYMBOLS, _optimized(i))
        out.append(Artifact(artifact_id=f"opt{i}", type=ArtifactType.CALCULATION, payload=Evidence(
            engine="fake", task="opt", level=LEVEL, start=final if i == 3 else geo, final=final,
            energy_hartree=i / 10, output=final.file, job_key=f"opt{i}")))
    return Manifest(run_id="r", stage_id="screen", created_at="t", artifacts=out + [
        Artifact(artifact_id=f"m{i}", type=MINIMUM, payload=MinimumRecord(
            minimum_id=f"m{i}", basin_id=f"b{i}", species_id=f"s{i}", tier="screen",
            composition_id="Ar_q0_m1" if i == 3 else CATION, level_key="k",
            opt_calc=f"opt{i}", freq_calc="f", energy_hartree=i / 10,
            state_label=state_label(["Ar"], np.zeros((1, 3))) if i == 3 else LABEL,
            members=(f"s{i}", "seed")[: 2 - min(i, 1)]))
        for i in range(4)])


def _explore(fake_runtime, tmp_run):
    """(run(rt) -> artifacts, the runtime, the fake): both sources reach one HNC basin by the
    1,2-H shift; C–N cuts test the window edges, C–H cuts the AFIR gating."""
    hnc = write_geometry(tmp_run, "fake/hnc.xyz", SYMBOLS, HNC)

    def script(source, trial, method, settings):
        if trial.mechanism == "afir":
            return Failure(kind=FailureKind.NONZERO_EXIT, reason="scripted")
        m0 = trial.source_minimum == "m0"
        act, reason = {("transfer", CH): (30.0, None), ("dissociation", CN): (
            50.0 if m0 else 50.01, None), ("dissociation", CH): (
            None if m0 else 20.0, NO_NT2_MAXIMUM if m0 else "same_as_source")}[
            (trial.kind, trial.dissociations)]
        found = reason is None
        return DiscoveryResult(
            outcome="product" if found else "negative", reason=reason,
            product=hnc if found else None, ts=hnc if act else None, ts_imag_cm1=-1400.0,
            dE_act_kcal=act, dE_rxn_kcal=40.0 if found else None,
            irc_connected_to_source=found, electronic_temperature_K=300.0, job_key="j")

    fake = FakeDiscovery(script)
    rt = fake_runtime(SystemConfig(system_id="t", species=[SpeciesInput(id="hcn", smiles="C#N")]),
                      {(Capability.DISCOVERY, "readuct"): fake},
                      methods={"gfn2": MethodSpec(id="gfn2", kind="xtb", gfn=2)})
    config = ExploreConfig(engine="readuct", method="gfn2")
    return lambda runtime: ExploreStage().run(_inputs(tmp_run), config, runtime), rt, fake


def test_explore_records_every_attempt(fake_runtime, tmp_run):
    run, rt, fake = _explore(fake_runtime, tmp_run)
    out = run(rt)
    found = {(a.payload.source_minimum, a.payload.trial and a.payload.trial.kind,
              a.payload.trial and a.payload.trial.dissociations, a.payload.mechanism): a
             for a in out if a.type == ArtifactType.DISCOVERY}
    got = {key: (a.payload.outcome, a.payload.reason) for key, a in found.items()}
    [species] = [a.payload for a in out if a.type == SPECIES]  # one HNC basin
    assert (species.source, species.state_label) == ("discovery", state_label(SYMBOLS, HNC))
    assert species.composition_id == CATION
    for source in ("m0", "m1"):
        shift = found[(source, "transfer", CH, "nt2")].payload
        assert (shift.outcome, shift.product_species) == ("product", species.species_id)
    assert got[("m0", "dissociation", CN, "nt2")] == ("product", None)  # 50 and 40 kcal/mol
    assert got[("m1", "dissociation", CN, "nt2")] == ("negative", "out_of_window")  # 50.01
    assert got[("m0", "dissociation", CH, "nt2")] == ("negative", NO_NT2_MAXIMUM)
    assert found[("m0", "dissociation", CH, "afir")].failure.kind == FailureKind.NONZERO_EXIT
    negative = found[("m1", "dissociation", CH, "nt2")].payload  # a TS: no AFIR, barrier kept
    assert (negative.reason, negative.dE_act_kcal) == ("same_as_source", 20.0)
    assert [c[1].source_minimum for c in fake.calls if c[1].mechanism == "afir"] == ["m0"]
    assert got[("m0", None, None, "relaxation")] == ("negative", f"collapsed_to:{LABEL}")
    assert {c[1].source_minimum for c in fake.calls} == {"m0", "m1"}  # m2: third lowest; m3: Ar
    for source, trial, *_ in fake.calls:  # from the opt final (bent: linear), not species.geometry
        i = int(trial.source_minimum[1:])
        assert np.allclose(source.xyz.coords, perturb_linear(_optimized(i)))


def test_parallel_units_are_recorded_as_serial_ones(fake_runtime, tmp_run):
    """Artifacts in input order and the same engine calls (the JobStore keys) on 1 or 4 cores."""
    run, rt, fake = _explore(fake_runtime, tmp_run)
    parallel, calls = run(rt), list(fake.calls)
    fake.calls.clear()
    serial = run(dataclasses.replace(rt, site=rt.site.model_copy(update={"cores": 1})))
    assert rt.site.cores == 4 and [a.model_dump() for a in parallel] == [
        a.model_dump() for a in serial]
    assert sorted(repr((m.fingerprint(), t)) for m, t, *_ in calls) == sorted(
        repr((m.fingerprint(), t)) for m, t, *_ in fake.calls)


def test_a_degenerate_product_never_joins_its_source_basin():
    """The NH3·HF double H exchange is the source's basin by permutation-invariant RMSD; its
    atom-indexed bonds keep it apart, and a second trial that reaches it joins it."""
    label, xyz = state_label(NH3_HF_SYMBOLS, NH3_HF), XYZ(list(NH3_HF_SYMBOLS), NH3_HF)
    source = _Basin("source", label, xyz)
    first = _Basin("p1", label, XYZ(list(NH3_HF_SYMBOLS), NH3_HF_EXCHANGED))
    second = _Basin("p2", label, XYZ(list(NH3_HF_SYMBOLS), NH3_HF_EXCHANGED + 0.01))
    assert _basin_of(first, [source], degenerate=False) == "source"
    assert _basin_of(first, [source], degenerate=True) is None
    assert _basin_of(second, [source, first], degenerate=True) == "p1"

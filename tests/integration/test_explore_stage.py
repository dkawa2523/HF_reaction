"""explore stage with FakeDiscovery (§8.2): one NT2 attempt per unit, the kcal/mol window,
failed attempts, one species per product basin as labelled (a relabelled product, degenerate
ones included, kept apart), each discovery's source species (X2), one relaxation product per
lost seed state (R6: a resolved bond change from its basin, X3), two lowest sources per state,
the screen-optimized start structure and the source's charge and multiplicity (HCN⁺•:
dissociations are trials), and parallel units recorded as serial ones."""

import dataclasses

import numpy as np
from fakes import (
    NH3_HF,
    NH3_HF_EXCHANGED,
    NH3_HF_SYMBOLS,
    FakeDiscovery,
    write_geometry,
    xyz_loader,
)

from hfauto.backends.protocols import Capability, DiscoveryResult
from hfauto.chemistry.elements import covalent_radius
from hfauto.chemistry.topology import BOND_TOLERANCE_A, state_label
from hfauto.chemistry.trials import perturb_linear
from hfauto.chemistry.xyz import XYZ
from hfauto.core.evidence import Evidence, Failure, FailureKind, Level
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec
from hfauto.core.records import ArtifactType, MinimumRecord, SpeciesRecord
from hfauto.core.system import SpeciesInput, SystemConfig
from hfauto.stages.explore import ExploreConfig, ExploreStage, _Basin, _basin_of, _relaxations

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
    geo, hnc, ar = (write_geometry(root, f"in/{n}.xyz", s, x) for n, s, x in (
        ("hcn", SYMBOLS, HCN), ("hnc", SYMBOLS, HNC), ("ar", ["Ar"], np.zeros((1, 3)))))
    out = [Artifact(artifact_id=sid, type=SPECIES, payload=SpeciesRecord(
        species_id=sid, composition_id=comp, charge=q, multiplicity=q + 1, geometry=g,
        source="conformer", state_label=state_label(g.symbols, x)))
        for sid, comp, q, g, x in (
            ("seed", CATION, 1, hnc, HNC), ("s0", CATION, 1, geo, HCN), ("s1", CATION, 1, geo, HCN),
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
    1,2-H shift; C–N cuts test the window edges, C–H cuts a failure and a kept barrier."""
    hnc = write_geometry(tmp_run, "fake/hnc.xyz", SYMBOLS, HNC)

    def script(source, trial, method, settings):
        m0 = trial.source_minimum == "m0"
        if m0 and (trial.kind, trial.dissociations) == ("dissociation", CH):
            return Failure(kind=FailureKind.NONZERO_EXIT, reason="scripted")
        act, reason = {("transfer", CH): (30.0, None), ("dissociation", CN): (
            50.0 if m0 else 50.01, None), ("dissociation", CH): (20.0, "same_as_source")}[
            (trial.kind, trial.dissociations)]
        found = reason is None
        return DiscoveryResult(
            outcome="product" if found else "negative", reason=reason,
            product=hnc if found else None, ts=hnc if act else None, ts_imag_cm1=-1400.0,
            dE_act_kcal=act, dE_rxn_kcal=40.0 if found else None,
            irc_connected_to_source=found, electronic_temperature_K=300.0)

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
    for source in ("m0", "m1"):  # from the source's representative, in whose order it ran
        shift = found[(source, "transfer", CH, "nt2")].payload
        assert (shift.outcome, shift.product_species) == ("product", species.species_id)
        assert shift.source_species == "s" + source[1:]
    assert got[("m0", "dissociation", CN, "nt2")] == ("product", None)  # 50 and 40 kcal/mol
    assert got[("m1", "dissociation", CN, "nt2")] == ("negative", "out_of_window")  # 50.01
    assert found[("m0", "dissociation", CH, "nt2")].failure.kind == FailureKind.NONZERO_EXIT
    negative = found[("m1", "dissociation", CH, "nt2")].payload  # a TS: the barrier is kept
    assert (negative.reason, negative.dE_act_kcal) == ("same_as_source", 20.0)
    assert len(fake.calls) == len(found) - 1  # one NT2 call per unit; the relaxation has none
    relaxed = found[("m0", None, None, "relaxation")].payload  # the lost HNC seed, unrelaxed
    assert (relaxed.outcome, relaxed.product_species, relaxed.ts) == ("product", "seed", None)
    assert relaxed.source_species == "spc_relax_seed"  # the seed refined as its own species
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


def _nh3_hf(nh):
    """C3v NH3·HF with N···H at r_thr + nh (Å) and H–F 0.95 Å."""
    r = covalent_radius("N") + covalent_radius("H") + BOND_TOLERANCE_A + nh
    return np.array([[0.0, 0.0, 0.0], [0.94, 0.0, -0.38], [-0.47, 0.814, -0.38],
                     [-0.47, -0.814, -0.38], [0.0, 0.0, r], [0.0, 0.0, r + 0.95]])


def test_one_relaxation_product_per_lost_seed_state(tmp_run):
    """R6 (X3): NH3 + HF seeds collapsed into N–H-bonded basins; no screen minimum keeps their
    state. It is lost through the seed of lowest low-level energy that a resolved bond change
    separates from its basin: 'near' (lowest) crosses the threshold by ±0.006 Å as S19's c01
    did, 'far' (+2.00 / −0.35 Å, S6's H–O) leaves from its own basin m1, 'late' has no energy;
    'deep' holds the kept state. Without m1 nothing is lost."""
    symbols = list(NH3_HF_SYMBOLS)
    final = {"m0": XYZ(symbols, _nh3_hf(-0.006)), "m1": XYZ(symbols, _nh3_hf(-0.35))}
    species = {sid: SpeciesRecord(
        species_id=sid, composition_id="FH4N_q0_m1", charge=0, multiplicity=1,
        geometry=write_geometry(tmp_run, f"in/{sid}.xyz", symbols, _nh3_hf(nh)),
        source="conformer", state_label=state_label(symbols, _nh3_hf(nh)),
        energy_hartree=energy)
        for sid, nh, energy in (("near", 0.006, -2.0), ("far", 2.0, -1.0), ("late", 2.0, None),
                                ("deep", -0.35, None))}
    minima = [MinimumRecord(minimum_id=mid, basin_id=mid, species_id=members[0], tier="screen",
                            composition_id="FH4N_q0_m1", level_key="k", opt_calc="o",
                            freq_calc="f", energy_hartree=0.0,
                            state_label=state_label(symbols, final[mid].coords), members=members)
              for mid, members in (("m0", ("deep", "near")), ("m1", ("far", "late")))]
    [record] = _relaxations(species, minima, final, xyz_loader(tmp_run))
    assert (record.source_minimum, record.product_species) == ("m1", "far")
    assert (record.mechanism, record.outcome) == ("relaxation", "product")
    assert _relaxations(species, minima[:1], final, xyz_loader(tmp_run)) == []


def test_a_product_joins_a_known_basin_only_as_labelled():
    """X2: the NH3·HF double H exchange is the source's basin by permutation-invariant RMSD, but
    with other atom-indexed bonds: it stays a species of its own (the discovery's ends keep its
    labelling), which a second trial that reaches it joins; the source's own labelling joins the
    source."""
    def basin(sid, coords):
        return _Basin(sid, XYZ(list(NH3_HF_SYMBOLS), coords))

    source, first = basin("source", NH3_HF), basin("p1", NH3_HF_EXCHANGED)
    assert _basin_of(first, [source]) is None
    assert _basin_of(basin("p2", NH3_HF_EXCHANGED + 0.01), [source, first]) == "p1"
    assert _basin_of(basin("p3", NH3_HF + 0.01), [source, first]) == "source"

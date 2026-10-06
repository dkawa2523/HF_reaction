"""explore stage with FakeDiscovery and scripted classes (chemistry.trials is tested apart):
edges over generations under one count budget, attempts that reach a known TS (S6: one edge),
edges without a TS (an attempt's start that relaxed, a seed state the screen lost: R6), ends
identified as labelled, unconnected edges and edges reached through their product end, and
parallel units recorded as serial ones."""

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
from hfauto.chemistry import trials
from hfauto.chemistry.elements import covalent_radius
from hfauto.chemistry.topology import BOND_TOLERANCE_A, state_label
from hfauto.chemistry.xyz import XYZ
from hfauto.core.evidence import Evidence, Failure, FailureKind, Level
from hfauto.core.hashing import sha256_text
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec
from hfauto.core.records import ArtifactType, MinimumRecord, SpeciesRecord
from hfauto.core.system import SpeciesInput, SystemConfig
from hfauto.stages.explore import ExploreConfig, ExploreStage, _Network, _relaxations

SPECIES, MINIMUM, DISCOVERY = ArtifactType.SPECIES, ArtifactType.MINIMUM, ArtifactType.DISCOVERY
HCN_SYMBOLS = ["H", "C", "N"]
HCN = np.array([[0, 0, -1.066], [0, 0, 0], [0, 0, 1.156]])
HNC = np.array([[0, 0, 2.17], [0, 0, 0], [0, 0, 1.17]])
H_CN = np.array([[0, 0, -3.0], [0, 0, 0], [0, 0, 1.156]])  # H apart: another state
NH4_F = np.vstack([NH3_HF[:4], [[0.31, -0.98, 0.0]], NH3_HF[5:]])  # the HF proton on N
NH2_H_FH = NH3_HF + np.outer([0, 1, 0, 0, 0, 0], [0, 0, 5.0])  # an H apart: two other states
NH3_H_F = NH4_F + np.outer([0, 1, 0, 0, 0, 0], [0, 0, 5.0])
LEVEL = Level(program="fake", version="0", method="gfn2", charge=1, multiplicity=2)


@dataclasses.dataclass(frozen=True)
class _Class:  # the attributes of a chemistry.trials class that explore reads
    key: str
    order: tuple
    conformer: int = 0
    formed: tuple = ((0, 2),)
    broken: tuple = ()
    start: np.ndarray | None = None


def _classes(monkeypatch, by_label):
    """trials.trials scripted per state label of the first structure; the start is that
    structure unless given. Returns the calls (label, number of structures)."""
    calls = []

    def fake(symbols, conformers, charge, multiplicity):
        label = state_label(symbols, conformers[0])
        calls.append((label, len(conformers)))
        return [dataclasses.replace(c, start=np.asarray(conformers[c.conformer]))
                if c.start is None else c for c in by_label.get(label, [])]

    monkeypatch.setattr(trials, "trials", fake, raising=False)
    return calls


def _screen(root, composition, symbols, minima, seeds=()):
    """Species, opt calculations and screen minima: (id, coords, energy) per minimum, each its
    own representative; ``seeds`` (id, coords) are members of the first."""
    out, q = [], 1 if composition.endswith("q1_m2") else 0
    for sid, x, energy in [*minima, *((s, x, None) for s, x in seeds)]:
        geo = write_geometry(root, f"in/{sid}.xyz", symbols, x)
        out.append(Artifact(artifact_id=sid, type=SPECIES, payload=SpeciesRecord(
            species_id=sid, composition_id=composition, charge=q, multiplicity=q + 1,
            geometry=geo, source="conformer", state_label=state_label(symbols, x))))
        if energy is None:
            continue
        out.append(Artifact(artifact_id=f"opt_{sid}", type=ArtifactType.CALCULATION,
                            payload=Evidence(engine="fake", task="opt", level=LEVEL, start=geo,
                                             final=geo, energy_hartree=energy, output=geo.file,
                                             job_key=f"opt_{sid}")))
        out.append(Artifact(artifact_id=f"m_{sid}", type=MINIMUM, payload=MinimumRecord(
            minimum_id=f"m_{sid}", basin_id=f"b_{sid}", species_id=sid, tier="screen",
            composition_id=composition, level_key="k", opt_calc=f"opt_{sid}", freq_calc="f",
            energy_hartree=energy, state_label=state_label(symbols, x),
            members=(sid, *(s for s, _ in seeds)) if sid == minima[0][0] else (sid,))))
    return out


def _stage(fake_runtime, artifacts, script, max_trials=3000):
    fake = FakeDiscovery(script)
    rt = fake_runtime(SystemConfig(system_id="t", species=[SpeciesInput(id="hcn", smiles="C#N",
                                                                multiplicity=1)]),
                      {(Capability.DISCOVERY, "readuct"): fake},
                      methods={"gfn2": MethodSpec(id="gfn2", kind="xtb", gfn=2)})
    config = ExploreConfig(engine="readuct", method="gfn2", max_trials=max_trials)
    inputs = Manifest(run_id="r", stage_id="screen", created_at="t", artifacts=artifacts)
    return (lambda runtime: ExploreStage().run(inputs, config, runtime)), rt, fake


def _result(ends=None, ts=None, reason=None, connected=True):
    return DiscoveryResult(outcome="product" if ends else "negative", reason=reason, ends=ends,
                           ts=ts, dE_act_kcal=30.0 if ts else None,
                           dE_rxn_kcal=10.0 if ends else None, irc_connected_to_source=connected)


def _hcn(monkeypatch, fake_runtime, tmp_run, max_trials=5):
    """HCN⁺• screen minima s0 (lowest, holding a seed with H apart that it lost), s1 (0.6
    kcal/mol up) and s2 (12.6 kcal/mol up: outside the window), and an Ar atom. Generation 1:
    'a' (on s1) reaches HNC, 'b' reaches its TS again, 'c' starts from a structure that relaxes
    into HNC, 'd' fails. Generation 2 from HNC: 'e' finds nothing, 'f' is over the budget."""
    geo = {n: write_geometry(tmp_run, f"fake/{n}.xyz", HCN_SYMBOLS, x)
           for n, x in (("hnc", HNC), ("ts", (HCN + HNC) / 2), ("start", HCN + 0.3))}
    artifacts = _screen(tmp_run, "CHN_q1_m2", HCN_SYMBOLS,
                        [("s0", HCN, 0.0), ("s1", HCN + 0.01, 0.001), ("s2", HCN + 0.02, 0.02)],
                        seeds=[("seed", H_CN)])
    artifacts += _screen(tmp_run, "Ar_q0_m1", ["Ar"], [("ar", np.zeros((1, 3)), -1.0)])
    calls = _classes(monkeypatch, {
        state_label(HCN_SYMBOLS, HCN): [
            _Class("d", (3, 0.1)), _Class("b", (2, 0.9)), _Class("a", (1, 0.5), conformer=1),
            _Class("c", (2, 1.2), start=HCN + 0.3)],
        state_label(HCN_SYMBOLS, HNC): [_Class("f", (2, 0.4)), _Class("e", (1, 0.3))]})

    def script(source, trial, method, settings):
        key = next(k for k in "abcdef" if trial.trial_id.endswith(_ids()[k]))
        return {"a": _result((geo["start"], geo["hnc"]), geo["ts"]),
                "b": _result((geo["start"], geo["hnc"]), geo["ts"]),
                "c": _result((geo["start"], geo["hnc"]), None, connected=False),
                "d": Failure(kind=FailureKind.NONZERO_EXIT, reason="scripted"),
                "e": _result(reason="no_nt2_maximum")}[key]

    run, rt, fake = _stage(fake_runtime, artifacts, script, max_trials=max_trials)
    return run, rt, fake, calls


def _ids():
    return {k: sha256_text(f"CHN_q1_m2|{k}") for k in "abcdef"}


def test_edges_over_generations_under_one_budget(monkeypatch, fake_runtime, tmp_run):
    run, rt, fake, calls = _hcn(monkeypatch, fake_runtime, tmp_run)
    out = run(rt)
    found = {a.artifact_id: a for a in out if a.type == DISCOVERY}
    new = {a.payload.species_id: a.payload for a in out if a.type == SPECIES}
    key = {f"disc_trial_{v}": k for k, v in _ids().items()}
    got = [(key.get(i, i), a.payload.outcome, a.payload.generation) for i, a in found.items()]
    assert got == [("relax_seed", "product", 1), ("a", "product", 1), ("b", "negative", 1),
                   ("c", "product", 1), ("d", "failed", 1), ("e", "negative", 2),
                   ("f", "not_attempted", 2)]  # the budget order; f is beyond max_trials = 5
    a, b, c = (found[f"disc_trial_{_ids()[k]}"].payload for k in "abc")
    assert (a.source_species, a.trial.kind, a.ts is not None) == ("s1", "f1b0", True)
    hnc = new[a.product_species]
    assert (hnc.source, hnc.state_label) == ("discovery", state_label(HCN_SYMBOLS, HNC))
    assert b.reason == f"same_edge:{a.discovery_id}"  # one TS, one edge
    assert c.ts is None and c.product_species == a.product_species  # joins the HNC species
    assert new[c.source_species].state_label == state_label(HCN_SYMBOLS, HCN)  # its own start
    assert found[f"disc_trial_{_ids()['d']}"].failure.kind == FailureKind.NONZERO_EXIT
    relax = found["relax_seed"].payload  # the lost seed, unrelaxed, as its own species, to s0
    assert (relax.mechanism, relax.ts, relax.product_species) == ("relaxation", None, "s0")
    assert new[relax.source_species].state_label == state_label(HCN_SYMBOLS, H_CN)
    # classes once per state: HCN on s0 and s1 (s2 is outside the window), HNC on its IRC end;
    # Ar has none; the start state of the seed has no minimum; c's start is no minimum
    assert sorted(calls) == sorted([(state_label(HCN_SYMBOLS, HCN), 2),
                                    (state_label(HCN_SYMBOLS, HNC), 1), (state_label(["Ar"],
                                                                         np.zeros((1, 3))), 1)])
    assert len(fake.calls) == 5
    started = {t.trial_id: s.xyz.coords for s, t, *_ in fake.calls}
    assert np.allclose(started[f"trial_{_ids()['a']}"], HCN + 0.01)  # its conformer, s1
    assert np.allclose(started[f"trial_{_ids()['c']}"], HCN + 0.3)  # its own start


def test_a_trial_that_raises_fails_alone(monkeypatch, fake_runtime, tmp_run):
    """X7-2: an exception in one attempt is contained by StageRuntime.contain: its discovery
    is failed (error:<type>), the stage counts one error (exit code 1) and the others stand."""
    monkeypatch.delenv("HFAUTO_STRICT")
    run, rt, fake, _ = _hcn(monkeypatch, fake_runtime, tmp_run)
    scripted, broken = fake._script, f"trial_{_ids()['d']}"

    def script(source, trial, *args):
        if trial.trial_id == broken:
            raise RuntimeError("parser broke")
        return scripted(source, trial, *args)

    fake._script = script
    found = {a.artifact_id: a for a in run(rt) if a.type == DISCOVERY}
    failed = found[f"disc_{broken}"]
    assert failed.payload.outcome == "failed" and failed.failure.kind == FailureKind.ERROR
    assert failed.payload.reason == "error:error:RuntimeError: parser broke"
    assert rt.errors == [broken] and found[f"disc_trial_{_ids()['a']}"].payload.outcome == (
        "product")


def test_a_state_left_once_the_budget_is_spent_is_not_attempted(monkeypatch, fake_runtime,
                                                                 tmp_run):
    """With 4 attempts generation 1 spends the budget: HNC is reached but not run from, and its
    classes are not enumerated."""
    run, rt, fake, calls = _hcn(monkeypatch, fake_runtime, tmp_run, max_trials=4)
    out = [a.payload for a in run(rt) if a.type == DISCOVERY]
    left = out[-1]
    assert (left.outcome, left.trial, left.generation) == ("not_attempted", None, 2)
    assert left.source_species == out[1].product_species  # the HNC species of 'a'
    assert len(fake.calls) == 4 and state_label(HCN_SYMBOLS, HNC) not in {c for c, _ in calls}


def test_parallel_units_are_recorded_as_serial_ones(monkeypatch, fake_runtime, tmp_run):
    """Artifacts in input order and the same engine calls (the JobStore keys) on 1 or 4 cores."""
    run, rt, fake, _ = _hcn(monkeypatch, fake_runtime, tmp_run)
    parallel, calls = run(rt), list(fake.calls)
    fake.calls.clear()
    serial = run(dataclasses.replace(rt, site=rt.site.model_copy(update={"cores": 1})))
    assert rt.site.cores == 4 and [a.model_dump() for a in parallel] == [
        a.model_dump() for a in serial]
    assert sorted(repr((m.fingerprint(), t)) for m, t, *_ in calls) == sorted(
        repr((m.fingerprint(), t)) for m, t, *_ in fake.calls)


def test_unconnected_edges_and_edges_reached_through_their_product(monkeypatch, fake_runtime,
                                                                    tmp_run):
    """R5a (S19, S6): 'x' ends in NH4+F- and in the source relabelled (two H swapped), neither
    with the source's bonds atom by atom: both ends are species of their own, and the edge runs
    from the relabelled source, the end nearer the start states. 'y' joins two states no edge
    reaches: unconnected, and neither is expanded."""
    symbols = list(NH3_HF_SYMBOLS)
    geo = {n: write_geometry(tmp_run, f"fake/{n}.xyz", symbols, x) for n, x in (
        ("nh4f", NH4_F), ("image", NH3_HF_EXCHANGED + [1.0, -2.0, 0.5]), ("a", NH2_H_FH),
        ("b", NH3_H_F), ("ts", (NH3_HF + NH4_F) / 2), ("ts2", (NH2_H_FH + NH3_H_F) / 2))}
    calls = _classes(monkeypatch, {state_label(symbols, NH3_HF): [
        _Class("x", (2, 1.0), formed=((0, 4),), broken=((4, 5),)), _Class("y", (2, 2.0))]})

    def script(source, trial, method, settings):
        x = trial.associations == ((0, 4),)
        ends = (geo["nh4f"], geo["image"]) if x else (geo["a"], geo["b"])
        return _result(ends, geo["ts"] if x else geo["ts2"], connected=False)

    run, rt, _ = _stage(fake_runtime, _screen(tmp_run, "FH4N_q0_m1", symbols,
                                              [("s0", NH3_HF, 0.0)]), script)
    out = run(rt)
    x, y = (a.payload for a in out if a.type == DISCOVERY)
    new = {a.payload.species_id: a.payload for a in out if a.type == SPECIES}
    assert (x.outcome, x.generation) == ("product", 1)
    assert new[x.source_species].state_label == state_label(symbols, NH3_HF)  # relabelled
    assert new[x.product_species].state_label == state_label(symbols, NH4_F)
    assert (x.dE_act_kcal, x.dE_rxn_kcal) == (20.0, -10.0)  # read from its source end
    assert (y.outcome, y.generation) == ("unconnected", None)
    assert {new[s].state_label for s in (y.source_species, y.product_species)} == {
        state_label(symbols, NH2_H_FH), state_label(symbols, NH3_H_F)}
    assert [label for label, _ in calls] == [state_label(symbols, NH3_HF),
                                             state_label(symbols, NH4_F)]  # not NH3 + HF


def _nh3_hf(nh):
    """C3v NH3·HF with N···H at r_thr + nh (Å) and H–F 0.95 Å."""
    r = covalent_radius("N") + covalent_radius("H") + BOND_TOLERANCE_A + nh
    return np.array([[0.0, 0.0, 0.0], [0.94, 0.0, -0.38], [-0.47, 0.814, -0.38],
                     [-0.47, -0.814, -0.38], [0.0, 0.0, r], [0.0, 0.0, r + 0.95]])


def test_one_edge_without_a_ts_per_lost_seed_state(tmp_run):
    """R6 (X3): NH3 + HF seeds collapsed into N–H-bonded basins; no screen minimum keeps their
    state. It is lost through the seed of lowest low-level energy that a resolved bond change
    separates from its basin: 'near' (lowest) crosses the threshold by ±0.006 Å as S19's c01
    did, 'far' (+2.00 / −0.35 Å, S6's H–O) leaves from its own basin m1, 'late' has no energy;
    'deep' holds the kept state. The edge runs from the seed, unrelaxed as a species of its
    own, to the basin's representative. Without m1 nothing is lost."""
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
    net = _Network(xyz_loader(tmp_run), set())
    _relaxations(net, species, minima, final)
    [(record, _)] = net.records
    assert (record.discovery_id, record.product_species) == ("relax_far", "far")
    assert (record.mechanism, record.outcome, record.ts) == ("relaxation", "product", None)
    [seed] = net.species
    assert (seed.species_id, record.source_species) == ("spc_relax_far_0", "spc_relax_far_0")
    assert (seed.geometry, seed.state_label) == (species["far"].geometry,
                                                 species["far"].state_label)
    empty = _Network(xyz_loader(tmp_run), set())
    _relaxations(empty, species, minima[:1], final)
    assert empty.records == []


def test_an_end_joins_a_known_minimum_only_as_labelled(tmp_run):
    """X2: the NH3·HF double H exchange is the source's basin by permutation-invariant RMSD, but
    with other atom-indexed bonds: it is a species of its own (an edge's ends keep its
    labelling), which a later end joins; the source's own labelling joins the source."""
    symbols = list(NH3_HF_SYMBOLS)
    net = _Network(xyz_loader(tmp_run), set())
    source = SpeciesRecord(species_id="source", composition_id="FH4N_q0_m1", charge=0,
                           multiplicity=1, source="conformer", state_label="x",
                           geometry=write_geometry(tmp_run, "in/source.xyz", symbols, NH3_HF))
    net.minima[("FH4N_q0_m1", "x")].append((source, XYZ(symbols, NH3_HF)))

    def end(name, coords):
        return net.end(source, write_geometry(tmp_run, f"in/{name}.xyz", symbols, coords), name)

    assert end("p1", NH3_HF_EXCHANGED) == "p1"
    assert end("p2", NH3_HF_EXCHANGED + 0.01) == "p1"
    assert end("p3", NH3_HF + 0.01) == "source"
    assert [s.species_id for s in net.species] == ["p1"]


def test_a_state_reached_in_two_atom_orders_is_enumerated_per_order(tmp_run):
    """W4 S6: OH + CH4 and CH3 + H2O are two declared compositions of one composition id, so a
    state holds structures in two atom orders. Each order is enumerated with its own symbols:
    the classes are those of either order alone, once each, and every start keeps the element
    of each atom (before the fix a structure was read with the other order's symbols)."""
    symbols = list(NH3_HF_SYMBOLS)
    order = [5, 4, 3, 2, 1, 0]  # the same structure, atoms reversed
    other = [symbols[i] for i in order]
    net = _Network(xyz_loader(tmp_run), set())
    state = ("FH4N_q0_m1", state_label(symbols, NH3_HF))
    for sid, syms, x in (("a", symbols, NH3_HF), ("b", other, NH3_HF[order])):
        record = SpeciesRecord(species_id=sid, composition_id=state[0], charge=0,
                               multiplicity=1, source="conformer", state_label=state[1],
                               geometry=write_geometry(tmp_run, f"in/{sid}.xyz", syms, x))
        net.minima[state].append((record, XYZ(syms, x)))
    units = net.units(state)
    alone = {t.key for t in trials.trials(symbols, [NH3_HF], 0, 1)}
    assert len(units) == len(alone) == len({u.trial.trial_id for u in units})
    assert {u.trial.trial_id for u in units} == {
        "trial_" + sha256_text(f"{state[0]}|{k}") for k in alone}
    for u in units:
        assert list(u.start.xyz.symbols) == list(u.source.geometry.symbols)


def test_states_beyond_the_budget_in_bottleneck_order_are_not_enumerated(monkeypatch, tmp_run):
    """W4 S19/S10: a generation enumerates its states lowest bottleneck first and stops once the
    budget is covered; the higher state is left (one not_attempted record, no enumeration)."""
    calls = _classes(monkeypatch, {state_label(HCN_SYMBOLS, HCN): [_Class("p", (1, 0.1)),
                                                                   _Class("q", (1, 0.2))],
                                   state_label(HCN_SYMBOLS, HNC): [_Class("r", (1, 0.0))]})
    net = _Network(xyz_loader(tmp_run), {("CHN_q1_m2", state_label(HCN_SYMBOLS, x))
                                         for x in (HCN, HNC)})
    for sid, x, level in (("low", HCN, 0.0), ("high", HNC, 5.0)):
        state = ("CHN_q1_m2", state_label(HCN_SYMBOLS, x))
        record = SpeciesRecord(species_id=sid, composition_id=state[0], charge=1, multiplicity=2,
                               source="conformer", state_label=state[1],
                               geometry=write_geometry(tmp_run, f"in/{sid}.xyz", HCN_SYMBOLS, x))
        net.minima[state].append((record, XYZ(HCN_SYMBOLS, x)))
        net.level[state] = (level, 0.0)
    nothing = DiscoveryResult(outcome="negative", reason="no_nt2_maximum", ends=None, ts=None,
                              dE_act_kcal=None, dE_rxn_kcal=None, irc_connected_to_source=False)
    net.expand(lambda units: [nothing] * len(units), 2)
    assert [label for label, _ in calls] == [state_label(HCN_SYMBOLS, HCN)]
    assert [(d.outcome, d.trial is None, d.source_species) for d, _ in net.records] == [
        ("negative", False, "low"), ("negative", False, "low"), ("not_attempted", True, "high")]

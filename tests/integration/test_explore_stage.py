"""explore stage with FakeDiscovery (§8.2): product, NT2 → AFIR fallback, out of window, failed
attempts, relaxation discovery (legacy relaxation K case), two lowest sources per state and the
screen-optimized start structure."""

import numpy as np
from fakes import FakeDiscovery, write_geometry

from hfauto.backends.protocols import Capability, DiscoveryResult
from hfauto.chemistry.topology import state_label
from hfauto.chemistry.trials import perturb_linear
from hfauto.core.evidence import Evidence, Failure, FailureKind, Level
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec
from hfauto.core.records import ArtifactType, MinimumRecord, SpeciesRecord
from hfauto.core.system import SpeciesInput, SystemConfig
from hfauto.stages.explore import ExploreConfig, ExploreStage

SYMBOLS, SPECIES, MINIMUM = ["H", "C", "N"], ArtifactType.SPECIES, ArtifactType.MINIMUM
HCN = np.array([[0, 0, -1.066], [0, 0, 0], [0, 0, 1.156]])
HNC = np.array([[0, 0, 2.17], [0, 0, 0], [0, 0, 1.17]])
LABEL = state_label(SYMBOLS, HCN)
LEVEL = Level(program="fake", version="0", method="gfn2", charge=0, multiplicity=1)


def _optimized(i):  # the screen-optimized HCN of m{i}: the input displaced by 0.02 (i + 1) Å
    return HCN + 0.02 * (i + 1)


def _inputs(root):
    """Three HCN screen minima of one state (each with its opt calculation); m0 also holds a
    seed that relaxed from HNC. m3 is an Ar atom: never a source."""
    geo, ar = (write_geometry(root, f"in/{n}.xyz", s, x) for n, s, x in (
        ("hcn", SYMBOLS, HCN), ("ar", ["Ar"], np.zeros((1, 3)))))
    out = [Artifact(artifact_id=sid, type=SPECIES, payload=SpeciesRecord(
        species_id=sid, composition_id=comp, charge=0, multiplicity=1, geometry=g,
        source="conformer", state_label=state_label(g.symbols, x)))
        for sid, comp, g, x in (("seed", "CHN_q0_m1", geo, HNC), ("s0", "CHN_q0_m1", geo, HCN),
                                ("s1", "CHN_q0_m1", geo, HCN), ("s2", "CHN_q0_m1", geo, HCN),
                                ("s3", "Ar_q0_m1", ar, np.zeros((1, 3))))]
    for i in range(4):
        final = ar if i == 3 else write_geometry(root, f"opt/m{i}.xyz", SYMBOLS, _optimized(i))
        out.append(Artifact(artifact_id=f"opt{i}", type=ArtifactType.CALCULATION, payload=Evidence(
            engine="fake", task="opt", level=LEVEL, start=final if i == 3 else geo, final=final,
            energy_hartree=i / 10, output=final.file, job_key=f"opt{i}")))
    return Manifest(run_id="r", stage_id="screen", created_at="t", artifacts=out + [
        Artifact(artifact_id=f"m{i}", type=MINIMUM, payload=MinimumRecord(
            minimum_id=f"m{i}", basin_id=f"b{i}", species_id=f"s{i}", tier="screen",
            composition_id="Ar_q0_m1" if i == 3 else "CHN_q0_m1", level_key="k",
            opt_calc=f"opt{i}", freq_calc="f", energy_hartree=i / 10,
            state_label=state_label(["Ar"], np.zeros((1, 3))) if i == 3 else LABEL,
            members=(f"s{i}", "seed")[: 2 - min(i, 1)]))
        for i in range(4)])


def test_explore_records_every_attempt(fake_runtime, tmp_run):
    hnc = write_geometry(tmp_run, "fake/hnc.xyz", SYMBOLS, HNC)

    def script(source, trial, method, settings):  # only the H shift finds a product
        if trial.mechanism == "afir":
            return Failure(kind=FailureKind.NONZERO_EXIT, reason="scripted")
        ok = trial.kind == "h_shift"
        return DiscoveryResult(
            outcome="product" if ok else "negative", reason=None if ok else "monotonic_uphill",
            product=hnc if ok else None, ts=hnc if ok else None, ts_imag_cm1=-1400.0,
            barrier_kj_mol=100.0 if trial.source_minimum == "m0" else 200.0, reaction_kj_mol=8.0,
            irc_connected_to_source=ok, electronic_temperature_K=300.0, job_key="j")

    fake = FakeDiscovery(script)
    rt = fake_runtime(SystemConfig(system_id="t", species=[SpeciesInput(id="hcn", smiles="C#N")]),
                      {(Capability.DISCOVERY, "readuct"): fake},
                      methods={"gfn2": MethodSpec(id="gfn2", kind="xtb", gfn=2)})
    out = ExploreStage().run(_inputs(tmp_run), ExploreConfig(engine="readuct", method="gfn2"), rt)
    found = {(a.payload.source_minimum, a.payload.trial and a.payload.trial.kind,
              a.payload.mechanism): a for a in out if a.type == ArtifactType.DISCOVERY}
    got = {key: (a.payload.outcome, a.payload.reason) for key, a in found.items()}
    assert got[("m0", "h_shift", "nt2")] == ("product", None)
    assert got[("m1", "h_shift", "nt2")] == ("negative", "out_of_window")  # 200 > 150 kJ/mol
    assert ("m0", "h_shift", "afir") not in got and ("m1", "h_shift", "afir") not in got
    assert got[("m0", "heavy_bond", "nt2")] == ("negative", "monotonic_uphill")
    assert found[("m0", "heavy_bond", "afir")].failure.kind == FailureKind.NONZERO_EXIT
    assert got[("m0", None, "relaxation")] == ("negative", f"collapsed_to:{LABEL}")
    assert {c[1].source_minimum for c in fake.calls} == {"m0", "m1"}  # m2: third lowest; m3: Ar
    assert all(c[1].perturbed for c in fake.calls)  # linear HCN starts bent
    for source, trial, *_ in fake.calls:  # from the opt final, not from species.geometry
        i = int(trial.source_minimum[1:])
        assert np.allclose(source.xyz.coords, perturb_linear(_optimized(i)))
    [species] = [a.payload for a in out if a.type == SPECIES]
    assert (species.source, species.state_label) == ("discovery", state_label(SYMBOLS, HNC))
    assert found[("m0", "h_shift", "nt2")].payload.product_species == species.species_id

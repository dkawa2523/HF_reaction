"""explore stage with FakeDiscovery (§8.2): product, NT2 → AFIR fallback, out of window, failed
attempts, relaxation discovery (legacy relaxation K case) and two lowest sources per state."""

import numpy as np
from fakes import FakeDiscovery, write_geometry

from hfauto.backends.protocols import Capability, DiscoveryResult
from hfauto.chemistry.topology import state_label
from hfauto.core.evidence import Failure, FailureKind
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec
from hfauto.core.records import ArtifactType, MinimumRecord, SpeciesRecord
from hfauto.core.system import SpeciesInput, SystemConfig
from hfauto.stages.explore import ExploreConfig, ExploreStage

SYMBOLS, SPECIES, MINIMUM = ["H", "C", "N"], ArtifactType.SPECIES, ArtifactType.MINIMUM
HCN = np.array([[0, 0, -1.066], [0, 0, 0], [0, 0, 1.156]])
HNC = np.array([[0, 0, 2.17], [0, 0, 0], [0, 0, 1.17]])
LABEL = state_label(SYMBOLS, HCN)


def _inputs(root):
    """Three HCN screen minima of one state; m0 also holds a seed that relaxed from HNC."""
    geo = write_geometry(root, "in/hcn.xyz", SYMBOLS, HCN)
    out = [Artifact(artifact_id=sid, type=SPECIES, payload=SpeciesRecord(
        species_id=sid, composition_id="CHN_q0_m1", charge=0, multiplicity=1,
        geometry=geo, source="conformer", state_label=state_label(SYMBOLS, x)))
        for sid, x in (("seed", HNC), ("s0", HCN), ("s1", HCN), ("s2", HCN))]
    return Manifest(run_id="r", stage_id="screen", created_at="t", artifacts=out + [
        Artifact(artifact_id=f"m{i}", type=MINIMUM, payload=MinimumRecord(
            minimum_id=f"m{i}", basin_id=f"b{i}", composition_id="CHN_q0_m1", species_id=f"s{i}",
            tier="screen", level_key="k", opt_calc="o", freq_calc="f", energy_hartree=i / 10,
            state_label=LABEL, members=(f"s{i}", "seed")[: 2 - min(i, 1)]))
        for i in range(3)])


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
    rt = fake_runtime(SystemConfig(system_id="t", species=[SpeciesInput(id="hcn")]),
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
    assert {c[1].source_minimum for c in fake.calls} == {"m0", "m1"}  # m2 is the third lowest
    assert all(c[1].perturbed for c in fake.calls)  # linear HCN starts bent
    [species] = [a.payload for a in out if a.type == SPECIES]
    assert (species.source, species.state_label) == ("discovery", state_label(SYMBOLS, HNC))
    assert found[("m0", "h_shift", "nt2")].payload.product_species == species.species_id

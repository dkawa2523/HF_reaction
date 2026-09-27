"""configs/pipelines/*.yaml end to end on fake engines (design §4.3, §7.4, §10.1; CH-05): the
double_well declared N–H–O reaction plus one composition, (HF)2, on 0.92 Å pair springs."""

from pathlib import Path

import fakes
import numpy as np
import pytest
import yaml
from scipy.spatial.distance import pdist
from typer.testing import CliRunner

from hfauto.backends.protocols import Capability as Cap
from hfauto.backends.protocols import ConformerEnsemble, DiscoveryResult, ThermoResult
from hfauto.chemistry.gates import zpe_hartree
from hfauto.chemistry.thermo import settings_sha, thermo_frequencies
from hfauto.cli.main import app
from hfauto.core import records as R
from hfauto.core.evidence import Evidence
from hfauto.core.manifest import load_manifest
from hfauto.pipeline.config import load
from hfauto.pipeline.layout import RunLayout
from hfauto.pipeline.runner import run_pipeline
from hfauto.stages import catalog

pytestmark = pytest.mark.integration
REPO = Path(__file__).resolve().parents[2]
T, PIPELINES = R.ArtifactType, REPO / "configs" / "pipelines"


@pytest.fixture
def pipeline(tmp_path, tmp_run, override_engine, monkeypatch):
    """``pipeline(name)`` runs configs/pipelines/<name>.yaml in tmp_run with every engine faked."""
    monkeypatch.chdir(REPO)  # the methods come from configs/methods, as for the CLI
    well = fakes.double_well()
    pes = fakes.PES(well.symbols, lambda x: well.energy(x) if np.size(x) == 9 else float(
        np.sum((pdist(np.reshape(x, (-1, 3))) - 0.92) ** 2)), well.points)
    geo = {p: fakes.write_geometry(tmp_run, f"in/{p}.xyz", well.symbols, x)
           for p, x in well.points.items()}  # also the endpoints' input xyz
    fakes.write_geometry(tmp_run, "in/hf.xyz", ["H", "F"], [[0, 0, 0], [0, 0, 0.92]])

    def search(mol, method, settings):  # the (HF)2 placement seed is the only conformer
        seed = fakes.write_geometry(tmp_run, "fake/crest.xyz", mol.xyz.symbols, mol.xyz.coords)
        return ConformerEnsemble(members=((seed, pes.energy(mol.xyz.coords)),), version="0",
                                 topology_removed=0, topology_stops=(), job_key="crest")

    def explore(source, trial, method, settings):  # N–H–O reaches the other well
        ok, x = source.xyz.coords.size == 9, source.xyz.coords
        return DiscoveryResult(
            outcome="product" if ok else "negative", reason=None if ok else "monotonic_uphill",
            product=geo["product" if x[1, 0] < 0 else "reactant"] if ok else None, job_key="rd",
            ts=geo["ts"] if ok else None, ts_imag_cm1=-1e3, barrier_kj_mol=30.0,
            reaction_kj_mol=10.0, irc_connected_to_source=ok, electronic_temperature_K=300.0)

    def goodvibes(freq, settings, temperatures, saddle):  # consistent: G = E + scaled ZPE
        e, nu = freq.energy_hartree, thermo_frequencies(freq.frequencies_cm1, saddle=saddle)
        z = zpe_hartree(nu, scale=settings[0].vib_scale)
        return [ThermoResult(settings_sha=settings_sha(s), T_K=t, E_hartree=e, H_hartree=e,
                             G_hartree=e + z, zpe_hartree=z, S_rot=1.0, notes=(), job_key="gv",
                             n_real=sum(f > 0 for f in nu))
                for s in settings for t in temperatures]

    qm, path, saddle = (c(tmp_run, pes) for c in (fakes.FakeQM, fakes.FakePath, fakes.FakeSaddle))
    fake = {(Cap.QM, "nwchem"): qm, (Cap.QM, "xtb"): qm, (Cap.SADDLE, "nwchem_saddle"): saddle,
            (Cap.PATH, "nwchem_string"): path, (Cap.PATH, "pysis_neb"): path,
            (Cap.CONFORMERS, "crest"): fakes.FakeConformers(search),
            (Cap.DISCOVERY, "readuct"): fakes.FakeDiscovery(explore),
            (Cap.THERMO, "goodvibes"): fakes.FakeThermo(goodvibes)}
    for (capability, name), engine in fake.items():
        override_engine(capability, name, engine)
    ends = [{"id": n, "xyz": f"run/in/{n}.xyz", "role": "endpoint"}
            for n in ("reactant", "product")]
    cfg = {"system": {"system_id": "dw", "species": [*ends, {"id": "hf", "xyz": "run/in/hf.xyz"}],
                      "compositions": [{"id": "hf2", "components": {"hf": 2}}],
                      "reactions": [{"id": "rx", "reactant": "reactant", "product": "product"}]},
           "site": {"site": "fake", "scratch_root": str(tmp_path / "scratch"), "cores": 1,
                    "engines": {n: {"version": "0"} for _, n in fake}}}
    for name, data in cfg.items():
        (tmp_path / f"{name}.yaml").write_text(yaml.safe_dump(data), encoding="utf-8")

    def run(name):
        resolved = load(PIPELINES / f"{name}.yaml", *(tmp_path / f"{f}.yaml" for f in cfg))
        layout = run_pipeline(resolved, tmp_run)
        for entry in resolved.pipeline.stages:  # types on disk stay within spec.produces
            produced = {a.type for a in load_manifest(layout.manifest_path(entry.id)).artifacts}
            assert produced <= set(catalog.get(entry.stage).spec.produces), entry.id
        return layout

    return run


def test_known_endpoints_then_method_panel_appended_to_the_run(pipeline):
    view = pipeline("known_endpoints").view()  # manifests read back from disk, typed (CH-05)
    rx = next(r for r in view.records(T.REACTION, R.ReactionRecord) if r.reaction_id == "rx")
    assert rx.source == "declared" and rx.outcome is R.CaseOutcome.ELEMENTARY_STEP
    (report,) = view.records(T.REPORT, R.ReportRecord)
    assert [row.reaction_id for row in report.rows if row.rankable] == ["rx"]
    layout = pipeline("method_panel")
    view = layout.view("panel_report")  # includes the known_endpoints stages
    panel = {e.level.basis for e in view.records(T.CALCULATION, Evidence) if e.task == "sp"}
    assert panel >= {"def2-tzvpd"} and "def2-svp" not in panel
    report = load_manifest(layout.manifest_path("panel_report")).records(T.REPORT, R.ReportRecord)
    assert [row.reaction_id for row in report[0].rows] == ["rx"]  # the paths stage's reaction
    [row] = report[0].rows  # ranked on panel_thermo's layer (no conditions copied)
    assert row.rankable and (row.T_K, row.energy_level) == (298.15, "wb97x-d3/def2-tzvpd")
    with pytest.raises(ValueError, match="already belongs"):
        pipeline("discover")  # structures, dft, ... are known_endpoints stage ids


def test_discover_flows_from_discovery_to_report(pipeline):
    view = pipeline("discover").view()
    found = {d.product_species for d in view.records(T.DISCOVERY, R.DiscoveryRecord)} - {None}
    dft = [m for m in view.records(T.MINIMUM, R.MinimumRecord) if m.tier == "dft"]
    assert found and found <= {s for m in dft for s in m.members}  # refined at the DFT level
    assert "F2H2_q0_m1" in {m.composition_id for m in dft}  # the (HF)2 composition
    steps = {r.reaction_id for r in view.records(T.REACTION, R.ReactionRecord)
             if r.outcome is R.CaseOutcome.ELEMENTARY_STEP}
    (report,) = view.records(T.REPORT, R.ReportRecord)  # rankable: thermo and report ran
    assert steps & {row.reaction_id for row in report.rows if row.rankable}


def test_a_run_whose_inputs_all_fail_still_writes_the_report(pipeline, tmp_path, tmp_run,
                                                            monkeypatch):
    monkeypatch.delenv("HFAUTO_STRICT")
    system = yaml.safe_load((tmp_path / "system.yaml").read_text(encoding="utf-8"))
    for species in system["species"]:  # an even electron count: no doublet
        species["multiplicity"] = 2
    (tmp_path / "doublets.yaml").write_text(yaml.safe_dump(system), encoding="utf-8")
    args = ["run", str(PIPELINES / "known_endpoints.yaml"), "--system",
            str(tmp_path / "doublets.yaml"), "--site", str(tmp_path / "site.yaml"),
            "--run-dir", str(tmp_run)]
    result = CliRunner().invoke(app, args)
    assert result.exit_code == 1, result.output
    states = {s.stage_id: s for s in RunLayout(tmp_run).read_state()}
    assert (states["structures"].n_ok, states["structures"].n_failed) == (0, 3)
    assert all((states[s].status, states[s].n_ok, states[s].n_failed) == ("done", 0, 0)
               for s in ("dft", "paths", "thermo"))
    assert (tmp_run / "report" / "report.html").is_file()

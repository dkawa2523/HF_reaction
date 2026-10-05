from __future__ import annotations

import pytest
from pydantic import ValidationError

from hfauto.core import records as r
from hfauto.core.evidence import Evidence, Failure, FileRef, Geometry, Level
from hfauto.core.files import write_atomic
from hfauto.core.manifest import Artifact, Manifest, load_manifest, save_manifest
from hfauto.core.records import ArtifactType as AT

REF = FileRef(path="jobs/ab/out.txt", sha256="f" * 64)
GEO = Geometry(file=REF, fingerprint="fp", symbols=("H", "C", "N"))
LEVEL = Level(program="nwchem", version="7.2.3", method="pbe0", basis="def2-svpd",
              dispersion="d3bj", charge=0, multiplicity=1, grid="fine", scf_tol=1e-7)
TERM = (r.StoichTerm(composition_id="CHN|0|1", coefficient=1),)
ELEMENTARY = r.CaseOutcome.ELEMENTARY_STEP
PAYLOADS = {
    AT.CALCULATION: Evidence(
        engine="nwchem", task="freq", level=LEVEL, start=GEO, final=GEO, energy_hartree=-93.4,
        frequencies_cm1=(-1131.6, 900.0, 2100.0), n_external=6,
        imaginary_modes=((0.1,) * 9,), hessian=REF, s2=None, output=REF, job_key="k"),
    AT.SPECIES: r.SpeciesRecord(
        species_id="s1", composition_id="CHN|0|1", charge=0, multiplicity=1,
        geometry=GEO, source="input", state_label="CHN:ab12cd34"),
    AT.MINIMUM: r.MinimumRecord(
        minimum_id="m1", basin_id="b1", composition_id="CHN|0|1", species_id="s1", tier="dft",
        level_key=LEVEL.full_key(), opt_calc="c1", freq_calc="c2", energy_hartree=-93.4,
        state_label="CHN:ab12cd34", members=("s1",), notes=("soft",)),
    AT.DISCOVERY: r.DiscoveryRecord(
        discovery_id="d1", mechanism="nt2", outcome="product",
        trial=r.ReactionTrial(trial_id="t1", kind="f1b0", associations=((0, 2),)),
        source_species="s1", product_species="s2", ts=GEO, dE_act_kcal=45.4, generation=2),
    AT.REACTION: r.ReactionRecord(
        reaction_id="r1", reactants=TERM, products=TERM, minima=("m1", "m2"),
        endpoints=("s1", "s2"), source="declared", ts_calc="c3",
        coordinate=(r.CoordinateTerm(kind="angle", atoms=(0, 1, 2)),),
        barrier=r.BarrierVerdict(verdict="single", source="string"),
        saddle=r.SaddleClaim(saddle_calc="c3", freq_calc="c4", imag_cm1=-1131.6,
                             energy_hartree=-93.3),
        connection=r.ConnectionClaim(side_calcs=("c5", "c6"), minima=("m1", "m2")),
        outcome=ELEMENTARY, log="cases/r1/log.jsonl"),
    AT.SPECIES_THERMO: r.SpeciesThermo(
        subject="m1", freq_calc="c2", T_K=298.15, G_hartree=None, H_hartree=None,
        zpe_hartree=None, settings_sha="abc", notes=("thermo_unavailable",)),
    AT.REACTION_THERMO: r.ReactionThermo(
        reaction_id="r1", T_K=298.15, standard_state="1M", dE_act_kcal=46.7, dE_rxn_kcal=14.0,
        dG_act_kcal=42.3, dG_rxn_kcal=12.7, band_kcal=(42.0, 42.6),
        dG_eff_kcal=42.3, notes=("submerged_barrier",)),
    AT.REPORT: r.ReportRecord(
        rows=(r.RankRow(reaction_id="r1", outcome=ELEMENTARY, tier="connected", rankable=True,
                        rank=1, dG_act_kcal=42.3, band_kcal=(42.0, 42.6), blockers=(),
                        dG_eff_kcal=42.3, torsional=True),),
        tables={"ranking.csv": REF}),
}


def art(artifact_id: str, kind: AT = AT.SPECIES, **kw) -> Artifact:
    return Artifact(artifact_id=artifact_id, type=kind, payload=PAYLOADS[kind], **kw)


def manifest(*artifacts: Artifact, stage: str = "s") -> Manifest:
    return Manifest(run_id="run", stage_id=stage, created_at="2026-09-25T00:00:00+00:00",
                    artifacts=list(artifacts))


def test_artifact_validator():
    failure = Failure(kind="timeout", reason="rc 124")
    with pytest.raises(ValidationError):
        Artifact(artifact_id="a", type=AT.SPECIES)
    with pytest.raises(ValidationError):
        Artifact(artifact_id="a", type=AT.MINIMUM, payload=PAYLOADS[AT.SPECIES])
    with pytest.raises(ValidationError):
        Artifact(artifact_id="a", type=AT.MINIMUM, status="failed")
    failed = Artifact(artifact_id="a", type=AT.MINIMUM, status="failed", failure=failure)
    assert failed.payload is None


def test_every_payload_round_trips(tmp_path):
    original = manifest(*(art(f"a_{kind}", kind) for kind in PAYLOADS),
                        Artifact(artifact_id="f", type=AT.CALCULATION, status="failed",
                                 failure=Failure(kind="nonzero_exit", reason="rc -15")))
    path = save_manifest(original, tmp_path / "stage" / "manifest.json")
    assert load_manifest(path) == original
    assert [p.name for p in path.parent.iterdir()] == ["manifest.json"]  # no temp file left
    path.write_text(path.read_text(encoding="utf-8").replace(".v2", ".v1"), encoding="utf-8")
    with pytest.raises(ValidationError):  # no compatibility with other schema versions
        load_manifest(path)


def test_an_older_discovery_still_parses():
    """The baselines (R2, B1) hold discoveries written before the edges."""
    trial = {"trial_id": "t", "source_minimum": "m1", "kind": "transfer",
             "associations": [[0, 2]], "dissociations": [[0, 1]]}
    old = {"kind": "discovery", "discovery_id": "disc_t_nt2", "source_minimum": "m1",
           "mechanism": "nt2", "trial": trial, "outcome": "product", "reason": None,
           "source_species": "s1", "product_species": "s2", "ts": None,
           "ts_imag_cm1": -300.0, "dE_act_kcal": 1.0, "dE_rxn_kcal": 0.5,
           "electronic_temperature_K": 300.0, "ts_calc": None}
    assert r.DiscoveryRecord.model_validate(old).generation is None


def test_an_atomic_write_leaves_the_old_file_or_the_new_one(tmp_path):
    target = write_atomic(tmp_path / "new" / "state.json", "one")
    with pytest.raises(TypeError):
        write_atomic(target, None)  # type: ignore[arg-type]
    assert target.read_text(encoding="utf-8") == "one"
    assert [p.name for p in target.parent.iterdir()] == ["state.json"]


def test_union_keeps_first_order_and_last_wins():
    newer = art("x", parents=("p",))
    union = Manifest.union([manifest(art("x"), art("y")), manifest(art("z"), newer)],
                           run_id="run", stage_id="view")
    assert [a.artifact_id for a in union.artifacts] == ["x", "y", "z"]
    assert union.get("x") == newer and union.stage_id == "view"


def test_queries():
    failed = Artifact(artifact_id="f", type=AT.SPECIES, status="failed",
                      failure=Failure(kind="gate_rejected", reason="no_collision_free_seed"))
    newer = art("s", parents=("p",))
    m = manifest(art("s"), art("c", AT.CALCULATION), failed, newer)
    assert m.get("s") == newer
    assert [a.artifact_id for a in m.of(AT.SPECIES)] == ["s", "s"]
    assert len(m.of(AT.SPECIES, ok_only=False)) == 3
    assert m.records(AT.SPECIES, r.SpeciesRecord)[0].species_id == "s1"
    with pytest.raises(TypeError):
        m.records(AT.SPECIES, r.MinimumRecord)
    assert m.evidence("c") == PAYLOADS[AT.CALCULATION]
    with pytest.raises(ValueError):
        m.evidence("s")
    with pytest.raises(KeyError):
        m.get("missing")

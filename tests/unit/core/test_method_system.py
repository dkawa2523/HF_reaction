from __future__ import annotations

from pathlib import Path

import pytest
from pydantic import ValidationError

from hfauto.core.evidence import Level
from hfauto.core.method import Deadline, EngineSite, MethodSpec, level_mismatches
from hfauto.core.system import load_system

DFT = MethodSpec(id="pbe0", kind="dft", functional="PBE0", basis="def2-SVPD", dispersion="d3bj",
                 grid="fine", scf_energy_tol=1e-7)
OBSERVED = Level(program="nwchem", version="7.2.3", method="pbe0", basis="def2-svpd",
                 dispersion="d3bj", charge=0, multiplicity=1, grid="fine", scf_tol=1e-7)


def observed(**update) -> Level:
    return OBSERVED.model_copy(update=update)


def test_dft_level_mismatches():
    assert level_mismatches(DFT, OBSERVED, version_pin="7.2.3") == []
    cart = level_mismatches(DFT, observed(basis="def2-svpd/cart"), version_pin="7.2.3")
    assert len(cart) == 1 and cart[0].startswith("basis")
    pinned = level_mismatches(DFT, OBSERVED, version_pin="7.2.2")
    assert len(pinned) == 1 and pinned[0].startswith("version")
    xfine = DFT.model_copy(update={"grid": "xfine"})
    assert level_mismatches(xfine, observed(grid="medium"), version_pin="7.2.3")
    no_grid = DFT.model_copy(update={"grid": None})
    assert level_mismatches(no_grid, observed(grid="medium"), version_pin="7.2.3") == []
    assert level_mismatches(DFT, observed(scf_tol=1e-7 * (1 + 1e-12)), version_pin="7.2.3") == []
    assert level_mismatches(DFT, observed(scf_tol=1e-6), version_pin="7.2.3")
    assert level_mismatches(DFT, observed(dispersion=None), version_pin="7.2.3")
    assert level_mismatches(DFT, observed(solvation="cosmo:78.4"), version_pin="7.2.3")


def test_xtb_and_wft_level_mismatches():
    xtb = MethodSpec(id="gfn2", kind="xtb", gfn=2)
    level = Level(program="xtb", version="6.7.1", method="gfn2", charge=0, multiplicity=1)
    assert level_mismatches(xtb, level, version_pin="6.7.1") == []
    assert level_mismatches(xtb, level.model_copy(update={"electronic_temperature_K": 300.0}),
                            version_pin="6.7.1") == []
    hot = level.model_copy(update={"electronic_temperature_K": 1000.0})
    assert level_mismatches(xtb, hot, version_pin="6.7.1")
    mp2 = MethodSpec(id="mp2", kind="wft", wft_method="mp2", basis="def2-TZVP")
    wft = observed(method="mp2", basis="def2-tzvp", dispersion=None)
    assert level_mismatches(mp2, wft, version_pin="7.2.3") == []
    assert level_mismatches(mp2, wft.model_copy(update={"method": "ccsd(t)"}), version_pin="7.2.3")


def test_method_signature_site_pin_and_deadline():
    assert "id" not in DFT.signature()
    assert DFT.model_copy(update={"id": "other"}).signature() == DFT.signature()
    with pytest.raises(ValidationError):
        EngineSite()  # type: ignore[call-arg]
    assert EngineSite(version="7.2.3").execution.timeout_s == 14_400
    assert not Deadline.after(100.0).expired()
    assert 0.0 < Deadline.after(100.0).remaining() <= 100.0
    past = Deadline.after(-1.0)
    assert past.expired() and past.remaining() == 0.0


SYSTEM = """
system_id: hcn
species:
  - {id: hcn, xyz: xyz/hcn.xyz, role: endpoint}
  - {id: hnc, xyz: xyz/hnc.xyz, role: endpoint}
  - {id: hf, smiles: F}
compositions:
  - {id: hcn_hf, components: {hcn: 1, hf: 1}}
reactions:
  - {id: iso, reactant: hcn, product: hnc, coordinate: [{kind: angle, atoms: [0, 1, 2]}]}
"""


def write_system(tmp_path: Path, text: str) -> Path:
    (tmp_path / "systems").mkdir()
    (tmp_path / "systems" / "hcn.yaml").write_text(text, encoding="utf-8")
    return tmp_path / "systems" / "hcn.yaml"


def test_load_system_resolves_relative_xyz(tmp_path):
    config = load_system(write_system(tmp_path, SYSTEM))
    hcn, _, hf = config.species
    assert hcn.xyz == (tmp_path / "systems" / "xyz" / "hcn.xyz").resolve()
    assert hf.xyz is None and hf.multiplicity is None and config.compositions[0].multiplicity is None
    assert config.reactions[0].coordinate[0].atoms == (0, 1, 2)


@pytest.mark.parametrize(("old", "new", "match"), [
    ("{id: hf, smiles: F}", "{id: hcn, smiles: F}", "duplicate species"),
    ("{id: hcn_hf,", "{id: hf,", "duplicate species/composition"),  # would collide in conformers
    ("reactions:\n", "reactions:\n  - {id: iso, reactant: hnc, product: hcn}\n", "duplicate reac"),
    ("product: hnc", "product: hf", "not an endpoint"),  # product is a monomer
    ("product: hnc", "product: missing", "not an endpoint"),  # unknown species
    ("{hcn: 1, hf: 1}", "{hcn: 1, hcl: 1}", "unknown species"),
    ("{hcn: 1, hf: 1}", "{hcn: 1, hf: 0}", "greater than 0"),
    ("{id: hf, smiles: F}", "{id: hf}", "exactly one"),
    ("{id: hf, smiles: F}", "{id: hf, smiles: F, xyz: xyz/hf.xyz}", "exactly one"),
    ("{id: hnc, xyz: xyz/hnc.xyz,", "{id: hnc, smiles: '[C-]#[NH+]',", "must be given as xyz"),
    ("atoms: [0, 1, 2]", "atoms: [0, 1]", "distinct atom"),  # an angle needs three atoms
    ("atoms: [0, 1, 2]", "atoms: [0, -1, 2]", "distinct atom"),
    ("atoms: [0, 1, 2]", "atoms: [0, 1, 1]", "distinct atom"),
])
def test_load_system_rejects_bad_input(tmp_path, old, new, match):
    assert old in SYSTEM
    with pytest.raises(ValidationError, match=match):
        load_system(write_system(tmp_path, SYSTEM.replace(old, new)))

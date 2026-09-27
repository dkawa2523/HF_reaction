import re
from html import unescape

import numpy as np
import pytest
from fakes import FakeQM, double_well, xyz_loader

from hfauto.core import records as rec
from hfauto.core.manifest import Artifact, Manifest
from hfauto.core.method import MethodSpec
from hfauto.reporting.html import MODE_AMPLITUDE_A, render

T = rec.ArtifactType


def view(*artifacts):
    return Manifest(run_id="run", stage_id="report", created_at="now", artifacts=list(artifacts))


def test_page_has_an_inline_energy_svg_and_an_animated_ts_mode(tmp_path):
    pes = double_well()
    freq = FakeQM(tmp_path, pes).frequencies(pes.molecule("ts"),
                                             MethodSpec(id="m", kind="dft", functional="pbe0"))
    saddle = rec.SaddleClaim(saddle_calc="ts", freq_calc="ts", imag_cm1=freq.frequencies_cm1[0],
                             energy_hartree=freq.energy_hartree)
    reaction = rec.ReactionRecord(reaction_id="r<1>", reactants=(), products=(), saddle=saddle,
                                  minima=("a", "b"), endpoints=("a", "b"), source="declared",
                                  outcome=rec.CaseOutcome.ELEMENTARY_STEP)
    thermo = rec.ReactionThermo(reaction_id="r<1>", T_K=298.15, standard_state="1atm",
                                dE_act_kcal=6.0, dE_rxn_kcal=2.0, dzpe_act_kcal=-1.0,
                                dG_act_kcal=5.0, dG_rxn_kcal=2.0, dG_assoc_kcal=-3.0,
                                dG_act_vs_separated_kcal=2.0, blockers=("spin_contaminated",),
                                dG_eff_kcal=5.25, notes=("submerged_barrier",))
    page = render(view(Artifact(artifact_id="ts", type=T.CALCULATION, payload=freq),
                       Artifact(artifact_id="r", type=T.REACTION, payload=reaction),
                       Artifact(artifact_id="t", type=T.REACTION_THERMO, payload=thermo)),
                  tmp_path / "out", load_xyz=xyz_loader(tmp_path)).read_text(encoding="utf-8")
    assert "<svg class='energy'" in page and ">separated<" in page and "spin_contaminated" in page
    assert "δG<sub>eff</sub>" in page and ">5.25<" in page and "submerged_barrier" in page
    assert "<script src='https://cdn.jsdelivr.net/npm/3dmol@" in page and ".vibrate(" in page
    assert "r&lt;1&gt;" in page and "r<1>" not in page  # ids are escaped
    [xyz] = [unescape(x) for x in re.findall(r"data-xyz='([^']*)'", page)]
    atoms = np.array([line.split()[1:] for line in xyz.splitlines()[2:]], dtype=float)
    assert atoms.shape == (3, 6)  # x y z dx dy dz per atom: 3Dmol's vibration columns
    assert np.linalg.norm(atoms[:, 3:], axis=1).max() == pytest.approx(MODE_AMPLITUDE_A, abs=1e-5)


def test_a_view_without_reactions_still_gives_a_page(tmp_path):
    path = render(view(), tmp_path / "empty", load_xyz=xyz_loader(tmp_path))
    page = path.read_text(encoding="utf-8")
    assert page.startswith("<!doctype html>") and "No reactions" in page and "3Dmol" not in page

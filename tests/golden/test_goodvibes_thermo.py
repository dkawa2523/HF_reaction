"""G02 / G06: GoodVibes CSVs of the old CLI path against the thermo gate and deltas."""

import csv

import pytest
from fakes import FakeQM, harmonic

from hfauto.chemistry import thermo as th
from hfauto.chemistry.gates import thermo_consistent
from hfauto.core import records as R
from hfauto.core.constants import HARTREE_TO_KCAL_MOL as H2K
from hfauto.core.method import MethodSpec

TERM = (R.StoichTerm(composition_id="c", coefficient=1),)
RX = R.ReactionRecord(
    reaction_id="r", reactants=TERM, products=TERM, minima=("a", "b"), endpoints=("a", "b"),
    source="declared", outcome=R.CaseOutcome.ELEMENTARY_STEP,
    saddle=R.SaddleClaim(saddle_calc="s", freq_calc="ts", imag_cm1=-680.0, energy_hartree=0.0))


def _st(golden, subject: str, name: str | None) -> R.SpeciesThermo:
    row = next(csv.DictReader(golden.text(f"goodvibes/{name}.csv").splitlines())) if name else {}
    G, zpe = (float(row[k]) if row else None for k in ("qh_gibbs_free_energy", "zpe"))
    return R.SpeciesThermo(subject=subject, freq_calc=subject, T_K=298.15, G_hartree=G,
                           H_hartree=None, zpe_hartree=zpe, settings_sha="csv")


def test_g02_doubled_zpe_fails_the_gate_so_24_kcal_is_not_reproduced(golden, tmp_run):
    ts, trans = _st(golden, "ts", "G02/ts"), _st(golden, "a", "G02/hono_trans")
    assert (ts.G_hartree - trans.G_hartree) * H2K == pytest.approx(24.11, abs=0.01)  # old value
    pes, e_ts = harmonic(), -205.323820905761  # G02/ts.csv scf_energy
    freq = FakeQM(tmp_run, pes).frequencies(pes.molecule("minimum"), MethodSpec(id="m", kind="dft"))
    freq = freq.model_copy(update={"energy_hartree": e_ts, "frequencies_cm1": (  # G01, last block
        -680.11, 630.77, 844.57, 1088.20, 1861.46, 3808.08)})
    sent = th.thermo_frequencies(freq.frequencies_cm1, saddle=True)
    gate = thermo_consistent(freq, frequencies_cm1=sent, gv_zpe_hartree=ts.zpe_hartree,
                             gv_energy_hartree=e_ts, gv_n_real=10, scale=0.985)
    assert not gate and "zpe_mismatch" in gate.reasons  # 0.0375115 = 0.985 x (both blocks)
    assert th.reaction_delta(RX, {"a": trans, "ts": _st(golden, "ts", None)})[0] is None


def test_g06_hcn_to_hnc_free_energies(golden):
    table = {k: _st(golden, k, f"G06/{n}") for k, n in (("a", "hcn"), ("b", "hnc"), ("ts", "ts"))}
    dG_act, dG_rxn, dzpe = th.reaction_delta(RX, table)
    assert dG_rxn == pytest.approx(12.742, abs=1e-3) and dG_act == pytest.approx(42.281, abs=1e-3)
    assert dzpe == pytest.approx((table["ts"].zpe_hartree - table["a"].zpe_hartree) * H2K)
    assert th.reaction_delta(RX.model_copy(update={"degenerate": True}), table)[1] == 0.0
    barrierless = RX.model_copy(update={"outcome": R.CaseOutcome.BARRIERLESS})
    assert th.reaction_delta(barrierless, table)[0] is None  # no connected TS, no barrier

"""Real ReaDuct (marker real, WSL): NT2 from bent HCN, TS counted on projected frequencies."""

import numpy as np
import pytest

from hfauto.backends.protocols import Capability, DiscoveryResult, DiscoverySettings
from hfauto.chemistry import trials
from hfauto.chemistry.xyz import XYZ, Molecule
from hfauto.core.method import MethodSpec
from hfauto.core.records import ReactionTrial


@pytest.mark.real
def test_nt2_on_bent_hcn_finds_hnc_through_a_first_order_saddle(real_engine):
    symbols, hcn = ["H", "C", "N"], np.array([[0, 0, -1.066], [0, 0, 0], [0, 0, 1.156]])
    [shift] = [t for t in trials.trials(symbols, [hcn], 0, 1)
               if (t.formed, t.broken) == (((0, 2),), ((0, 1),))]
    trial = ReactionTrial(trial_id="t", kind="f1b1", associations=shift.formed,
                          dissociations=shift.broken)
    result = real_engine(Capability.DISCOVERY, "readuct").explore(
        Molecule(XYZ(symbols, shift.start), 0, 1), trial, MethodSpec(id="gfn2", kind="xtb", gfn=2),
        DiscoverySettings())
    assert isinstance(result, DiscoveryResult), result
    assert (result.outcome, result.irc_connected_to_source) == ("product", True)
    assert result.ts is not None and result.ends is not None and result.dE_act_kcal > 0

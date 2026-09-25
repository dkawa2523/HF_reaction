"""Real ReaDuct (marker real, WSL): NT2 from bent HCN, TS counted on projected frequencies."""

import numpy as np
import pytest

from hfauto.backends.protocols import Capability, DiscoveryResult, DiscoverySettings
from hfauto.chemistry import trials
from hfauto.chemistry.xyz import XYZ, Molecule
from hfauto.core.method import MethodSpec


@pytest.mark.real
def test_nt2_on_bent_hcn_finds_hnc_through_a_first_order_saddle(real_engine):
    symbols, hcn = ["H", "C", "N"], np.array([[0, 0, -1.066], [0, 0, 0], [0, 0, 1.156]])
    start, drives = trials.generate("hcn", symbols, hcn)
    [shift] = [t for t in drives if t.kind == "h_shift"]
    result = real_engine(Capability.DISCOVERY, "readuct").explore(
        Molecule(XYZ(symbols, start), 0, 1), shift, MethodSpec(id="gfn2", kind="xtb", gfn=2),
        DiscoverySettings())
    assert isinstance(result, DiscoveryResult), result
    assert (result.outcome, result.irc_connected_to_source) == ("product", True)
    assert result.ts is not None and result.ts_imag_cm1 < -50.0 and result.barrier_kj_mol > 0

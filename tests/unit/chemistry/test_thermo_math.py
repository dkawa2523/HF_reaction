"""chemistry.thermo: textbook values, the electronic term from NIST levels, the fixed rule for
negative modes, and the points a reaction reads (one chain with a common zero)."""

import importlib.metadata
import math
from pathlib import Path

import numpy as np
import pytest

from hfauto.chemistry import thermo as th
from hfauto.chemistry.electronic_state import atom_levels, check_electronic_state
from hfauto.chemistry.topology import state_label
from hfauto.chemistry.vibrations import external_basis
from hfauto.chemistry.xyz import XYZ
from hfauto.core import records as R
from hfauto.core.constants import HARTREE_TO_KCAL_MOL as H2K
from hfauto.core.constants import R_KCAL_MOL_K
from hfauto.core.evidence import FileRef, Geometry
from hfauto.core.system import CompositionInput, SpeciesInput, load_system

RT = R_KCAL_MOL_K * 298.15
OH = XYZ(["O", "H"], np.array([[0.0, 0, 0], [0, 0, 0.97]]))
CH4 = XYZ(["C", "H", "H", "H", "H"], np.array(
    [[0.0, 0, 0], [0.63, 0.63, 0.63], [-0.63, -0.63, 0.63], [-0.63, 0.63, -0.63],
     [0.63, -0.63, -0.63]]))
CH3 = XYZ(["C", "H", "H", "H"], CH4.coords[:4] * [1, 1, 0])
H2O = XYZ(["O", "H", "H"], np.array([[0.0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]]))


def _species(sid, xyz, multiplicity=1):
    geo = Geometry(file=FileRef(path=f"{sid}.xyz", sha256="0"), fingerprint=sid,
                   symbols=tuple(xyz.symbols))
    return R.SpeciesRecord(species_id=sid, composition_id=sid.upper(), charge=0,
                           multiplicity=multiplicity, geometry=geo, source="input",
                           state_label=state_label(xyz.symbols, xyz.coords))


def _complex(a: XYZ, b: XYZ, gap=3.0) -> XYZ:  # b moved along z clear of a
    shift = a.coords[:, 2].max() - b.coords[:, 2].min() + gap
    return XYZ(a.symbols + b.symbols, np.vstack([a.coords, b.coords + [0, 0, shift]]))


def test_standard_states_per_molecule_count():
    assert th.standard_state_shift(1, 298.15, "1M") == pytest.approx(1.894, abs=1e-3)
    assert th.standard_state_shift(2, 298.15, "1M") == pytest.approx(
        2 * th.standard_state_shift(1, 298.15, "1M"))  # n molecules at one point
    assert th.standard_state_shift(2, 298.15, "1atm") == 0.0


@pytest.mark.parametrize(("G_ts", "G_P", "expected"), [
    (12.0, 3.0, 12.0),  # normal: the TS is the highest point
    (-0.5, 3.0, 3.0),  # dG_act < 0 (a submerged TS): the product state is the bottleneck
    (2.0, 5.0, 5.0),  # dG_act < dG_rxn: the TS lies below the product
    (None, 4.0, 4.0),  # a submerged TS dropped out: max(dG_rxn, 0)
    (None, -6.0, 0.0),
])
def test_a_unimolecular_chain_is_the_highest_point_above_the_reactant(G_ts, G_P, expected):
    for G_R in (0.0, -7.0):  # relative to G_R
        chain = [(G_R, True), *(() if G_ts is None else [(G_ts + G_R, False)]), (G_P + G_R, True)]
        assert th.chain_barrier(chain) == pytest.approx(expected)
        assert th.chain_barrier(chain) == pytest.approx(  # the old effective barrier
            max(G_R, G_P + G_R, G_R if G_ts is None else G_ts + G_R) - G_R)


def test_the_chain_refers_to_the_lowest_well_before_each_point():
    """S6: a complex 4.4 kcal/mol above its separated monomers (1 atm) reads the TS from the
    monomers; an SN2-like bound complex keeps its own zero; a TS below the complex is no
    bottleneck; a deep product well before an uphill dissociation counts from that well."""
    assert th.chain_barrier([(0.0, True), (4.4, True), (10.4, False), (-10.6, True)]
                            ) == pytest.approx(10.4)
    assert th.chain_barrier([(0.0, True), (-5.2, True), (9.3, False), (-5.2, True),
                             (0.0, True)]) == pytest.approx(14.5)
    assert th.chain_barrier([(0.0, True), (-3.0, True), (-1.0, False), (-20.0, True)]
                            ) == pytest.approx(2.0)
    assert th.chain_barrier([(0.0, True), (5.0, False), (-20.0, True), (-2.0, True)]
                            ) == pytest.approx(18.0)


def test_declared_monomers_keep_every_composition_of_one_formula():
    """oh_ch4 declares OH + CH4 and CH3 + H2O, both CH5O: neither replaces the other."""
    spc = [_species("oh", OH, 2), _species("ch4", CH4), _species("ch3", CH3, 2),
           _species("h2o", H2O)]
    comps = [{"id": "a", "components": {"oh": 1, "ch4": 1}},
             {"id": "b", "components": {"ch3": 1, "h2o": 1}},
             {"id": "c", "components": {"h2o": 2}}, {"id": "lost", "components": {"x": 1, "oh": 1}},
             {"id": "one", "components": {"oh": 1}}]
    monomers = th.declared_monomers(spc, [CompositionInput.model_validate(c) for c in comps])
    by_id = {s.species_id: (s.composition_id, s.state_label) for s in spc}
    assert monomers == {("CH5O", 0): [(by_id["oh"], by_id["ch4"]), (by_id["ch3"], by_id["h2o"])],
                        ("H4O2", 0): [(by_id["h2o"], by_id["h2o"])]}


def test_separated_states_need_the_fragments_to_be_the_declared_monomers():
    spc = [_species("oh", OH, 2), _species("ch4", CH4), _species("ch3", CH3, 2),
           _species("h2o", H2O)]
    monomers = th.declared_monomers(spc, [CompositionInput.model_validate({"id": "a", "components": {"oh": 1, "ch4": 1}}),
                                          CompositionInput.model_validate({"id": "b", "components": {"ch3": 1, "h2o": 1}})])
    by_id = {s.species_id: (s.composition_id, s.state_label) for s in spc}
    reactant, product = _complex(CH4, OH), _complex(CH3, H2O)
    assert th.separated_states(reactant, 0, monomers) == (by_id["oh"], by_id["ch4"])
    assert th.separated_states(product, 0, monomers) == (by_id["ch3"], by_id["h2o"])
    assert th.separated_states(reactant, -1, monomers) == ()  # another charge
    bonded = XYZ(product.symbols, np.vstack([CH3.coords, H2O.coords + [0, 0, 1.4]]))
    assert th.separated_states(bonded, 0, monomers) == ()  # C-O bonded: one fragment


def _rx(minima=("r", "p"), outcome=R.CaseOutcome.ELEMENTARY_STEP, source="discovery", **kw):
    saddle = R.SaddleClaim(saddle_calc="s", freq_calc="ts", imag_cm1=-900.0, energy_hartree=0.0)
    return R.ReactionRecord(reaction_id="rx", reactants=(), products=(), minima=minima,
                            endpoints=("a", "b"), source=source, outcome=outcome,
                            saddle=saddle, **kw)


STATES = {m: (m.upper(), m) for m in ("r", "p", "a", "b", "c", "d")}
SEP = {"r": (STATES["a"], STATES["b"]), "p": (STATES["c"], STATES["d"])}


def test_reaction_points_make_one_chain_with_its_separated_ends():
    pts = th.reaction_points(_rx(), STATES, SEP)
    assert pts.chain == (th.Point(SEP["r"], "r"), th.Point((STATES["r"],), "r"),
                         th.Point((), "ts"), th.Point((STATES["p"],), "p"),
                         th.Point(SEP["p"], "p"))
    assert [len(p.states) for p in pts.chain] == [2, 1, 0, 1, 2]  # molecules (a TS: 1)
    assert pts.own == (("r",), ("p",), ("ts",))
    assert (pts.separated, pts.complex) == (pts.chain[0], pts.chain[1])
    plain = th.reaction_points(_rx(), STATES, {})  # no declared monomers: [R, TS, P]
    assert [p.subject for p in plain.chain] == ["r", "ts", "p"] and plain.separated is None


def test_a_split_child_and_a_non_connected_outcome_have_no_chain_ends():
    child = th.reaction_points(_rx(source="split"), STATES, SEP)  # inside its parent's chain
    assert [p.subject for p in child.chain] == ["r", "ts", "p"]
    assert child.separated is None and child.complex is None
    flat = th.reaction_points(_rx(outcome=R.CaseOutcome.BARRIERLESS), STATES, SEP)
    assert flat.chain == () and flat.own == (("r",), ("p",), ())  # no TS read
    assert flat.complex == th.Point((STATES["r"],), "r")  # dG_assoc: an auxiliary point
    lost = th.reaction_points(_rx(minima=("r", "")), STATES, SEP)
    assert (lost.chain, lost.separated) == ((), None)


def test_an_association_reads_its_monomers_and_its_complex_unless_it_collapsed():
    """The monomers are the reactant side; the complex is R (a G_ref candidate) only where it
    is a minimum of its own, not the adduct it relaxed into."""
    assoc = _rx(monomers=("a", "b"), minima=("r", "p"))
    pts = th.reaction_points(assoc, STATES, {})
    assert pts.own == (("a", "b"), ("p",), ("ts",))
    assert [p.subject for p in pts.chain] == ["r", "r", "ts", "p"]
    assert pts.chain[0].states == (STATES["a"], STATES["b"]) and pts.complex == pts.chain[1]
    collapsed = th.reaction_points(assoc.model_copy(update={"minima": ("p", "p")}), STATES, {})
    assert [p.states for p in collapsed.chain] == [(STATES["a"], STATES["b"]), (), (STATES["p"],)]
    assert collapsed.complex is None


def test_a_chiral_structure_gains_minus_rt_ln2():  # m = 2 (HONO TS: 12.23 -> 11.82)
    assert th.chiral_G(298.15) * H2K == pytest.approx(-0.4107, abs=1e-4)
    assert th.chiral_G(298.15) == pytest.approx(-RT * math.log(2) / H2K)


def test_minimum_keeps_every_negative_mode_as_its_magnitude():  # noise (-4) and soft (-30)
    assert th.thermo_frequencies((500.0, -4.0, -30.0), saddle=False) == (4.0, 30.0, 500.0)


def test_saddle_drops_the_lowest_mode_and_flips_the_others():  # incl. an unannotated -8
    assert th.thermo_frequencies((600.0, -8.0, -700.0, -30.0), saddle=True) == (
        8.0, 30.0, 600.0)


def test_positive_or_empty_modes_pass_unchanged():
    assert th.thermo_frequencies((900.0, 300.0), saddle=False) == (300.0, 900.0)
    assert th.thermo_frequencies((), saddle=False) == th.thermo_frequencies((), saddle=True) == ()


def test_thermal_modes_follow_the_point_group_not_the_rank():
    """S18: H...H2 0.006 Å off the axis has rank 6 (3 modes, 2.27 kcal/mol in G); as C∞v it has
    3N - 5 = 4 modes at its own structure, a saddle 3."""
    x = np.array([[0.0, 0.0, -1.8], [0.0, 0.006, 0.0], [0.0, 0.0, 0.74]])
    h = 0.3 * np.eye(9)
    assert external_basis(["H"] * 3, x).shape[1] == 6
    assert len(th.thermal_modes(["H"] * 3, x, h, linear=True, saddle=False)) == 4
    assert len(th.thermal_modes(["H"] * 3, x, h, linear=True, saddle=True)) == 3
    assert len(th.thermal_modes(["H"] * 3, x, h, linear=False, saddle=False)) == 3
    assert th.thermal_modes(["H"], np.zeros((1, 3)), np.zeros((3, 3)), linear=False,
                            saddle=False) == ()


def _q_el(levels, T=298.15):  # exp(-G_el / RT)
    return math.exp(-th.electronic(levels, T) * H2K / (R_KCAL_MOL_K * T))


def test_electronic_levels_against_nist_hand_values():
    """O(3P) q_el(298) 5 + 3 e^(-158.3 hc/kT) + e^(-227.0 hc/kT) = 6.73; Cl(2P) E_SO
    -2/6 x 882.4 cm-1 = -0.84 kcal/mol; OH(2Pi) as oh_ch4.yaml declares it: q_el 3.02, E_SO
    -0.20 kcal/mol (the spin multiplet alone gave 2 and 0)."""
    oxygen, chlorine = atom_levels("O", 0, 3), atom_levels("Cl", 0, 2)
    assert _q_el(oxygen) == pytest.approx(6.73, abs=0.01)
    assert th.spin_orbit(oxygen) * H2K == pytest.approx(-0.223, abs=1e-3)
    assert th.spin_orbit(chlorine) * H2K == pytest.approx(-0.841, abs=1e-3)
    assert _q_el(chlorine) == pytest.approx(4.0 + 2.0 * math.exp(-882.4 / 207.22), abs=1e-3)
    system = load_system(Path(__file__).parents[3] / "configs" / "systems" / "oh_ch4.yaml")
    oh = next(s.electronic_levels for s in system.species if s.id == "oh")
    assert _q_el(oh) == pytest.approx(3.02, abs=0.01)
    assert th.spin_orbit(oh) * H2K == pytest.approx(-0.199, abs=1e-3)


def test_a_spin_multiplet_alone_is_the_old_electronic_entropy():
    assert th.electronic(((2, 0.0),), 298.15) == pytest.approx(-RT * math.log(2) / H2K, rel=1e-12)
    assert th.electronic(((1, 0.0),), 298.15) == 0.0 and th.spin_orbit(((3, 0.0),)) == 0.0


def test_declared_levels_hold_whole_spin_multiplets():
    levels = ((2, 0.0), (2, 139.21))
    assert SpeciesInput(id="oh", smiles="[OH]", multiplicity=2,
                        electronic_levels=levels).electronic_levels == levels
    with pytest.raises(ValueError, match="multiple of the multiplicity"):
        SpeciesInput(id="oh", smiles="[OH]", multiplicity=2, electronic_levels=((3, 0.0),))


def test_the_atom_table_holds_only_for_its_ground_term():
    assert atom_levels("Cl", -1, 1) is None and atom_levels("N", 0, 4) is None  # not split
    with pytest.raises(ValueError, match="ground term"):
        check_electronic_state(["O"], 0, 1)  # O(1D) is not the table's term
    check_electronic_state(["O"], 0, 3)


def test_pymsym_is_the_pinned_version():  # R14: another libmsym build may find another sigma
    pytest.importorskip("pymsym")
    assert importlib.metadata.version("pymsym") == "0.3.5"

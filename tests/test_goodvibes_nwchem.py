from pathlib import Path

from hfauto.backends.thermo.goodvibes import GoodVibesEngine, _prepare_goodvibes_input
from hfauto.core.hashing import sha256_file


def test_goodvibes_nwchem_copy_only_normalizes_parser_breaking_charge(
    tmp_path: Path,
) -> None:
    source = tmp_path / "nwchem.out"
    original = (
        "Northwest Computational Chemistry Package (NWChem)\n"
        "          Charge           :     0\n"
        "  charge          =   0.00\n"
        " Rotational Constants\n"
        " A= ********** cm-1  (********** K)\n"
        " B=   1.508433 cm-1  (  2.170251 K)\n"
        " C=   1.508433 cm-1  (  2.170251 K)\n"
        "         Total DFT energy =     -93.241363195833\n"
        " P.Frequency       55.46       90.80\n"
    )
    source.write_text(original, encoding="utf-8")
    original_hash = sha256_file(source)

    evidence = _prepare_goodvibes_input(source, tmp_path / "work")
    prepared = Path(evidence["prepared_path"])

    assert source.read_text(encoding="utf-8") == original
    assert sha256_file(source) == original_hash
    assert "charge 0\n" in prepared.read_text(encoding="utf-8")
    assert "A=   1.508433 cm-1  (  2.170251 K)" in prepared.read_text(
        encoding="utf-8"
    )
    assert "C*V symmetry detected" in prepared.read_text(encoding="utf-8")
    assert "Total DFT energy =     -93.241363195833" in prepared.read_text(
        encoding="utf-8"
    )
    assert evidence["compatibility_edits"] == 3
    assert evidence["compatibility_rules"] == [
        "nwchem_decimal_population_charge",
        "nwchem_linear_rotational_constant_mapping",
        "nwchem_linear_point_group_marker",
    ]
    assert evidence["scientific_values_changed"] is False


def test_goodvibes_v4_csv_fields_are_normalized(tmp_path: Path) -> None:
    csv_path = tmp_path / "Goodvibes.csv"
    csv_path.write_text(
        "file,name,zpe,qh_enthalpy,qh_gibbs_free_energy\n"
        "x.out,x,0.0108,-93.1524,-93.1772\n",
        encoding="utf-8",
    )

    row = GoodVibesEngine().parse_goodvibes_csv(csv_path).rows["x"]

    assert row["goodvibes_zpe_hartree"] == 0.0108
    assert row["goodvibes_enthalpy_hartree"] == -93.1524
    assert row["goodvibes_gibbs_hartree"] == -93.1772

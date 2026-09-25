from hfauto.stages.thermo import THERMAL_FREQUENCY_TASKS


def test_nwchem_saddle_frequency_is_a_ts_thermochemistry_source() -> None:
    assert "saddle_freq" in THERMAL_FREQUENCY_TASKS

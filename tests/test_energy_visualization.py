from hfauto_viz.renderers.plotly_energy import energy_points


def test_energy_profile_uses_generic_reaction_fields() -> None:
    labels, values = energy_points(
        {
            "delta_G_activation_standard_kcal_mol": 42.28,
            "delta_G_reaction_standard_kcal_mol": 12.74,
        }
    )

    assert labels == ["Reactants", "TS", "Products"]
    assert values == [0.0, 42.28, 12.74]


def test_energy_profile_does_not_invent_missing_states() -> None:
    labels, values = energy_points({})

    assert labels == ["Reactants"]
    assert values == [0.0]

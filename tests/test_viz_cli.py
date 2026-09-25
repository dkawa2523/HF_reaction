def test_visualization_cli_imports() -> None:
    from hfauto_viz.cli.main import app

    assert app.info.help == "Visualization companion CLI for completed hfauto runs"

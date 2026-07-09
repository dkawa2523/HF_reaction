from __future__ import annotations
from pathlib import Path
from hfauto_viz.data.loaders import RunData

def load_run(run_dir: str | Path) -> RunData:
    return RunData(run_dir)

__all__ = ['RunData', 'load_run']

from __future__ import annotations

from pathlib import Path
from typing import Any

from hfauto_viz.core.html import dataframe_to_html, write_html
from hfauto_viz.renderers.mol3d_3dmol import html_for_structure
from hfauto_viz.structure.bond_changes import annotation_shapes_for_frame
from hfauto_viz.structure.modes import (
    annotations_for_mode,
    approximate_proton_transfer_mode,
    mode_vectors_table,
    write_mode_vectors,
)
from hfauto_viz.structure.xyz import Frame


def write_imaginary_mode_view(
    ts: Frame | None,
    reaction_data: dict[str, Any] | None,
    out_path: str | Path,
    reactant: Frame | None = None,
    product: Frame | None = None,
    library_mode: str = "cdn",
    asset_prefix: str = "../../assets",
    title: str | None = None,
) -> dict[str, Path]:
    import pandas as pd

    out = Path(out_path); out.parent.mkdir(parents=True, exist_ok=True)
    vectors = approximate_proton_transfer_mode(ts, reaction_data, reactant=reactant, product=product)
    vec_paths = write_mode_vectors(vectors, out.with_name("imaginary_mode_vectors"))
    annotations = []
    annotations.extend(annotation_shapes_for_frame(ts, reaction_data))
    annotations.extend(annotations_for_mode(vectors))
    rows = mode_vectors_table(vectors)
    title = title or "Approximate imaginary-mode viewer"
    note = (
        "<div class='warn'><b>Mode-vector note:</b> This viewer can show true normal-mode arrows when parser data are available. "
        "For offline/fallback runs, arrows are an approximate proton-transfer direction based on q = r(B-H)-r(H-F), not a substitute for Hessian-derived normal modes.</div>"
    )
    viewer = html_for_structure(ts, "imag_mode_viewer", title, library_mode=library_mode, asset_prefix=asset_prefix, annotations=annotations)
    table = dataframe_to_html(pd.DataFrame(rows)) if rows else "<p>No mode vectors available.</p>"
    html = note + "<div class='card'>" + viewer + "</div><div class='card'><h2>Mode vector table</h2>" + table + "</div>"
    html_path = write_html(out, title, html)
    return {"imaginary_mode_html": html_path, **vec_paths}

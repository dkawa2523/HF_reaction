from __future__ import annotations

import csv
import json
from pathlib import Path

from hfauto_viz.core.html import esc, write_html


def write_publication_bundle(reaction_dir: str | Path, out_path: str | Path | None = None) -> dict[str, Path]:
    """Create a small index of static / reusable figure files for reports and papers."""
    reaction_dir = Path(reaction_dir)
    out_path = Path(out_path) if out_path else reaction_dir / "publication_figures.html"
    candidates = []
    for name, description in [
        ("energy_profile.svg", "Static relative free-energy diagram"),
        ("energy_profile.csv", "Energy profile source data"),
        ("reaction_coordinate.html", "Interactive reaction-coordinate diagnostics"),
        ("reaction_coordinate.csv", "Reaction-coordinate source data"),
        ("bond_change_annotations.csv", "Bond-change distances for B-H / H-F"),
        ("imaginary_mode_vectors.csv", "Imaginary-mode vector data"),
        ("path.xyz", "Reaction path frames"),
    ]:
        p = reaction_dir / name
        if p.exists():
            candidates.append({"file": name, "description": description, "bytes": p.stat().st_size})
    csv_path = reaction_dir / "publication_figure_index.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["file", "description", "bytes"])
        writer.writeheader(); writer.writerows(candidates)
    body = "<div class='card'><h1>Publication / review figure bundle</h1><p>Static and source-data files generated from the reaction dossier.</p><table class='data-table'><tr><th>File</th><th>Description</th><th>Bytes</th></tr>"
    for row in candidates:
        body += f"<tr><td><a href='{esc(row['file'])}'>{esc(row['file'])}</a></td><td>{esc(row['description'])}</td><td>{row['bytes']}</td></tr>"
    body += "</table></div>"
    html_path = write_html(out_path, "Publication figure bundle", body)
    json_path = reaction_dir / "publication_figure_index.json"
    json_path.write_text(json.dumps(candidates, ensure_ascii=False, indent=2), encoding="utf-8")
    return {"publication_figures_html": html_path, "publication_figures_csv": csv_path, "publication_figures_json": json_path}

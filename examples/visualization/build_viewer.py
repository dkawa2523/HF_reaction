"""Package saved scientific examples for offline browsing; no calculation is launched."""

from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def attach_profile_coordinates(data: dict) -> None:
    """Embed exported scan coordinates so file:// browsing does not need fetch()."""
    for example in data["examples"]:
        for profile in example.get("electronic_profiles", []):
            for point in profile["points"]:
                filename = point.get("geometry_file")
                if not filename:
                    continue
                source = HERE / "data" / filename
                xyz = source.read_text(encoding="utf-8")
                lines = xyz.splitlines()
                atoms = []
                for index, line in enumerate(lines[2 : 2 + int(lines[0])]):
                    symbol, x, y, z = line.split()[:4]
                    atoms.append({"index": index, "element": symbol,
                                  "x": float(x), "y": float(y), "z": float(z)})
                point.update(atoms=atoms, xyz=xyz)


def main() -> None:
    data = json.loads((HERE / "data" / "data.json").read_text(encoding="utf-8"))
    attach_profile_coordinates(data)
    references = json.loads((HERE / "references.json").read_text(encoding="utf-8"))
    library = json.loads((HERE / "vendor" / "source.json").read_text(encoding="utf-8"))
    payload = {"data": data, "references": references, "library": library}
    template = (HERE / "viewer.html").read_text(encoding="utf-8")
    embedded = json.dumps(payload, ensure_ascii=False).replace("</", "<\\/")
    output = template.replace("__EXAMPLE_PAYLOAD__", embedded)
    (HERE / "index.html").write_text(output, encoding="utf-8")
    print(f"Built {HERE / 'index.html'} ({len(data['examples'])} saved examples)")


if __name__ == "__main__":
    main()

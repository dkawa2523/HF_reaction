"""G20 (CH-26): eight ReaDuct NT2 TS frequency blocks (rotations and translations projected out)
each have exactly one imaginary mode; the legacy unprojected count rejected all of them."""

import json
import re

from hfauto.backends.readuct.worker import imaginary_count


def test_g20_blocks_have_one_projected_imaginary_mode(golden):
    text = golden.text("readuct/G20/pipeline.stdout.txt")
    blocks = [[float(v) for v in re.findall(r"^\s*\d+\s+([-+]\d+\.\d+)\s*$", block, re.MULTILINE)]
              for block in text.split("Vib. Frequencies:")[1:]]
    assert len(blocks) == 8
    assert [imaginary_count(freqs, 50.0) for freqs in blocks] == [1] * 8
    assert max(freqs[0] for freqs in blocks) == -239.2
    rows = golden.text("readuct/G20/reaction_discovery_attempts.jsonl").splitlines()
    assert all(json.loads(row)["imaginary_mode_count"] > 1 for row in rows)

"""Integrity of the golden fixture set: SOURCES.json, data files, size budget and naming rules."""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.golden

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
ENTRIES = json.loads((HERE / "SOURCES.json").read_text(encoding="utf-8"))
FILES = [f for e in ENTRIES for f in e["files"]]


def test_ids_are_g01_to_g24() -> None:
    assert [e["id"] for e in ENTRIES] == [f"G{i:02d}" for i in range(1, 25)]
    assert all(e["expect"] and isinstance(e["refs"], list) for e in ENTRIES)


def test_every_listed_file_exists_or_is_explained() -> None:
    for f in FILES:
        if f.get("missing"):
            assert f.get("reason"), f
        else:
            assert (DATA / f["golden"]).is_file(), f["golden"]
            assert re.fullmatch(r"[0-9a-f]{64}", f["source_sha256"]), f["golden"]
            assert f["method"] and f["source"]


def test_data_holds_only_listed_files_within_budget() -> None:
    stored = {p.relative_to(DATA).as_posix() for p in DATA.rglob("*") if p.is_file()}
    assert stored == {f["golden"] for f in FILES if not f.get("missing")}
    assert sum((DATA / name).stat().st_size for name in stored) <= 2 * 1024 * 1024


def test_names_avoid_gitignore_patterns() -> None:
    assert not list(HERE.rglob("__init__.py"))
    assert not [p for p in DATA.rglob("*") if p.is_dir() and p.name in {"runs", "cache"}]
    assert not list(DATA.rglob("*.log"))


def test_golden_fixture_reads_plain_and_gzipped_files(golden) -> None:
    assert "P.Frequency" in golden.text("nwchem/G01/nwchem.out")
    assert (
        golden.path("nwchem/G01/nwchem.out.gz")
        .read_text(encoding="utf-8")
        .count("Vibrational analysis via the FX method")
        == 4
    )
    assert golden.path("nwchem/G24/command_result.json") == DATA / "nwchem/G24/command_result.json"

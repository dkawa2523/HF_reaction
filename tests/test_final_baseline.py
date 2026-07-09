from pathlib import Path

from hfauto.core.artifact_types import STAGE_CONTRACTS


def test_current_docs_exist():
    root = Path('docs/current')
    required = [
        'index.md',
        'user_guide.md',
        'science_gate.md',
        'stage_io_reference.md',
        'output_reference.md',
        'production_checklist.md',
        'extension_points.md',
        'future_methods.md',
        'stage_future_extensions.md',
        'code_maintenance.md',
        'final_review.md',
    ]
    for name in required:
        path = root / name
        assert path.exists(), name
        assert path.read_text(encoding='utf-8').strip(), name


def test_stage_contracts_are_compact_and_cover_core_flow():
    for stage in [
        'ingest', 'enrich', 'detect-sites', 'conformers', 'build-hf', 'preopt',
        'dft-minima', 'ts-search', 'irc', 'sp', 'thermo', 'kinetics', 'rank', 'viz', 'ops'
    ]:
        assert stage in STAGE_CONTRACTS
        assert STAGE_CONTRACTS[stage]['out']


def test_rank_stage_no_duplicate_scientific_artifact_registration():
    text = Path('hfauto/stages/rank.py').read_text(encoding='utf-8')
    assert text.count('("rank_scientific", "scientific", scientific_path, len(scientific))') == 1

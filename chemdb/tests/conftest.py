import sys
from pathlib import Path

import pytest


REPORT_TEST_PREFIXES = (
    "test_heatmap",
    "test_report",
    "test_run_report",
)
INTEGRATION_TEST_FILES = {
    "test_ligand_run_modes.py",
    "test_main_full_run.py",
    "test_modes_integration.py",
}
EXTERNAL_TOOL_TEST_FILES = {
    "test_protein_prep_refactor_receptor_protonation.py",
    "test_protonation_reduce.py",
    "test_receptor_prep_unit.py",
    "test_tool_runners_commands.py",
}
FAST_TEST_FILES = {
    "test_cli_doctor.py",
    "test_cli_help.py",
    "test_cli_init_new_run.py",
    "test_config_inline_comments_parsing.py",
    "test_config_reference_docs.py",
    "test_path_router.py",
    "test_record_data_csv.py",
}
SLOW_TEST_FILES = {
    "test_main_full_run.py",
    "test_modes_integration.py",
}


def pytest_configure(config):
    # Calculate repo root relative to this conftest file
    # chemdb/tests/conftest.py -> ../.. -> repo_root
    repo_root = Path(__file__).resolve().parents[2]
    src_root = repo_root / "src"

    # Ensure repo_root is in sys.path
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    if str(src_root) not in sys.path:
        # Insert after repo_root to avoid shadowing existing modules during migration
        try:
            idx = sys.path.index(str(repo_root))
            sys.path.insert(idx + 1, str(src_root))
        except ValueError:
            sys.path.append(str(src_root))


def pytest_collection_modifyitems(config, items):
    for item in items:
        filename = Path(str(item.path)).name
        if filename.startswith(REPORT_TEST_PREFIXES):
            item.add_marker(pytest.mark.report)
        if filename in INTEGRATION_TEST_FILES:
            item.add_marker(pytest.mark.integration)
        if filename in EXTERNAL_TOOL_TEST_FILES:
            item.add_marker(pytest.mark.external_tool)
        if filename in FAST_TEST_FILES:
            item.add_marker(pytest.mark.fast)
        if filename in SLOW_TEST_FILES:
            item.add_marker(pytest.mark.slow)

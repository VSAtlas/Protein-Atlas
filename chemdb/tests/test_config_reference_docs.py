import re
from pathlib import Path


def test_config_reference_documents_template_keys() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    docs_text = (repo_root / "docs" / "config_reference.md").read_text(encoding="utf-8")
    template_text = "\n".join(
        [
            (repo_root / "config.example.txt").read_text(encoding="utf-8"),
            (repo_root / "config.full.example.txt").read_text(encoding="utf-8"),
        ]
    )
    expected_keys = {
        "OVERALL_DIR",
        "MGLTOOLS_PATH",
        "ADFRSUITE_BIN",
        "P2RANK_PATH",
        "PHENIX_DIR",
        "PHENIX_LIB_PATH",
        "SCORCH_ENV_PREFIX",
        "SCORCH_ENV",
        "MMGBSA_AMBERTOOLS_PREFIX",
        "AMBERTOOLS_PREFIX",
        "TEST_LIBRARY_MAP",
        "TOOL_VERIFY_ON_START",
    }
    documented = set(re.findall(r"`([A-Z][A-Z0-9_]+)`", docs_text))
    assert expected_keys <= documented
    for key in expected_keys:
        assert re.search(rf"^\s*#?\s*{key}\s*=", template_text, re.M), key

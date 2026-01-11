import sys
from pathlib import Path


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

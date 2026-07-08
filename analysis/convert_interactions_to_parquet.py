"""Compatibility wrapper for moved module."""

import sys
import warnings
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))
if str(_REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT / "src"))

from analysis.cli import convert_interactions_to_parquet as _impl  # noqa: E402

warnings.warn(
    f"{__name__} is deprecated; use analysis.cli.convert_interactions_to_parquet. Removal after 2 release cycles.",
    DeprecationWarning,
    stacklevel=2,
)

if __name__ == "__main__":
    raise SystemExit(_impl.main())

sys.modules[__name__] = _impl

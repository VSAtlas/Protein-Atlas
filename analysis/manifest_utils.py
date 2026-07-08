"""Compatibility wrapper for moved module."""

import sys
import warnings

from analysis.reporting import manifest_utils as _impl

warnings.warn(
    f"{__name__} is deprecated; use analysis.reporting.manifest_utils. Removal after 2 release cycles.",
    DeprecationWarning,
    stacklevel=2,
)

sys.modules[__name__] = _impl

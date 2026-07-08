"""Compatibility wrapper for moved module."""

import sys
import warnings

from analysis.reporting import fda_name_map as _impl

warnings.warn(
    f"{__name__} is deprecated; use analysis.reporting.fda_name_map. Removal after 2 release cycles.",
    DeprecationWarning,
    stacklevel=2,
)

sys.modules[__name__] = _impl

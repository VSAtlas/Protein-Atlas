"""Compatibility wrapper for moved module."""

import sys
import warnings

from analysis.dud_eval_core import reporting as _impl

warnings.warn(
    f"{__name__} is deprecated; use analysis.dud_eval_core.reporting. Removal after 2 release cycles.",
    DeprecationWarning,
    stacklevel=2,
)

sys.modules[__name__] = _impl

# -*- coding: utf-8 -*-
# ruff: noqa: F401
import csv
import datetime
import gzip
import html
import json
import logging
import math
import os
import re
import shutil
import sys
import sys as _sys
import yaml  # type: ignore[import-untyped]
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from analysis.reporting.fda_name_map import (
    resolve_ligand_display_name,
    resolve_mapping_csv_path,
    try_load_fda_index,
)
from analysis.dud_eval_core.discovery import derive_target_name
from analysis.reporting.fdr_methods import compute_score_fdr
from analysis.reporting.heatmap_html import render_interactive_heatmap_html
from analysis.reporting.manifest_utils import extract_pocket, load_run_manifest
from analysis.reporting.publication_exports import write_publication_exports
from analysis.reporting.report_config import (
    read_config_key,
    report_config_candidates,
    resolve_config_key,
)
from analysis.reporting.statistical_formula_html import render_statistical_formula_block
from analysis.reporting.test_mode_tokens import (
    TEST_MODE_OFF_VALUES,
    TEST_MODE_ON_VALUES,
    resolve_test_mode_tokens,
)
from analysis.reporting.value_utils import (
    as_float,
    as_report_row_bool,
    clean_report_text,
    median_floats,
    strip_quotes,
)
from analysis.statistics import ranked_binary_metrics
from analysis.target_ids import build_target_id, parse_target_id, pdb_id_from_target_id
from analysis import pathway_resolver
from config.output_paths import run_output_dir
from chemdb.decoy_control import (
    ligand_filename_from_row,
    row_has_explicit_decoy,
    row_is_control,
    row_is_decoy_by_role,
)

COMPONENT = "[run-report]"
DECOY_PREFIX_KEY = "DECOY_PREFIX"
DUD_PREFIX_KEY = "DUD_PREFIX"
DECOY_PREFIX_DEFAULT = "dud"

try:
    from pandas.errors import EmptyDataError as _PandasEmptyDataError
    from pandas.errors import ParserError as _PandasParserError
except ImportError:  # pragma: no cover
    class _PandasCsvReadError(Exception):
        """Placeholder when pandas exception types are unavailable."""

    _PandasEmptyDataError = _PandasCsvReadError
    _PandasParserError = _PandasCsvReadError

_REPORT_IO_ERRORS: tuple[type[BaseException], ...] = (
    OSError,
    UnicodeDecodeError,
    csv.Error,
    ValueError,
    TypeError,
    _PandasEmptyDataError,
    _PandasParserError,
)

_REPORT_PUBLICATION_PIPELINE_ERRORS: tuple[type[BaseException], ...] = _REPORT_IO_ERRORS + (
    MemoryError,
    KeyError,
)

_REPORT_HTML_BUILD_ERRORS: tuple[type[BaseException], ...] = _REPORT_IO_ERRORS + (
    KeyError,
    AttributeError,
    MemoryError,
)

_REPORT_MAIN_ERRORS: tuple[type[BaseException], ...] = _REPORT_HTML_BUILD_ERRORS + (
    RuntimeError,
    json.JSONDecodeError,
    yaml.YAMLError,
)

MIN_DECOYS_FOR_FDR = 200
MIN_UNIQUE_DECOY_SCORES = 10
REPORT_HTML_WARN_BYTES_DEFAULT = 10 * 1024 * 1024
_INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_RDK_PLACEHOLDER_RE = re.compile(r"rdk[_-]?\d+", re.IGNORECASE)
_VINA_RESULT_RE = re.compile(r"REMARK\s+VINA\s+RESULT:\s*([-+0-9.eE]+)")
_TARGET_NAME_CACHE: Dict[str, str] = {}

__all__ = sorted(
    name for name in _sys.modules[__name__].__dict__ if not name.startswith("__")
)

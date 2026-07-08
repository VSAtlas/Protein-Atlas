# -*- coding: utf-8 -*-
import csv
import logging
from pathlib import Path
from typing import Dict, List, Optional

from analysis.reporting.master_schema.constants import COMPONENT, _SCHEMA_CSV_ERRORS

def _configure_logging(verbose: bool) -> logging.Logger:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(message)s")
    return logging.getLogger("master-export")


def _try_read_utf8_lines_logged(
    path: Path, logger: logging.Logger, *, action: str
) -> Optional[List[str]]:
    """Return ``None`` on read failure so callers can distinguish from an empty file."""
    try:
        with path.open("r", encoding="utf-8") as handle:
            return handle.readlines()
    except OSError as exc:
        logger.warning(
            "%s action=%s status=read_failed path=%s error=%s",
            COMPONENT,
            action,
            path,
            exc,
        )
        return None


def _read_csv_rows(path: Path) -> List[Dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            return list(reader)
    except _SCHEMA_CSV_ERRORS:
        return []

# -*- coding: utf-8 -*-
"""Run report generation package (split from run_report_core)."""

from analysis.reporting.run_report.build_report import build_report
from analysis.reporting.run_report.config_resolve import _configure_logging
from analysis.reporting.run_report.constants import (
    COMPONENT,
    DECOY_PREFIX_DEFAULT,
    DECOY_PREFIX_KEY,
    DUD_PREFIX_KEY,
    MIN_DECOYS_FOR_FDR,
    MIN_UNIQUE_DECOY_SCORES,
    REPORT_HTML_WARN_BYTES_DEFAULT,
)
from analysis.reporting.run_report.export_io import write_yaml
from analysis.reporting.run_report.heatmap_csv import _write_heatmap_input_csv
from analysis.reporting.run_report.html_report import _write_html_report
from analysis.reporting.run_report.main import main
from analysis.reporting.run_report.pathway_reports import _write_html_reports
from analysis.reporting.run_report.target_stats import _load_scorch_stats

__all__ = [
    "COMPONENT",
    "DECOY_PREFIX_DEFAULT",
    "DECOY_PREFIX_KEY",
    "DUD_PREFIX_KEY",
    "MIN_DECOYS_FOR_FDR",
    "MIN_UNIQUE_DECOY_SCORES",
    "REPORT_HTML_WARN_BYTES_DEFAULT",
    "_configure_logging",
    "_load_scorch_stats",
    "_write_heatmap_input_csv",
    "_write_html_report",
    "_write_html_reports",
    "build_report",
    "main",
    "write_yaml",
]

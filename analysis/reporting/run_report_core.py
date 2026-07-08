# ruff: noqa: F403, F405
# -*- coding: utf-8 -*-
"""Backward-compatible re-export shim for analysis.reporting.run_report."""

from analysis.reporting.run_report import *  # noqa: F401,F403
from analysis.reporting.heatmap_html import render_interactive_heatmap_html
from analysis.reporting.run_report import html_report as _html_report_module


def _write_html_report(*args, **kwargs):
    original_renderer = _html_report_module.render_interactive_heatmap_html
    _html_report_module.render_interactive_heatmap_html = render_interactive_heatmap_html
    try:
        return _html_report_module._write_html_report(*args, **kwargs)
    finally:
        _html_report_module.render_interactive_heatmap_html = original_renderer

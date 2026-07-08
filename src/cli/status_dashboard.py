from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from analysis.cli.summarize_run_errors import build_error_summary
from analysis.reporting.manifest_utils import load_run_manifest
from cli.status_dashboard_eta import (
    ETA_HISTORY_NAME,
    estimate_eta,
    eta_history_summary,
    refresh_eta_history,
)
from cli.status_dashboard_explain import explain_status, next_steps
from cli.status_dashboard_panels import (
    format_runs_table,
    format_status_dashboard,
    render_status_html,
    write_status_html,
)
from cli.status_dashboard_progress import (
    distributed_progress,
    lifecycle,
    overall_progress,
    report_only_lifecycle,
    report_only_progress,
    scorch_failure_summary,
    stage_progress,
)
from cli.status_dashboard_slurm import parse_sacct_output, parse_squeue_output, slurm_status
from cli.status_dashboard_utils import path_age_seconds, run_candidates
from config.output_paths import output_root, run_output_dir

# Re-export public API for backward compatibility.
__all__ = (
    "ETA_HISTORY_NAME",
    "build_runs_index",
    "build_status_dashboard",
    "default_status_html_path",
    "explain_status",
    "format_runs_table",
    "format_status_dashboard",
    "parse_sacct_output",
    "parse_squeue_output",
    "render_status_html",
    "write_status_html",
)


def build_status_dashboard(
    root: Path,
    run_id: str,
    *,
    include_slurm: bool = True,
    include_errors: bool = False,
    deep_errors: bool = False,
) -> dict[str, Any]:
    manifest, manifest_path = load_run_manifest(root, run_id)
    paths, existing_paths, distributed_dir = _collect_status_paths(root, run_id, manifest_path)
    status = _base_status_payload(
        run_id,
        manifest,
        manifest_path,
        paths,
        existing_paths,
        distributed_dir,
    )
    _apply_report_only_fallback(status, manifest, paths["report_html"])
    status["scorch_failures"] = scorch_failure_summary(manifest)
    if include_slurm:
        status["slurm"] = slurm_status(manifest)
    eta_history = refresh_eta_history(root)
    status["eta_history"] = eta_history_summary(eta_history)
    status["eta"] = estimate_eta(status, eta_history)
    if include_errors:
        status["errors"] = build_error_summary(
            root=root,
            run_id=run_id,
            deep=deep_errors,
            top=5,
            max_examples=2,
        )
    status["next"] = next_steps(status)
    status["explanations"] = explain_status(status)
    return status


def _collect_status_paths(
    root: Path,
    run_id: str,
    manifest_path: Path | None,
) -> tuple[dict[str, Path | None], dict[str, str], Path]:
    data_dir = run_output_dir(root, "data", run_id)
    manifests_dir = run_output_dir(root, "manifests", run_id)
    paths: dict[str, Path | None] = {
        "manifest": manifest_path,
        "report_html": data_dir / "report.html",
        "report_yaml": data_dir / "report.yaml",
        "heatmap_csv": data_dir / "heatmap_input.csv",
        "docked": run_output_dir(root, "docked", run_id),
        "post_docked": run_output_dir(root, "post_docked", run_id),
        "logs": output_root(root, "logs"),
        "configs": run_output_dir(root, "configs", run_id),
        "distributed": manifests_dir / "distributed",
        "scheduler_summary": manifests_dir / "scheduler_summary.json",
        "run_efficiency": manifests_dir / "run_efficiency.json",
    }
    existing_paths = {
        key: str(path)
        for key, path in paths.items()
        if path is not None and Path(path).exists()
    }
    return paths, existing_paths, manifests_dir / "distributed"


def _base_status_payload(
    run_id: str,
    manifest: Mapping[str, Any] | None,
    manifest_path: Path | None,
    paths: dict[str, Path | None],
    existing_paths: dict[str, str],
    distributed_dir: Path,
) -> dict[str, Any]:
    return {
        "run_id": str(run_id),
        "found_any": bool(existing_paths or manifest),
        "manifest": str(manifest_path) if manifest_path else "",
        "lifecycle": lifecycle(manifest),
        "progress": overall_progress(manifest),
        "stages": stage_progress(manifest),
        "distributed": distributed_progress(distributed_dir),
        "slurm": {},
        "eta": {},
        "errors": {},
        "paths": {key: str(value) for key, value in paths.items() if value is not None},
        "existing_paths": existing_paths,
    }


def _apply_report_only_fallback(
    status: dict[str, Any],
    manifest: Mapping[str, Any] | None,
    report_html: Path | None,
) -> None:
    if report_html is None or not report_html.exists():
        return
    if isinstance(manifest, Mapping) and manifest:
        return
    status["lifecycle"] = report_only_lifecycle(report_html)
    status["progress"] = report_only_progress(report_html)
    status["stages"] = {}


def build_runs_index(root: Path, *, limit: int = 20) -> list[dict[str, Any]]:
    rows = []
    for _mtime, run_id in run_candidates(root)[: max(1, int(limit))]:
        manifest, _manifest_path = load_run_manifest(root, run_id)
        lifecycle_data = lifecycle(manifest)
        progress = overall_progress(manifest)
        report_html = run_output_dir(root, "data", run_id) / "report.html"
        manifest_path = run_output_dir(root, "manifests", run_id) / "run_manifest.yaml"
        if (not isinstance(manifest, Mapping) or not manifest) and report_html.exists():
            lifecycle_data = report_only_lifecycle(report_html)
            progress = report_only_progress(report_html)
        rows.append(
            {
                "run_id": run_id,
                "status": lifecycle_data.get("status", "unknown"),
                "percent_done": progress.get("percent_done", 0.0),
                "completed": progress.get("completed", 0),
                "failed": progress.get("failed", 0),
                "running": progress.get("running", 0),
                "scheduled": progress.get("scheduled", 0),
                "age_sec": path_age_seconds(report_html if report_html.exists() else manifest_path),
                "report_html": str(report_html) if report_html.exists() else "",
                "manifest": str(manifest_path) if manifest_path.exists() else "",
            }
        )
    return rows


def default_status_html_path(root: Path, run_id: str) -> Path:
    return output_root(root, "data") / str(run_id) / "status.html"

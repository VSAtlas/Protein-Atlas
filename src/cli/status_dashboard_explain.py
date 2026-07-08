from __future__ import annotations

from typing import Any, Mapping

from cli.status_dashboard_utils import STAGE_ORDER, as_int


def explain_status(status: Mapping[str, Any]) -> list[str]:
    run_id = str(status.get("run_id") or "<RUN_ID>")
    if not status.get("found_any"):
        return [f"No run artifacts were found for {run_id}. Check the run id with `atlas runs`."]
    explanations: list[str] = []
    _append_slurm_explanations(explanations, status.get("slurm"))
    _append_stage_explanations(explanations, status.get("stages"))
    _append_scorch_explanations(explanations, status.get("scorch_failures"), run_id)
    _append_error_explanations(explanations, status.get("errors"), run_id)
    _append_report_explanations(explanations, status.get("existing_paths"), run_id)
    if not explanations:
        explanations.append("No immediate operator action is required from the dashboard state.")
    return explanations


def next_steps(status: Mapping[str, Any]) -> list[str]:
    paths = status.get("existing_paths")
    paths = paths if isinstance(paths, Mapping) else {}
    progress = status.get("progress")
    progress = progress if isinstance(progress, Mapping) else {}
    run_id = str(status.get("run_id") or "<RUN_ID>")
    steps: list[str] = []
    if paths.get("report_html"):
        steps.append(f"Open report: {paths['report_html']}")
    elif paths.get("heatmap_csv"):
        steps.append(f"Generate report: atlas report {run_id}")
    elif paths.get("manifest") or paths.get("docked") or paths.get("post_docked"):
        steps.append(f"Generate report: atlas report {run_id} --vina")
    scorch = status.get("scorch_failures")
    if isinstance(scorch, Mapping) and (as_int(scorch.get("failed")) or 0) > 0:
        steps.append(f"Inspect Scorch failures: atlas debug {run_id} --deep")
    if (as_int(progress.get("failed")) or 0) > 0:
        steps.append(f"Inspect root causes: atlas status {run_id} --errors")
    if not steps:
        steps.append("Create a demo: atlas demo")
    return steps


def _append_slurm_explanations(explanations: list[str], slurm: Any) -> None:
    if not isinstance(slurm, Mapping):
        return
    for job in slurm.get("jobs") or []:
        if not isinstance(job, Mapping):
            continue
        state = str(job.get("state") or "").upper()
        reason = str(job.get("reason_or_nodelist") or "").strip()
        if state == "PENDING":
            explanations.append(
                f"Slurm job {job.get('job_id')} is pending"
                + (f" because {reason.strip('()')}." if reason else ".")
            )
        elif state in {"FAILED", "OUT_OF_MEMORY", "TIMEOUT", "CANCELLED", "NODE_FAIL"}:
            detail = job.get("exit_code") or reason
            explanations.append(
                f"Slurm job {job.get('job_id')} ended with state {state}"
                + (f" ({detail})." if detail else ".")
            )


def _append_stage_explanations(explanations: list[str], stages: Any) -> None:
    if not isinstance(stages, Mapping):
        return
    for stage in STAGE_ORDER:
        row = stages.get(stage)
        if not isinstance(row, Mapping):
            continue
        failed = as_int(row.get("failed")) or 0
        running = as_int(row.get("running")) or 0
        pending = as_int(row.get("pending")) or 0
        if failed:
            examples = row.get("failed_examples") or []
            suffix = f" Examples: {', '.join(str(x) for x in examples[:3])}." if examples else ""
            explanations.append(f"{stage} has {failed} failed item(s).{suffix}")
        elif running:
            explanations.append(f"{stage} is currently running for {running} item(s).")
        elif pending:
            explanations.append(f"{stage} has {pending} pending item(s).")


def _append_scorch_explanations(
    explanations: list[str], scorch: Any, run_id: str
) -> None:
    if not isinstance(scorch, Mapping):
        return
    failed = as_int(scorch.get("failed")) or 0
    if failed <= 0:
        return
    errors = scorch.get("errors")
    top_error = ""
    if isinstance(errors, Mapping) and errors:
        top_error = str(next(iter(errors)))
    returncodes = scorch.get("returncodes")
    top_returncode = ""
    if isinstance(returncodes, Mapping) and returncodes:
        top_returncode = str(next(iter(returncodes)))
    examples = scorch.get("examples") or []
    example_text = ""
    if examples:
        combos = [str(row.get("combo")) for row in examples[:3] if isinstance(row, Mapping)]
        if combos:
            example_text = f" Examples: {', '.join(combos)}."
    detail = top_error or "postprocessing_failed"
    if top_returncode and f"returncode={top_returncode}" not in detail:
        detail = f"{detail}; returncode={top_returncode}"
    explanations.append(
        f"SCORCH postprocessing failed for {failed} combo(s) ({detail})."
        + example_text
        + f" Run atlas debug {run_id} --deep for the manifest-backed failure summary."
    )


def _append_error_explanations(
    explanations: list[str], errors: Any, run_id: str
) -> None:
    if not isinstance(errors, Mapping):
        return
    root_counts = errors.get("root_cause_counts")
    if isinstance(root_counts, Mapping) and root_counts:
        top_cause = next(iter(root_counts))
        explanations.append(
            f"The most common error bucket is {top_cause}; run `atlas status {run_id} --errors --deep` if the examples are not enough."
        )


def _append_report_explanations(
    explanations: list[str], paths: Any, run_id: str
) -> None:
    paths = paths if isinstance(paths, Mapping) else {}
    if not paths.get("report_html") and (paths.get("manifest") or paths.get("docked")):
        explanations.append(f"No HTML report exists yet; run `atlas report {run_id} --vina`.")

from __future__ import annotations

import datetime as _dt
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from cli.status_dashboard_utils import (
    STAGE_ORDER,
    as_int,
    esc,
    eta_label,
    format_seconds,
    mapping,
)


def format_status_dashboard(status: Mapping[str, Any], *, explain: bool = False) -> str:
    lines = [f"run_id: {status.get('run_id')}"]
    if not status.get("found_any"):
        return "\n".join([*lines, "status: not found"])
    lifecycle = mapping(status.get("lifecycle"))
    lines.append(f"status: {lifecycle.get('status', 'unknown')}")
    progress = mapping(status.get("progress"))
    lines.append(
        "progress: {done}/{scheduled} done ({pct}%) completed={completed} failed={failed} running={running}".format(
            done=progress.get("done", 0),
            scheduled=progress.get("scheduled", 0),
            pct=progress.get("percent_done", 0.0),
            completed=progress.get("completed", 0),
            failed=progress.get("failed", 0),
            running=progress.get("running", 0),
        )
    )
    append_slurm(lines, status.get("slurm"))
    append_eta(lines, status.get("eta"))
    append_stage_table(lines, status.get("stages"))
    append_scorch_failures(lines, status.get("scorch_failures"))
    append_distributed(lines, status.get("distributed"))
    append_errors(lines, status.get("errors"))
    if explain:
        append_explanations(lines, status.get("explanations"))
    failed_examples = progress.get("failed_examples")
    if failed_examples:
        lines.append("failed_examples: " + ", ".join(str(x) for x in failed_examples))
    existing = status.get("existing_paths")
    if isinstance(existing, Mapping) and existing:
        lines.append("outputs:")
        for key in sorted(existing):
            lines.append(f"  {key}: {existing[key]}")
    for step in status.get("next") or []:
        lines.append(f"next: {step}")
    return "\n".join(lines)


def append_scorch_failures(lines: list[str], scorch: Any) -> None:
    if not isinstance(scorch, Mapping):
        return
    failed = as_int(scorch.get("failed")) or 0
    if failed <= 0:
        return
    returncodes = mapping(scorch.get("returncodes"))
    errors = mapping(scorch.get("errors"))
    rc_text = ", ".join(f"{key}:{value}" for key, value in returncodes.items()) or "none"
    top_error = next(iter(errors), "postprocessing_failed")
    lines.append(
        f"scorch_failures: failed={failed} top_error={top_error} returncodes={rc_text}"
    )
    for row in (scorch.get("examples") or [])[:3]:
        if not isinstance(row, Mapping):
            continue
        lines.append(
            "  scorch_failure: combo={combo} error={error} returncode={returncode} mode={mode} phase={phase} device={device} wall={wall}".format(
                combo=row.get("combo", ""),
                error=row.get("error", ""),
                returncode=row.get("returncode", ""),
                mode=row.get("mode", ""),
                phase=row.get("phase", ""),
                device=row.get("device", ""),
                wall=row.get("wall_time_sec", ""),
            )
        )


def format_runs_table(rows: Sequence[Mapping[str, Any]]) -> str:
    if not rows:
        return "No runs found. Start with `atlas demo` or `atlas --pdb <ID> --fast`."
    lines = ["run_id              status      progress  failed  age      report"]
    for row in rows:
        report = "yes" if row.get("report_html") else "no"
        lines.append(
            "{run_id:<19} {status:<11} {pct:>6.1f}% {failed:>7} {age:<8} {report}".format(
                run_id=str(row.get("run_id", ""))[:19],
                status=str(row.get("status", "unknown"))[:11],
                pct=float(row.get("percent_done") or 0.0),
                failed=int(row.get("failed") or 0),
                age=format_seconds(float(row.get("age_sec") or 0.0)),
                report=report,
            )
        )
    return "\n".join(lines)


def render_status_html(status: Mapping[str, Any]) -> str:
    run_id = str(status.get("run_id") or "")
    lifecycle = mapping(status.get("lifecycle"))
    progress = mapping(status.get("progress"))
    eta = mapping(status.get("eta"))
    slurm = mapping(status.get("slurm"))
    errors = mapping(status.get("errors"))
    title = f"Atlas Status {run_id}".strip()
    status_text = str(lifecycle.get("status") or "unknown")
    pct = float(progress.get("percent_done") or 0.0)
    payload = json.dumps(status, indent=2, sort_keys=True, default=str)
    generated_at = _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat()
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{esc(title)}</title>
  <style>
    :root {{
      color-scheme: light;
      --bg: #f7f7f4;
      --fg: #18201d;
      --muted: #5b635f;
      --line: #d8ddd7;
      --panel: #ffffff;
      --accent: #1f7a5a;
      --warn: #a35d00;
      --bad: #b3261e;
    }}
    body {{
      margin: 0;
      background: var(--bg);
      color: var(--fg);
      font: 14px/1.45 system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    }}
    header {{
      padding: 24px clamp(16px, 4vw, 48px) 16px;
      border-bottom: 1px solid var(--line);
      background: var(--panel);
    }}
    main {{
      padding: 20px clamp(16px, 4vw, 48px) 36px;
      max-width: 1180px;
    }}
    h1 {{ margin: 0 0 8px; font-size: 24px; letter-spacing: 0; }}
    h2 {{ margin: 0 0 10px; font-size: 16px; letter-spacing: 0; }}
    .meta {{ color: var(--muted); }}
    .grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 12px;
      margin-bottom: 16px;
    }}
    section, .metric {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 14px;
    }}
    .metric strong {{ display: block; font-size: 22px; margin-top: 4px; }}
    .bar {{ height: 10px; background: #e6ebe7; border-radius: 999px; overflow: hidden; }}
    .bar span {{ display: block; height: 100%; width: {max(0.0, min(100.0, pct)):.1f}%; background: var(--accent); }}
    table {{ width: 100%; border-collapse: collapse; }}
    th, td {{ text-align: left; padding: 8px; border-bottom: 1px solid var(--line); vertical-align: top; }}
    th {{ color: var(--muted); font-weight: 600; }}
    code, pre {{ font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; }}
    pre {{ overflow: auto; background: #111814; color: #e8f3ec; padding: 12px; border-radius: 8px; }}
    .stack {{ display: grid; gap: 12px; }}
    .pill {{ display: inline-block; border: 1px solid var(--line); border-radius: 999px; padding: 2px 8px; }}
    .bad {{ color: var(--bad); }}
    .warn {{ color: var(--warn); }}
  </style>
</head>
<body>
  <header>
    <h1>{esc(title)}</h1>
    <div class="meta">Generated {generated_at} from <code>atlas status --json</code></div>
  </header>
  <main class="stack">
    <div class="grid">
      <div class="metric"><span>Status</span><strong>{esc(status_text)}</strong></div>
      <div class="metric"><span>Progress</span><strong>{pct:.1f}%</strong><div class="bar"><span></span></div></div>
      <div class="metric"><span>ETA</span><strong>{esc(eta_label(eta))}</strong><span class="meta">{esc(str(eta.get("basis") or "unknown"))}</span></div>
      <div class="metric"><span>Errors</span><strong>{int(errors.get("raw_hits") or 0)}</strong><span class="meta">{int(errors.get("unique_signatures") or 0)} unique</span></div>
    </div>
    {html_explanations(status.get("explanations"))}
    {html_stage_table(status.get("stages"))}
    {html_slurm(slurm)}
    {html_errors(errors)}
    {html_paths(status.get("existing_paths"))}
    <section>
      <h2>Status JSON</h2>
      <pre>{esc(payload)}</pre>
    </section>
  </main>
</body>
</html>
"""


def write_status_html(status: Mapping[str, Any], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_status_html(status), encoding="utf-8")
    return path


def append_explanations(lines: list[str], explanations: Any) -> None:
    if not isinstance(explanations, list) or not explanations:
        return
    lines.append("explain:")
    for item in explanations:
        lines.append(f"  - {item}")


def append_stage_table(lines: list[str], stages: Any) -> None:
    if not isinstance(stages, Mapping) or not stages:
        return
    if not any(
        (as_int(row.get("total")) or 0) > 0 for row in stages.values() if isinstance(row, Mapping)
    ):
        return
    lines.append("stages:")
    for stage in STAGE_ORDER:
        row = stages.get(stage)
        if not isinstance(row, Mapping):
            continue
        lines.append(
            "  {stage}: {pct}% done completed={completed} failed={failed} running={running} pending={pending}".format(
                stage=stage,
                pct=row.get("percent_done", 0.0),
                completed=row.get("completed", 0),
                failed=row.get("failed", 0),
                running=row.get("running", 0),
                pending=row.get("pending", 0),
            )
        )


def append_distributed(lines: list[str], distributed: Any) -> None:
    if not isinstance(distributed, Mapping) or not distributed.get("enabled"):
        return
    lines.append(
        "distributed: chunks={done}/{planned} ({pct}%) failed={failed} claims={claims}".format(
            done=distributed.get("completed_chunks", 0),
            planned=distributed.get("planned_chunks", 0),
            pct=distributed.get("percent_done", 0.0),
            failed=distributed.get("failed_chunks", 0),
            claims=distributed.get("active_claims", 0),
        )
    )
    phases = distributed.get("phase_markers")
    if isinstance(phases, Mapping) and phases:
        lines.append(
            "  phases: "
            + ", ".join(f"{key}={value}" for key, value in sorted(phases.items()))
        )


def append_slurm(lines: list[str], slurm: Any) -> None:
    if not isinstance(slurm, Mapping) or not slurm:
        return
    source = slurm.get("source") or "unknown"
    if not slurm.get("available"):
        lines.append(f"slurm: unavailable ({source})")
        return
    counts = slurm.get("state_counts")
    if isinstance(counts, Mapping) and counts:
        lines.append(
            "slurm: "
            + ", ".join(f"{key}={value}" for key, value in sorted(counts.items()))
        )
    else:
        lines.append(f"slurm: no active/accounting rows ({source})")


def append_eta(lines: list[str], eta: Any) -> None:
    if not isinstance(eta, Mapping):
        return
    eta_sec = eta.get("eta_sec")
    if isinstance(eta_sec, (int, float)):
        lines.append(
            f"eta: {format_seconds(float(eta_sec))} confidence={eta.get('confidence')} basis={eta.get('basis')}"
        )
    else:
        lines.append(f"eta: unknown confidence={eta.get('confidence', 'low')}")


def append_errors(lines: list[str], errors: Any) -> None:
    if not isinstance(errors, Mapping) or not errors:
        return
    lines.append(
        "errors: scanned={scanned} hits={hits} unique={unique}".format(
            scanned=errors.get("scanned_files", 0),
            hits=errors.get("raw_hits", 0),
            unique=errors.get("unique_signatures", 0),
        )
    )
    root_counts = errors.get("root_cause_counts")
    if isinstance(root_counts, Mapping) and root_counts:
        lines.append(
            "  root_causes: "
            + ", ".join(f"{key}={value}" for key, value in list(root_counts.items())[:5])
        )
    for row in (errors.get("signatures") or [])[:3]:
        if not isinstance(row, Mapping):
            continue
        first_line = str(row.get("signature") or "").splitlines()[0]
        lines.append(
            f"  [{row.get('root_cause', 'unknown')}] count={row.get('count', 0)} {first_line[:160]}"
        )


def html_explanations(explanations: Any) -> str:
    if not isinstance(explanations, list) or not explanations:
        return ""
    items = "\n".join(f"<li>{esc(str(item))}</li>" for item in explanations)
    return f"<section><h2>Explain</h2><ul>{items}</ul></section>"


def html_stage_table(stages: Any) -> str:
    if not isinstance(stages, Mapping):
        return ""
    if not any((as_int(mapping(row).get("total")) or 0) > 0 for row in stages.values()):
        return ""
    rows = []
    for stage in STAGE_ORDER:
        row = mapping(stages.get(stage))
        rows.append(
            "<tr>"
            f"<td>{esc(stage)}</td>"
            f"<td>{float(row.get('percent_done') or 0.0):.1f}%</td>"
            f"<td>{int(row.get('completed') or 0)}</td>"
            f"<td>{int(row.get('failed') or 0)}</td>"
            f"<td>{int(row.get('running') or 0)}</td>"
            f"<td>{int(row.get('pending') or 0)}</td>"
            "</tr>"
        )
    return (
        "<section><h2>Stages</h2><table><thead><tr>"
        "<th>Stage</th><th>Done</th><th>Completed</th><th>Failed</th><th>Running</th><th>Pending</th>"
        "</tr></thead><tbody>"
        + "\n".join(rows)
        + "</tbody></table></section>"
    )


def html_slurm(slurm: Mapping[str, Any]) -> str:
    if not slurm:
        return ""
    counts = slurm.get("state_counts")
    if isinstance(counts, Mapping) and counts:
        summary = ", ".join(f"{esc(str(k))}={int(v)}" for k, v in sorted(counts.items()))
    else:
        summary = f"unavailable ({esc(str(slurm.get('source') or 'unknown'))})"
    rows = []
    for job in slurm.get("jobs") or []:
        if not isinstance(job, Mapping):
            continue
        rows.append(
            "<tr>"
            f"<td>{esc(str(job.get('job_id') or ''))}</td>"
            f"<td>{esc(str(job.get('state') or ''))}</td>"
            f"<td>{esc(str(job.get('reason_or_nodelist') or job.get('exit_code') or ''))}</td>"
            "</tr>"
        )
    table = ""
    if rows:
        table = "<table><thead><tr><th>Job</th><th>State</th><th>Reason/Exit</th></tr></thead><tbody>" + "\n".join(rows) + "</tbody></table>"
    return f"<section><h2>Slurm</h2><p>{summary}</p>{table}</section>"


def html_errors(errors: Mapping[str, Any]) -> str:
    if not errors:
        return ""
    root_counts = errors.get("root_cause_counts")
    chips = ""
    if isinstance(root_counts, Mapping):
        chips = " ".join(
            f"<span class=\"pill\">{esc(str(key))}: {int(value)}</span>"
            for key, value in root_counts.items()
        )
    rows = []
    for row in errors.get("signatures") or []:
        if not isinstance(row, Mapping):
            continue
        rows.append(
            "<tr>"
            f"<td>{esc(str(row.get('root_cause') or 'unknown'))}</td>"
            f"<td>{int(row.get('count') or 0)}</td>"
            f"<td><code>{esc(str(row.get('signature') or '').splitlines()[0][:180])}</code></td>"
            "</tr>"
        )
    table = ""
    if rows:
        table = "<table><thead><tr><th>Root cause</th><th>Count</th><th>Signature</th></tr></thead><tbody>" + "\n".join(rows) + "</tbody></table>"
    return f"<section><h2>Errors</h2><p>{chips}</p>{table}</section>"


def html_paths(paths: Any) -> str:
    if not isinstance(paths, Mapping) or not paths:
        return ""
    rows = "\n".join(
        f"<tr><td>{esc(str(key))}</td><td><code>{esc(str(value))}</code></td></tr>"
        for key, value in sorted(paths.items())
    )
    return f"<section><h2>Outputs</h2><table><tbody>{rows}</tbody></table></section>"

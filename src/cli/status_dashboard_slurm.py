from __future__ import annotations

import os
import shutil
import subprocess
from collections import Counter
from typing import Any, Mapping

from cli.status_dashboard_utils import as_float, as_int


def parse_squeue_output(text: str) -> list[dict[str, Any]]:
    jobs = []
    for line in str(text or "").splitlines():
        if not line.strip():
            continue
        parts = line.split("|")
        if len(parts) < 6:
            continue
        jobs.append(
            {
                "job_id": parts[0].strip(),
                "state": parts[1].strip(),
                "elapsed": parts[2].strip(),
                "time_limit": parts[3].strip(),
                "nodes": parts[4].strip(),
                "reason_or_nodelist": parts[5].strip(),
            }
        )
    return jobs


def parse_sacct_output(text: str) -> list[dict[str, Any]]:
    jobs = []
    for line in str(text or "").splitlines():
        if not line.strip():
            continue
        parts = line.split("|")
        if len(parts) < 5:
            continue
        jobs.append(
            {
                "job_id": parts[0].strip(),
                "state": parts[1].strip(),
                "elapsed_sec": as_float(parts[2].strip()),
                "alloc_cpus": as_int(parts[3].strip()),
                "exit_code": parts[4].strip(),
            }
        )
    return jobs


def slurm_job_ids(manifest: Mapping[str, Any] | None) -> list[str]:
    ids: list[str] = []
    for key in ("SLURM_ARRAY_JOB_ID", "SLURM_JOB_ID"):
        raw = str(os.environ.get(key, "") or "").strip()
        if raw and raw not in ids:
            ids.append(raw)
    resources = manifest.get("resources") if isinstance(manifest, Mapping) else None
    slurm = resources.get("slurm") if isinstance(resources, Mapping) else None
    if isinstance(slurm, Mapping):
        for key in ("array_job_id", "job_id"):
            raw = str(slurm.get(key) or "").strip()
            if raw and raw not in ids:
                ids.append(raw)
    return ids


def slurm_status(manifest: Mapping[str, Any] | None) -> dict[str, Any]:
    job_ids = slurm_job_ids(manifest)
    out: dict[str, Any] = {
        "available": bool(shutil.which("squeue") or shutil.which("sacct")),
        "job_ids": job_ids,
        "source": "unavailable",
        "jobs": [],
        "state_counts": {},
    }
    if not job_ids:
        out["source"] = "no_job_id"
        return out
    jobs: list[dict[str, Any]] = []
    if shutil.which("squeue"):
        out["source"] = "squeue"
        for job_id in job_ids:
            jobs.extend(_query_squeue(job_id))
    if not jobs and shutil.which("sacct"):
        out["source"] = "sacct"
        for job_id in job_ids:
            jobs.extend(_query_sacct(job_id))
    counts = Counter(str(job.get("state") or "UNKNOWN") for job in jobs)
    out["jobs"] = jobs
    out["state_counts"] = dict(counts)
    return out


def _query_squeue(job_id: str) -> list[dict[str, Any]]:
    cmd = [
        "squeue",
        "--array",
        "--noheader",
        "--jobs",
        str(job_id),
        "--format=%i|%T|%M|%l|%D|%R",
    ]
    proc = _run_command(cmd)
    if proc is None or proc.returncode != 0:
        return []
    return parse_squeue_output(proc.stdout)


def _query_sacct(job_id: str) -> list[dict[str, Any]]:
    cmd = [
        "sacct",
        "-X",
        "-j",
        str(job_id),
        "--parsable2",
        "--noheader",
        "--format=JobIDRaw,State,ElapsedRaw,AllocCPUS,ExitCode",
    ]
    proc = _run_command(cmd)
    if proc is None or proc.returncode != 0:
        return []
    return parse_sacct_output(proc.stdout)


def _run_command(cmd: list[str]) -> subprocess.CompletedProcess[str] | None:
    try:
        return subprocess.run(cmd, check=False, capture_output=True, text=True, timeout=10)
    except Exception:
        return None

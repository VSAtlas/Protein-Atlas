from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from analysis.ml.external_tracking import log_external_trackers

try:
    from ml.registry import append_registry_index, collect_env_versions, get_git_sha, write_registry_record
except Exception:  # pragma: no cover - optional packaging path in some developer contexts
    append_registry_index = None  # type: ignore[assignment]
    collect_env_versions = None  # type: ignore[assignment]
    get_git_sha = None  # type: ignore[assignment]
    write_registry_record = None  # type: ignore[assignment]


def _append_jsonl(path: Path, row: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(row, sort_keys=True, default=str) + "\n")


def track_model_run(
    *,
    model_dir: str | Path,
    repo_root: str | Path,
    params: dict[str, Any],
    metrics: dict[str, Any],
    dataset_manifest: dict[str, Any],
    artifact_manifest: dict[str, Any],
    claim_readiness: dict[str, Any],
) -> dict[str, Any]:
    root = Path(model_dir)
    repo = Path(repo_root)
    record: dict[str, Any] = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "model_dir": str(root),
        "params": params,
        "metrics": metrics,
        "dataset_version_id": dataset_manifest.get("dataset_version_id"),
        "dataset_path": dataset_manifest.get("dataset_path"),
        "artifact_manifest": str(root / "model_artifacts.json"),
        "artifact_count": len(artifact_manifest.get("artifacts", [])),
        "claim_status": claim_readiness.get("overall_status"),
        "git_sha": get_git_sha(repo) if get_git_sha else "unknown",
        "env": collect_env_versions() if collect_env_versions else {},
    }
    record["external_tracking"] = log_external_trackers(
        model_dir=root,
        params=params,
        metrics=metrics,
        dataset_version_id=str(record["dataset_version_id"]) if record.get("dataset_version_id") else None,
        claim_status=str(record["claim_status"]) if record.get("claim_status") else None,
    )
    _append_jsonl(root / "experiment_runs.jsonl", record)
    if write_registry_record is not None:
        write_registry_record(root, record)
    if append_registry_index is not None:
        append_registry_index(
            root.parent / "model_registry_index.csv",
            {
                "model_dir": str(root),
                "dataset_version_id": record["dataset_version_id"],
                "claim_status": record["claim_status"],
                "selected_model": params.get("selected_model"),
                "split_mode": params.get("split_mode"),
                "AUROC": metrics.get("AUROC"),
                "AUPRC": metrics.get("AUPRC"),
            },
        )
    return record

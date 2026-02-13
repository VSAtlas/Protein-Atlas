from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


def _safe_version(package_name: str) -> str:
    try:
        return importlib.metadata.version(package_name)
    except Exception:
        return "unknown"


def _hash_frame_rows(df: pd.DataFrame, columns: Iterable[str]) -> bytes:
    selected_columns = [col for col in columns if col in df.columns]
    if not selected_columns:
        selected_columns = [col for col in ("lig_smiles", "active", "pocket", "ex_rec_pdb") if col in df.columns]
    if not selected_columns:
        selected_columns = list(df.columns[: min(8, len(df.columns))])
    subset = df.loc[:, selected_columns].copy()
    subset = subset.fillna("")
    row_hashes = pd.util.hash_pandas_object(subset, index=True).values
    return row_hashes.tobytes()


def compute_dataset_hash(
    train_df: pd.DataFrame,
    holdout_df: pd.DataFrame,
    feature_columns: list[str],
    config_snapshot: dict[str, Any],
) -> str:
    digest = hashlib.sha256()
    digest.update(json.dumps(config_snapshot, sort_keys=True, default=str).encode("utf-8"))
    digest.update(json.dumps(sorted(feature_columns), sort_keys=True).encode("utf-8"))

    digest.update(f"train_rows={len(train_df)}".encode("utf-8"))
    digest.update(_hash_frame_rows(train_df, feature_columns))

    digest.update(f"holdout_rows={len(holdout_df)}".encode("utf-8"))
    digest.update(_hash_frame_rows(holdout_df, feature_columns))
    return digest.hexdigest()


def get_git_sha(repo_root: Path) -> str:
    try:
        proc = subprocess.run(
            ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
            check=False,
            capture_output=True,
            text=True,
        )
    except Exception:
        return "unknown"
    if proc.returncode != 0:
        return "unknown"
    value = proc.stdout.strip()
    return value if value else "unknown"


def collect_env_versions() -> dict[str, Any]:
    packages = {
        "numpy": _safe_version("numpy"),
        "pandas": _safe_version("pandas"),
        "scipy": _safe_version("scipy"),
        "scikit-learn": _safe_version("scikit-learn"),
        "rdkit": _safe_version("rdkit"),
        "lightgbm": _safe_version("lightgbm"),
        "xgboost": _safe_version("xgboost"),
        "matplotlib": _safe_version("matplotlib"),
        "optuna": _safe_version("optuna"),
    }
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "packages": packages,
    }


def write_registry_record(run_dir: Path, record: dict[str, Any]) -> Path:
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    if "timestamp_utc" not in record:
        record["timestamp_utc"] = datetime.now(timezone.utc).isoformat()
    out_path = run_dir / "registry.json"
    out_path.write_text(json.dumps(record, indent=2, sort_keys=True, default=str), encoding="utf-8")
    return out_path


def append_registry_index(index_csv: Path, summary_row: dict[str, Any]) -> Path:
    index_path = Path(index_csv)
    index_path.parent.mkdir(parents=True, exist_ok=True)

    row = dict(summary_row)
    if "timestamp_utc" not in row:
        row["timestamp_utc"] = datetime.now(timezone.utc).isoformat()

    existing_header: list[str] | None = None
    if index_path.exists():
        with index_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle)
            existing_header = next(reader, None)

    if existing_header:
        fieldnames = list(existing_header)
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
        existing_rows: list[dict[str, str]] = []
        with index_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for existing in reader:
                existing_rows.append(dict(existing))
        with index_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for existing in existing_rows:
                writer.writerow({name: existing.get(name, "") for name in fieldnames})
            writer.writerow({name: row.get(name, "") for name in fieldnames})
        return index_path

    fieldnames = list(row.keys())
    with index_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(row)
    return index_path

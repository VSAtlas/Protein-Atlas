from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.feature_metadata import refresh_tables
from analysis.ml.labels import binary_label_series
from analysis.ml.audit_utils import source_holdout
from analysis.ml.preflight import discover_default_expert_datasets
from analysis.ml.split_manifest import load_locked_split_manifest, write_split_manifest
from analysis.ml.splits import make_split, split_overlap_summary

DEFAULT_LOCK_SPLITS = [
    "random",
    "drug_holdout",
    "target_holdout",
    "scaffold_holdout",
    "chemical_cluster_holdout",
    "target_family_holdout",
    "temporal_holdout",
    "source_holdout",
]


def lock_run_splits(
    *,
    run_id: str,
    run_dir: str | Path,
    out_dir: str | Path,
    experts: list[str] | None = None,
    split_modes: list[str] | None = None,
    seed: int = 42,
    validation_fraction: float = 0.15,
) -> dict[str, Any]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    selected_experts = set(experts or [])
    modes = split_modes or DEFAULT_LOCK_SPLITS
    results: list[dict[str, Any]] = []
    for item in discover_default_expert_datasets(run_dir):
        expert = str(item["expert"])
        if selected_experts and expert not in selected_experts:
            continue
        dataset = item.get("dataset")
        if not dataset:
            results.append({"expert": expert, "status": "skipped", "reason": "dataset_not_found", "searched": item.get("searched", [])})
            continue
        label_col = str(item.get("label") or "")
        try:
            refresh = refresh_tables(
                [Path(str(dataset))],
                out_dir=out / "_feature_metadata" / expert,
                in_place=False,
                run_dir=Path(run_dir),
                chemical_cluster="auto",
                target_family="auto",
                source_lineage="auto",
                drop_columns=None,
            )
            outputs = refresh.get("outputs") or []
            lock_dataset = Path(outputs[0]) if outputs else Path(str(dataset))
            frame = pd.read_csv(lock_dataset, low_memory=False)
            frame["_atlas_observed_label"] = binary_label_series(frame[label_col])
            data = frame.loc[frame["_atlas_observed_label"].notna()].copy()
            data[label_col] = data["_atlas_observed_label"].astype(int)
        except Exception as exc:
            results.append({"expert": expert, "dataset": dataset, "status": "failed", "error": str(exc)})
            continue
        for split_mode in modes:
            target = out / expert / split_mode
            try:
                if split_mode == "source_holdout":
                    source_col = next((col for col in ["label_source", "source_family", "upstream_source"] if col in data.columns and data[col].notna().any()), None)
                    if source_col is None:
                        raise ValueError("source_holdout requires label_source, source_family, or upstream_source metadata")
                    train_idx, test_idx, summary = source_holdout(data, source_col, seed=seed, test_fraction=0.2)
                else:
                    train_idx, test_idx = make_split(data, split_mode=split_mode, seed=seed)
                    train = data.loc[train_idx]
                    test = data.loc[test_idx]
                    summary = split_overlap_summary(train, test, split_mode)
                train = data.loc[train_idx]
                test = data.loc[test_idx]
                summary["locked_dataset"] = str(lock_dataset)
                summary["source_dataset"] = str(dataset)
                if summary.get("overlaps") and not summary.get("passes_holdout"):
                    raise ValueError(f"holdout split leakage detected: {summary}")
                validation_idx = pd.Index([])
                if validation_fraction > 0 and len(train) >= 5:
                    validation_idx = train.sample(frac=validation_fraction, random_state=seed).index
                    train_idx = train.index.difference(validation_idx)
                write_split_manifest(
                    data,
                    train_idx,
                    test_idx,
                    target,
                    label_col=label_col,
                    split_mode=split_mode,
                    split_summary=summary,
                    validation_idx=validation_idx,
                )
                results.append(
                    {
                        "expert": expert,
                        "dataset": str(lock_dataset),
                        "source_dataset": dataset,
                        "label": label_col,
                        "split_mode": split_mode,
                        "status": "locked",
                        "split_manifest": str(target / "split_manifest.csv"),
                        "split_metadata": str(target / "split_manifest.json"),
                    }
                )
            except Exception as exc:
                results.append({"expert": expert, "dataset": dataset, "label": label_col, "split_mode": split_mode, "status": "skipped", "reason": str(exc)})
    manifest = {"status": "written", "run_id": run_id, "run_dir": str(run_dir), "out_dir": str(out), "seed": seed, "validation_fraction": validation_fraction, "results": results}
    (out / "split_lock_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True, default=str), encoding="utf-8")
    return manifest


def validate_run_splits(*, run_dir: str | Path, split_root: str | Path) -> dict[str, Any]:
    root = Path(split_root)
    results: list[dict[str, Any]] = []
    dataset_by_expert = {item["expert"]: item for item in discover_default_expert_datasets(run_dir)}
    for metadata_path in sorted(root.glob("*/*/split_manifest.json")):
        expert = metadata_path.parent.parent.name
        item = dataset_by_expert.get(expert)
        if not item or not item.get("dataset"):
            results.append({"expert": expert, "split_mode": metadata_path.parent.name, "status": "failed", "error": "expert dataset not found"})
            continue
        try:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            locked_dataset = Path(str((metadata.get("split_summary") or {}).get("locked_dataset", "")))
            if locked_dataset.exists():
                validate_dataset = locked_dataset
            else:
                refresh = refresh_tables(
                    [Path(str(item["dataset"]))],
                    out_dir=root / "_feature_metadata_validate" / expert,
                    in_place=False,
                    run_dir=Path(run_dir),
                    chemical_cluster="auto",
                    target_family="auto",
                    source_lineage="auto",
                    drop_columns=None,
                )
                outputs = refresh.get("outputs") or []
                validate_dataset = Path(outputs[0]) if outputs else Path(str(item["dataset"]))
            frame = pd.read_csv(validate_dataset, low_memory=False)
            labels = binary_label_series(frame[str(item.get("label"))])
            frame["_atlas_observed_label"] = labels
            data = frame.loc[labels.notna()].copy()
            data[str(item.get("label"))] = labels.loc[labels.notna()].astype(int)
            loaded = load_locked_split_manifest(data, metadata_path.parent, label_col=str(item.get("label")), split_mode=metadata_path.parent.name)
            results.append({"expert": expert, "split_mode": metadata_path.parent.name, "status": "valid", "n_train": len(loaded["train_idx"]), "n_test": len(loaded["test_idx"])})
        except Exception as exc:
            results.append({"expert": expert, "split_mode": metadata_path.parent.name, "status": "failed", "error": str(exc)})
    status = "valid" if results and all(row.get("status") == "valid" for row in results) else "failed"
    return {"split_root": str(root), "status": status, "results": results}

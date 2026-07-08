from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.labels import binary_label_series
from analysis.ml.split_manifest import dataframe_content_hash, load_locked_split_manifest, write_split_manifest


def _allocate_counts(groups: pd.Series, n_rows: int) -> dict[str, int]:
    counts = groups.value_counts(dropna=False).to_dict()
    total = int(sum(counts.values()))
    if total <= n_rows:
        return {str(key): int(value) for key, value in counts.items()}
    allocated: dict[str, int] = {}
    for key, value in counts.items():
        share = max(1, int(round((int(value) / total) * n_rows)))
        allocated[str(key)] = min(int(value), share)
    while sum(allocated.values()) > n_rows:
        key = max(allocated, key=lambda item: allocated[item])
        if allocated[key] > 1:
            allocated[key] -= 1
        else:
            break
    while sum(allocated.values()) < n_rows:
        for key, value in counts.items():
            skey = str(key)
            if allocated.get(skey, 0) < int(value):
                allocated[skey] = allocated.get(skey, 0) + 1
                if sum(allocated.values()) >= n_rows:
                    break
        else:
            break
    return allocated


def materialize_sampled_dataset(
    dataset: str | Path,
    *,
    label_col: str,
    out_dir: str | Path,
    sample_rows: int | None,
    seed: int = 42,
    strategy: str = "stratified",
) -> dict[str, Any]:
    """Write a deterministic sampled copy of a model-ready table.

    The sample is intentionally file-backed so every downstream artifact records a
    concrete dataset path and hash, just like full-data training.
    """

    source = Path(dataset)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(source, low_memory=False)
    source_hash = dataframe_content_hash(frame)
    requested = int(sample_rows or 0)
    if requested <= 0 or len(frame) <= requested:
        sample = frame.copy()
        sampled = False
    elif strategy == "stratified" and label_col in frame.columns:
        labels = binary_label_series(frame[label_col]).fillna("unlabeled").astype(str)
        allocations = _allocate_counts(labels, requested)
        chunks = []
        for value, group in frame.groupby(labels, dropna=False):
            take = allocations.get(str(value), 0)
            if take <= 0:
                continue
            chunks.append(group.sample(n=min(take, len(group)), random_state=seed))
        sample = pd.concat(chunks, axis=0).sample(frac=1.0, random_state=seed).reset_index(drop=True) if chunks else frame.head(0)
        sampled = True
    else:
        sample = frame.sample(n=requested, random_state=seed).reset_index(drop=True)
        sampled = True
    sample_path = out / "sampled_dataset.csv"
    sample.to_csv(sample_path, index=False)
    manifest = {
        "source_dataset": str(source),
        "sampled_dataset": str(sample_path),
        "source_dataset_hash": source_hash,
        "sample_dataset_hash": dataframe_content_hash(sample),
        "source_rows": int(len(frame)),
        "sample_rows": int(len(sample)),
        "requested_sample_rows": requested,
        "sampled": sampled,
        "label_col": label_col,
        "strategy": strategy,
        "seed": int(seed),
    }
    if label_col in sample.columns:
        labels = binary_label_series(sample[label_col])
        manifest["sample_label_counts"] = {str(k): int(v) for k, v in labels.value_counts(dropna=False).items()}
    (out / "sample_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


def materialize_split_locked_sampled_dataset(
    dataset: str | Path,
    *,
    label_col: str,
    split_manifest: str | Path,
    out_dir: str | Path,
    sample_rows: int | None,
    seed: int = 42,
    strategy: str = "stratified",
    split_mode: str | None = None,
) -> dict[str, Any]:
    """Write a deterministic sampled table plus a matching locked split manifest.

    Full-run split manifests intentionally hash the full metadata-refreshed table.
    When a user asks for a fast sample, the original manifest cannot be reused
    directly because row identities and dataset hashes change. This helper samples
    within the existing train/validation/test assignment and writes a derived
    manifest for the sampled copy.
    """

    source = Path(dataset)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(source, low_memory=False)
    frame["_atlas_observed_label"] = binary_label_series(frame[label_col])
    data = frame.loc[frame["_atlas_observed_label"].notna()].copy()
    data[label_col] = data["_atlas_observed_label"].astype(int)
    locked = load_locked_split_manifest(
        data,
        split_manifest,
        label_col=label_col,
        split_mode=split_mode,
    )
    split_by_index = pd.Series("unused", index=data.index, dtype=object)
    split_by_index.loc[locked["train_idx"]] = "train"  # type: ignore[index]
    split_by_index.loc[locked["test_idx"]] = "test"  # type: ignore[index]
    split_by_index.loc[locked["validation_idx"]] = "validation"  # type: ignore[index]

    requested = int(sample_rows or 0)
    source_hash = dataframe_content_hash(data)
    if requested <= 0 or len(data) <= requested:
        selected = data.index
        sampled = False
    else:
        group_frame = pd.DataFrame(
            {
                "split": split_by_index,
                "label": data[label_col].astype(str),
            },
            index=data.index,
        )
        if strategy == "stratified":
            groups = group_frame["split"].astype(str) + "|" + group_frame["label"].astype(str)
        else:
            groups = group_frame["split"].astype(str)
        allocations = _allocate_counts(groups, requested)
        chunks = []
        for value, group in data.groupby(groups, dropna=False):
            take = allocations.get(str(value), 0)
            if take <= 0:
                continue
            chunks.append(group.sample(n=min(take, len(group)), random_state=seed))
        if not chunks:
            raise ValueError("sampled locked split produced no rows")
        selected = pd.concat(chunks, axis=0).sample(frac=1.0, random_state=seed).index
        sampled = True

    sample = data.loc[selected].copy().reset_index(drop=True)
    selected_split = split_by_index.loc[selected].reset_index(drop=True)
    train_idx = sample.index[selected_split == "train"]
    validation_idx = sample.index[selected_split == "validation"]
    test_idx = sample.index[selected_split == "test"]
    if len(train_idx) == 0 or len(test_idx) == 0:
        raise ValueError(
            "sampled locked split produced empty train/test: "
            f"train={len(train_idx)} test={len(test_idx)}; increase --sample-rows"
        )
    sample_path = out / "sampled_dataset.csv"
    sample.to_csv(sample_path, index=False)
    sampled_split_dir = out / "split_manifest"
    split_summary = dict(locked.get("split_summary") or {})
    split_summary["sampled_from_locked_split"] = True
    split_summary["source_locked_split_manifest"] = str(split_manifest)
    write_split_manifest(
        sample,
        train_idx,
        test_idx,
        sampled_split_dir,
        label_col=label_col,
        split_mode=str(split_mode or split_summary.get("split_mode") or "locked_sample"),
        split_summary=split_summary,
        validation_idx=validation_idx,
    )
    manifest = {
        "source_dataset": str(source),
        "sampled_dataset": str(sample_path),
        "source_dataset_hash": source_hash,
        "sample_dataset_hash": dataframe_content_hash(sample),
        "source_rows": int(len(data)),
        "sample_rows": int(len(sample)),
        "requested_sample_rows": requested,
        "sampled": sampled,
        "label_col": label_col,
        "strategy": strategy,
        "seed": int(seed),
        "split_locked": True,
        "source_split_manifest": str(split_manifest),
        "sampled_split_manifest": str(sampled_split_dir),
        "split_counts": {str(k): int(v) for k, v in selected_split.value_counts(dropna=False).items()},
    }
    labels = binary_label_series(sample[label_col])
    manifest["sample_label_counts"] = {str(k): int(v) for k, v in labels.value_counts(dropna=False).items()}
    (out / "sample_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


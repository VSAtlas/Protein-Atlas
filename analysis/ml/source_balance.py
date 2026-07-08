from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


def build_source_balanced_dataset(
    dataset_path: str | Path,
    out_path: str | Path,
    *,
    label_col: str,
    source_col: str = "negative_source",
    negative_label: int = 0,
    seed: int = 42,
    max_per_source: int | None = None,
    min_per_source: int = 1,
) -> pd.DataFrame:
    df = pd.read_csv(dataset_path, low_memory=False)
    labels = pd.to_numeric(df[label_col], errors="coerce")
    labeled = df[labels.notna()].copy()
    labeled["_balance_label"] = labels.loc[labeled.index].astype(int)
    non_negative = labeled[labeled["_balance_label"].ne(int(negative_label))].copy()
    negatives = labeled[labeled["_balance_label"].eq(int(negative_label))].copy()
    if negatives.empty or source_col not in negatives.columns:
        out = labeled.drop(columns=["_balance_label"], errors="ignore")
    else:
        source_counts = negatives[source_col].fillna("missing").astype(str).value_counts()
        cap = int(max_per_source) if max_per_source is not None else int(source_counts.min())
        cap = max(int(min_per_source), cap)
        sampled = []
        for source, group in negatives.groupby(negatives[source_col].fillna("missing").astype(str), sort=True):
            n = min(len(group), cap)
            sampled.append(group.sample(n=n, random_state=seed) if n < len(group) else group)
        balanced_negatives = pd.concat(sampled, ignore_index=False) if sampled else negatives.iloc[0:0]
        out = pd.concat([non_negative, balanced_negatives], ignore_index=False)
        out = out.sample(frac=1.0, random_state=seed).drop(columns=["_balance_label"], errors="ignore")
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    manifest = {
        "input": str(dataset_path),
        "output": str(path),
        "label_col": label_col,
        "source_col": source_col,
        "negative_label": negative_label,
        "seed": seed,
        "max_per_source": max_per_source,
        "rows_in": int(len(df)),
        "rows_out": int(len(out)),
        "label_counts_out": {
            str(k): int(v)
            for k, v in pd.to_numeric(out[label_col], errors="coerce").value_counts(dropna=False).items()
        },
        "source_counts_out": {
            str(k): int(v) for k, v in out.get(source_col, pd.Series(dtype=object)).value_counts(dropna=False).items()
        },
    }
    path.with_suffix(".manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return out


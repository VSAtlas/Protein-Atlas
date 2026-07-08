from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd


MANIFEST_ID_COLS = [
    "canonical_pair_key",
    "drug_id",
    "target_id",
    "pdb_id",
    "scaffold_key",
    "chemical_cluster",
    "target_family",
    "protein_class",
    "label_source",
    "source_family",
    "upstream_source",
    "assay_type",
    "endpoint_type",
    "activity_type",
]


def dataframe_content_hash(df: pd.DataFrame) -> str:
    """Return a stable hash for the current row/column content."""

    cols = sorted(df.columns)
    normalized = df.reindex(columns=cols).astype("object").where(pd.notna(df.reindex(columns=cols)), "<NA>").astype(str)
    digest = hashlib.sha256()
    digest.update("|".join(cols).encode("utf-8"))
    row_hashes = pd.util.hash_pandas_object(normalized, index=True).astype(str)
    digest.update("\n".join(row_hashes.tolist()).encode("utf-8"))
    return digest.hexdigest()


def row_identity_hash(row: pd.Series, cols: list[str]) -> str:
    values = [str(row.get(col, "") if pd.notna(row.get(col, "")) else "") for col in cols]
    return hashlib.sha256("|".join(values).encode("utf-8")).hexdigest()


def _split_manifest_hash(table: pd.DataFrame) -> str:
    digest = hashlib.sha256()
    cols = sorted(table.columns)
    normalized = table.reindex(columns=cols).astype("object").where(pd.notna(table.reindex(columns=cols)), "<NA>").astype(str)
    digest.update("|".join(cols).encode("utf-8"))
    row_hashes = pd.util.hash_pandas_object(normalized, index=False).astype(str)
    digest.update("\n".join(row_hashes.tolist()).encode("utf-8"))
    return digest.hexdigest()


def write_split_manifest(
    data: pd.DataFrame,
    train_idx: pd.Index,
    test_idx: pd.Index,
    out_dir: str | Path,
    *,
    label_col: str,
    split_mode: str,
    split_summary: dict[str, Any],
    validation_idx: pd.Index | None = None,
) -> pd.DataFrame:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ids = [col for col in MANIFEST_ID_COLS if col in data.columns]
    table = data[ids].copy() if ids else pd.DataFrame(index=data.index)
    table["_row_index"] = data.index.astype(str)
    table["split"] = "unused"
    table.loc[train_idx, "split"] = "train"
    table.loc[test_idx, "split"] = "test"
    if validation_idx is not None and len(validation_idx):
        table.loc[validation_idx, "split"] = "validation"
    table["fold"] = split_mode
    table["label_col"] = label_col
    if label_col in data.columns:
        table[label_col] = data[label_col]
    hash_cols = [col for col in ["canonical_pair_key", "drug_id", "target_id", "pdb_id", label_col] if col in table.columns]
    if not hash_cols:
        hash_cols = ["_row_index"]
    table["row_identity_hash"] = table.apply(lambda row: row_identity_hash(row, hash_cols), axis=1)
    table.to_csv(out / "split_manifest.csv", index=False)
    manifest_hash = _split_manifest_hash(table)
    metadata = {
        "dataset_hash": dataframe_content_hash(data),
        "split_manifest_hash": manifest_hash,
        "label_col": label_col,
        "split_mode": split_mode,
        "n_rows": int(len(data)),
        "n_train": int((table["split"] == "train").sum()),
        "n_validation": int((table["split"] == "validation").sum()),
        "n_test": int((table["split"] == "test").sum()),
        "split_summary": split_summary,
        "manifest": str(out / "split_manifest.csv"),
    }
    (out / "split_manifest.json").write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")
    return table



def _manifest_paths(path: str | Path) -> tuple[Path, Path]:
    root = Path(path)
    if root.is_dir():
        return root / "split_manifest.csv", root / "split_manifest.json"
    if root.name.endswith(".csv"):
        return root, root.with_suffix(".json")
    if root.name.endswith(".json"):
        return root.with_name("split_manifest.csv"), root
    return root / "split_manifest.csv", root / "split_manifest.json"


def load_locked_split_manifest(
    data: pd.DataFrame,
    manifest_path: str | Path,
    *,
    label_col: str,
    split_mode: str | None = None,
) -> dict[str, object]:
    csv_path, json_path = _manifest_paths(manifest_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"split manifest CSV not found: {csv_path}")
    if not json_path.exists():
        raise FileNotFoundError(f"split manifest JSON not found: {json_path}")
    table = pd.read_csv(csv_path, dtype=str).fillna("")
    metadata = json.loads(json_path.read_text(encoding="utf-8"))
    dataset_hash = dataframe_content_hash(data)
    if metadata.get("dataset_hash") != dataset_hash:
        raise ValueError(
            "split manifest dataset hash mismatch: "
            f"manifest={metadata.get('dataset_hash')} current={dataset_hash}"
        )
    if label_col and str(metadata.get("label_col")) != str(label_col):
        raise ValueError(f"split manifest label mismatch: manifest={metadata.get('label_col')} current={label_col}")
    if split_mode and str(metadata.get("split_mode")) != str(split_mode):
        raise ValueError(
            "split manifest mode mismatch: "
            f"manifest={metadata.get('split_mode')} requested={split_mode}"
        )
    if split_mode and "fold" in table.columns:
        fold_modes = {str(value) for value in table["fold"].dropna().unique() if str(value)}
        if fold_modes and fold_modes != {str(split_mode)}:
            raise ValueError(
                "split manifest fold mismatch: "
                f"manifest={sorted(fold_modes)} requested={split_mode}"
            )
    expected_rows = {str(idx) for idx in data.index.astype(str)}
    manifest_rows = set(table.get("_row_index", pd.Series(dtype=str)).astype(str))
    if expected_rows != manifest_rows:
        missing = len(expected_rows - manifest_rows)
        extra = len(manifest_rows - expected_rows)
        raise ValueError(f"split manifest row identity mismatch: missing={missing} extra={extra}")
    if metadata.get("split_manifest_hash") and _split_manifest_hash(table) != metadata.get("split_manifest_hash"):
        raise ValueError("split manifest hash mismatch; manifest table may have been edited")
    split_by_row = table.set_index("_row_index")["split"].astype(str)
    row_index = data.index.astype(str)
    train_idx = data.index[row_index.map(split_by_row).isin(["train"])]
    validation_idx = data.index[row_index.map(split_by_row).isin(["validation"])]
    test_idx = data.index[row_index.map(split_by_row).isin(["test"])]
    if len(train_idx) == 0 or len(test_idx) == 0:
        raise ValueError(f"locked split produced empty train/test: train={len(train_idx)} test={len(test_idx)}")
    return {
        "metadata": metadata,
        "table": table,
        "train_idx": train_idx,
        "validation_idx": validation_idx,
        "test_idx": test_idx,
        "split_summary": metadata.get("split_summary") or {},
    }

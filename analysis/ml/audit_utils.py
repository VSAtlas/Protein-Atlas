from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd


SOURCE_COLS = [
    "source_name",
    "source",
    "source_objective",
    "evidence_sources",
    "parent_sources",
    "label_source",
    "dataset",
    "source_family",
    "upstream_source",
]

OVERLAP_FIELDS = [
    "drug_id",
    "target_id",
    "pdb_id",
    "ligand_chemotype",
    "scaffold_key",
    "chemical_cluster",
    "ligand_cluster",
    "target_family",
    "protein_family",
    "protein_class",
    "canonical_pair_key",
    "dedup_drug_key",
    "dedup_target_key",
]

HIGH_RISK_COLUMN_TOKENS = [
    "label",
    "threshold",
    "assay_id",
    "assay_count",
    "activity",
    "ac50",
    "cmax",
    "exposure_margin",
    "source",
    "status",
    "available",
    "rescored",
]

LINEAGE_COLS = [
    "source_family",
    "upstream_source",
    "parent_sources",
    "source_release",
    "source_version",
    "label_source",
]

ASSAY_CONTEXT_COLS = [
    "assay_type",
    "assay_mode",
    "endpoint_type",
    "activity_type",
    "standard_type",
    "relation_domain",
    "source_objective",
]


def clean_value(value: Any) -> str:
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().lower().split())


def clean_values(series: pd.Series) -> set[str]:
    return {
        clean_value(value)
        for value in series.dropna().tolist()
        if clean_value(value) and clean_value(value) != "nan"
    }


def split_source_tokens(value: Any) -> list[str]:
    text = clean_value(value).replace("|", ";").replace(",", ";")
    return [token.strip() for token in text.split(";") if token.strip()]


def first_col(df: pd.DataFrame, names: list[str]) -> str | None:
    for name in names:
        if name in df.columns:
            return name
    return None


def load_table(path: str | Path) -> pd.DataFrame:
    table = Path(path)
    suffix = table.suffix.lower()
    if suffix in {".tsv", ".tab"}:
        return pd.read_csv(table, sep="\t", low_memory=False)
    if suffix == ".parquet":
        return pd.read_parquet(table)
    return pd.read_csv(table, low_memory=False)


def parse_named_path(value: str) -> tuple[str, Path]:
    if "=" in value:
        name, path_text = value.split("=", 1)
        return name.strip(), Path(path_text.strip())
    parsed_path = Path(value)
    return parsed_path.stem, parsed_path


def ensure_canonical_pair_key(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "canonical_pair_key" not in out.columns and {"drug_id", "target_id"}.issubset(out.columns):
        drug = out["drug_id"].fillna("").astype(str).str.strip().str.lower()
        target = out["target_id"].fillna("").astype(str).str.strip().str.upper()
        out["canonical_pair_key"] = "drug_target|" + drug + "|" + target
    return out


def canonical_keys(df: pd.DataFrame) -> pd.Series:
    if "canonical_pair_key" in df.columns:
        return df["canonical_pair_key"].map(clean_value)
    if {"drug_id", "target_id"}.issubset(df.columns):
        return "drug_target|" + df["drug_id"].map(clean_value) + "|" + df["target_id"].map(clean_value)
    adr_col = first_col(df, ["adr_id", "adr_term", "meddra_term", "outcomeName"])
    if adr_col is not None and "drug_id" in df.columns:
        return "drug_adr|" + df["drug_id"].map(clean_value) + "|" + df[adr_col].map(clean_value)
    return pd.Series("", index=df.index, dtype="object")


def source_tokens(df: pd.DataFrame) -> set[str]:
    tokens: set[str] = set()
    for col in SOURCE_COLS:
        if col not in df.columns:
            continue
        for value in df[col].dropna().astype(str):
            tokens.update(split_source_tokens(value))
    return tokens


def source_lineage_tokens(df: pd.DataFrame) -> set[str]:
    tokens: set[str] = set()
    for col in LINEAGE_COLS:
        if col not in df.columns:
            continue
        for value in df[col].dropna().astype(str):
            tokens.update(split_source_tokens(value))
    return tokens


def id_set(df: pd.DataFrame, col: str) -> set[str]:
    if col not in df.columns:
        return set()
    return {value for value in df[col].map(clean_value).tolist() if value}


def label_balance(df: pd.DataFrame, label_col: str) -> dict[str, Any]:
    if label_col not in df.columns:
        return {"available": False}
    labels = pd.to_numeric(df[label_col], errors="coerce")
    observed = labels.dropna()
    return {
        "available": True,
        "n_labeled": int(observed.shape[0]),
        "positive": int((observed == 1).sum()),
        "negative": int((observed == 0).sum()),
        "excluded_or_ambiguous": int((observed == -1).sum()),
        "unknown_or_missing": int(labels.isna().sum()),
        "positive_rate": float((observed == 1).mean()) if len(observed) else None,
    }


def overlap_row(train: pd.DataFrame, test: pd.DataFrame, field: str, split_mode: str) -> dict[str, Any]:
    if field not in train.columns or field not in test.columns:
        return {
            "split_mode": split_mode,
            "field": field,
            "available": False,
            "train_unique": 0,
            "test_unique": 0,
            "overlap_unique": 0,
            "test_overlap_fraction": None,
        }
    train_values = clean_values(train[field])
    test_values = clean_values(test[field])
    overlap = train_values & test_values
    return {
        "split_mode": split_mode,
        "field": field,
        "available": True,
        "train_unique": len(train_values),
        "test_unique": len(test_values),
        "overlap_unique": len(overlap),
        "test_overlap_fraction": len(overlap) / len(test_values) if test_values else 0.0,
    }


def duplicate_label_conflicts(data: pd.DataFrame, label_col: str) -> dict[str, Any]:
    if "canonical_pair_key" not in data.columns or label_col not in data.columns:
        return {"available": False}
    observed = data.dropna(subset=["canonical_pair_key", label_col]).copy()
    labels = pd.to_numeric(observed[label_col], errors="coerce")
    observed = observed.loc[labels.notna()].copy()
    observed["_audit_label"] = labels.loc[labels.notna()].astype(int)
    counts = observed.groupby("canonical_pair_key")["_audit_label"].nunique()
    conflicts = counts[counts > 1]
    return {
        "available": True,
        "duplicate_pair_rows": int(observed.duplicated("canonical_pair_key", keep=False).sum()),
        "conflicting_pair_keys": int(conflicts.shape[0]),
    }


def source_holdout(
    data: pd.DataFrame,
    source_col: str,
    seed: int,
    test_fraction: float,
) -> tuple[pd.Index, pd.Index, dict[str, Any]]:
    if source_col not in data.columns:
        raise ValueError(f"source_holdout requires {source_col}")
    exploded = data[source_col].map(split_source_tokens)
    sources = sorted({source for tokens in exploded for source in tokens})
    if len(sources) < 2:
        raise ValueError(f"source_holdout requires at least two distinct {source_col} values")
    held_out = set(pd.Series(sources).sample(frac=1.0, random_state=seed).tolist()[: max(1, int(round(len(sources) * test_fraction)))])
    test_mask = exploded.map(lambda tokens: bool(set(tokens) & held_out))
    train_mask = ~test_mask
    return data.index[train_mask], data.index[test_mask], {
        "split_mode": "source_holdout",
        "source_col": source_col,
        "held_out_sources": sorted(held_out),
        "n_train": int(train_mask.sum()),
        "n_test": int(test_mask.sum()),
        "passes_holdout": True,
        "overlaps": {source_col: 0},
    }


def feature_missingness(df: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    rows = []
    for feature in features:
        if feature not in df.columns:
            rows.append({"feature": feature, "available": False, "missing_fraction": None, "n_unique": 0})
            continue
        values = df[feature]
        rows.append(
            {
                "feature": feature,
                "available": True,
                "missing_fraction": float(values.isna().mean()),
                "n_unique": int(values.nunique(dropna=True)),
                "dtype": str(values.dtype),
            }
        )
    return pd.DataFrame(rows)


def source_label_balance(df: pd.DataFrame, label_col: str, source_col: str) -> pd.DataFrame:
    if source_col not in df.columns or label_col not in df.columns:
        return pd.DataFrame(columns=[source_col, "n", "positive", "negative", "positive_rate"])
    rows: list[dict[str, Any]] = []
    work = df.copy()
    work["_audit_sources"] = work[source_col].map(split_source_tokens)
    labels = pd.to_numeric(work[label_col], errors="coerce")
    work = work.loc[labels.notna()].copy()
    work[label_col] = labels.loc[labels.notna()].astype(int)
    for source, group in work.explode("_audit_sources").groupby("_audit_sources", dropna=False):
        if not source:
            continue
        rows.append(
            {
                source_col: source,
                "n": int(len(group)),
                "positive": int((group[label_col] == 1).sum()),
                "negative": int((group[label_col] == 0).sum()),
                "positive_rate": float((group[label_col] == 1).mean()) if len(group) else None,
            }
        )
    return pd.DataFrame(rows).sort_values(["n", source_col], ascending=[False, True]) if rows else pd.DataFrame(rows)


def grouped_label_balance(df: pd.DataFrame, label_col: str, group_cols: list[str]) -> pd.DataFrame:
    labels = pd.to_numeric(df[label_col], errors="coerce") if label_col in df.columns else pd.Series(dtype=float)
    data = df.loc[labels.notna()].copy()
    if data.empty:
        return pd.DataFrame(columns=[*group_cols, "n", "positive", "negative", "positive_rate"])
    data[label_col] = labels.loc[labels.notna()].astype(int)
    available = [col for col in group_cols if col in data.columns]
    if not available:
        return pd.DataFrame(columns=[*group_cols, "n", "positive", "negative", "positive_rate"])
    rows: list[dict[str, Any]] = []
    for key, group in data.groupby(available, dropna=False):
        if not isinstance(key, tuple):
            key = (key,)
        row = {col: value for col, value in zip(available, key, strict=False)}
        row.update(
            {
                "n": int(len(group)),
                "positive": int((group[label_col] == 1).sum()),
                "negative": int((group[label_col] == 0).sum()),
                "positive_rate": float((group[label_col] == 1).mean()) if len(group) else None,
            }
        )
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["n", *available], ascending=[False, *([True] * len(available))])


def label_evidence_composition(df: pd.DataFrame, label_col: str) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    labels = pd.to_numeric(df[label_col], errors="coerce") if label_col in df.columns else pd.Series(pd.NA, index=df.index)
    benchmark_only = (
        df["benchmark_only"].fillna(False).astype(str).str.lower().isin({"1", "true", "yes"})
        if "benchmark_only" in df.columns
        else pd.Series(False, index=df.index)
    )
    sample_weight = pd.to_numeric(df["_sample_weight"], errors="coerce") if "_sample_weight" in df.columns else pd.Series(pd.NA, index=df.index)
    weak_negative = (
        df["negative_evidence_type"].fillna("").astype(str).str.contains("weak|faers|nonsignal|reliable", case=False, regex=True)
        if "negative_evidence_type" in df.columns
        else pd.Series(False, index=df.index)
    )
    rows.extend(
        [
            {"category": "strict_positive", "rows": int(labels.eq(1).sum())},
            {"category": "strict_negative", "rows": int((labels.eq(0) & ~weak_negative & ~benchmark_only).sum())},
            {"category": "weak_negative", "rows": int((labels.eq(0) & weak_negative & ~benchmark_only).sum())},
            {"category": "PU_background", "rows": int((labels.isna() & sample_weight.notna()).sum())},
            {"category": "unknown_unlabeled", "rows": int(labels.isna().sum())},
            {"category": "excluded_or_ambiguous", "rows": int(labels.eq(-1).sum())},
            {"category": "benchmark_only_excluded", "rows": int(benchmark_only.sum())},
        ]
    )
    total = max(1, len(df))
    for row in rows:
        row["fraction"] = row["rows"] / total
    return pd.DataFrame(rows)


def high_risk_column_scan(columns: list[str], label_col: str, features: list[str]) -> pd.DataFrame:
    feature_set = set(features)
    rows = []
    for column in columns:
        lower = column.lower()
        reasons = [token for token in HIGH_RISK_COLUMN_TOKENS if token in lower]
        if not reasons:
            continue
        rows.append(
            {
                "column": column,
                "is_selected_feature": column in feature_set,
                "is_label": column == label_col,
                "risk_tokens": ";".join(reasons),
            }
        )
    return pd.DataFrame(rows)

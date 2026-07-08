from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.audit_utils import (
    SOURCE_COLS,
    canonical_keys,
    clean_value,
    id_set,
    load_table,
    parse_named_path,
    source_tokens,
    source_lineage_tokens,
)


def audit_dataset_independence(
    train_path: str | Path,
    candidates: list[str],
    out_dir: str | Path,
    *,
    strict_pair_overlap: bool = True,
    strict_source_overlap: bool = False,
) -> dict[str, Any]:
    """Audit whether benchmark/calibration tables are independent of training rows.

    The hard default is zero canonical pair overlap. Source overlap is reported but
    not a default failure, because within-source sensitivity analyses are sometimes
    intentional and should be named as such rather than silently blocked.
    """

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    train = load_table(train_path)
    train_keys = set(canonical_keys(train).dropna())
    train_keys.discard("")
    train_sources = source_tokens(train)
    train_lineage = source_lineage_tokens(train)
    train_drugs = id_set(train, "drug_id")
    train_targets = id_set(train, "target_id")

    rows: list[dict[str, Any]] = []
    for item in candidates:
        name, path = parse_named_path(item)
        data = load_table(path)
        keys = set(canonical_keys(data).dropna())
        keys.discard("")
        sources = source_tokens(data)
        lineage = source_lineage_tokens(data)
        drugs = id_set(data, "drug_id")
        targets = id_set(data, "target_id")
        pair_overlap = train_keys & keys
        source_overlap = train_sources & sources
        lineage_overlap = train_lineage & lineage
        drug_overlap = train_drugs & drugs
        target_overlap = train_targets & targets
        fails_pair = strict_pair_overlap and bool(pair_overlap)
        fails_source = strict_source_overlap and bool(source_overlap | lineage_overlap)
        rows.append(
            {
                "candidate": name,
                "path": str(path),
                "rows": int(len(data)),
                "candidate_pairs": int(len(keys)),
                "train_pair_overlap": int(len(pair_overlap)),
                "train_pair_overlap_fraction": len(pair_overlap) / len(keys) if keys else 0.0,
                "source_overlap_count": int(len(source_overlap)),
                "source_overlap_values": ";".join(sorted(source_overlap)),
                "source_lineage_overlap_count": int(len(lineage_overlap)),
                "source_lineage_overlap_values": ";".join(sorted(lineage_overlap)),
                "drug_overlap_count": int(len(drug_overlap)),
                "target_overlap_count": int(len(target_overlap)),
                "passes_strict_independence": not (fails_pair or fails_source),
                "interpretation": _interpret(bool(pair_overlap), bool(source_overlap), strict_source_overlap),
            }
        )
    table = pd.DataFrame(rows)
    table_path = out / "dataset_independence_audit.csv"
    table.to_csv(table_path, index=False)
    manifest = {
        "train_path": str(train_path),
        "train_rows": int(len(train)),
        "train_pairs": int(len(train_keys)),
        "strict_pair_overlap": strict_pair_overlap,
        "strict_source_overlap": strict_source_overlap,
        "all_candidates_passed": bool(table["passes_strict_independence"].all()) if not table.empty else True,
        "outputs": {"table": str(table_path)},
        "warning": (
            "Zero canonical-pair overlap is required for pure external benchmark/calibration claims. "
            "Source overlap is reported separately because within-source sensitivity analyses may be intentional."
        ),
    }
    manifest_path = out / "dataset_independence_audit.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


def build_independent_training_table(
    train_path: str | Path,
    heldout: list[str],
    out_path: str | Path,
    *,
    exclude_source_tokens: list[str] | None = None,
) -> pd.DataFrame:
    train = load_table(train_path)
    train_keys = canonical_keys(train)
    heldout_keys: set[str] = set()
    for item in heldout:
        _, path = parse_named_path(item)
        data = load_table(path)
        keys = set(canonical_keys(data).dropna())
        keys.discard("")
        heldout_keys.update(keys)
    keep = ~train_keys.isin(heldout_keys)
    tokens = {clean_value(token) for token in (exclude_source_tokens or []) if clean_value(token)}
    if tokens:
        source_text = pd.Series("", index=train.index, dtype="object")
        for col in SOURCE_COLS:
            if col in train.columns:
                source_text = source_text + ";" + train[col].fillna("").astype(str).str.lower()
        for token in tokens:
            keep &= ~source_text.str.contains(token, regex=False, na=False)
    out = train.loc[keep].copy()
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    manifest = {
        "input_train_path": str(train_path),
        "heldout": [parse_named_path(item)[0] for item in heldout],
        "exclude_source_tokens": sorted(tokens),
        "input_rows": int(len(train)),
        "output_rows": int(len(out)),
        "removed_rows": int(len(train) - len(out)),
        "output": str(path),
    }
    path.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return out


def _interpret(pair_overlap: bool, source_overlap: bool, strict_source_overlap: bool) -> str:
    if pair_overlap:
        return "not_independent_same_canonical_pairs"
    if source_overlap and strict_source_overlap:
        return "not_independent_same_source_family"
    if source_overlap:
        return "pair_independent_but_source_overlap_sensitivity_only"
    return "pair_and_source_independent"

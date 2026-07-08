from __future__ import annotations

import hashlib
import json
import random
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.mechanism_panel_audit import filter_mechanism_panel
from analysis.ml.split_manifest import dataframe_content_hash, row_identity_hash


DEFAULT_SPLIT_MODES = ("drug_holdout", "target_holdout", "source_holdout", "random")
LABEL_COL = "mechanism_pu_label"


def _clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _safe_token(value: Any, fallback: str) -> str:
    text = _clean(value)
    return text if text else fallback


def _node_id(kind: str, value: Any, fallback: str) -> str:
    return f"{kind}:{_safe_token(value, fallback)}"


def _source_text(row: pd.Series) -> str:
    for col in ("mechanism_pu_label_source_all", "label_source", "source_family", "upstream_source"):
        text = _clean(row.get(col, ""))
        if text:
            return text
    return ""


def _edge_key(row: pd.Series, panel: str) -> str:
    key = _clean(row.get("canonical_pair_key", ""))
    if key:
        return key
    drug = _clean(row.get("drug_id", row.get("_pu_drug_key", ""))).lower()
    target = _clean(row.get("target_id", row.get("_pu_target_key", ""))).upper()
    adr = _clean(row.get("adr_id", row.get("adr_term", panel))).lower()
    return f"{drug}|{target}|{adr or panel}"


def _event_id(edge_key: str) -> str:
    digest = hashlib.sha1(edge_key.encode("utf-8")).hexdigest()[:16]
    return f"mechanism_event:{digest}"


def _label_frame(df: pd.DataFrame, panel: str, label_col: str) -> pd.DataFrame:
    work = filter_mechanism_panel(df, panel) if panel != "all" else df.copy()
    labels = pd.to_numeric(work.get(label_col, pd.Series(pd.NA, index=work.index)), errors="coerce")
    work = work.loc[labels.isin([0, 1])].copy()
    work[label_col] = labels.loc[work.index].astype(int)
    if work.empty:
        return work
    work["_graph_edge_key"] = work.apply(lambda row: _edge_key(row, panel), axis=1)
    work["_graph_event_node"] = work["_graph_edge_key"].map(_event_id)
    work["_graph_drug_node"] = [
        _node_id("drug", row.get("drug_id", row.get("_pu_drug_key", "")), f"row_{idx}")
        for idx, row in work.iterrows()
    ]
    work["_graph_target_node"] = [
        _node_id("target", row.get("target_id", row.get("_pu_target_key", "")), f"row_{idx}")
        for idx, row in work.iterrows()
    ]
    work["_graph_adr_node"] = [
        _node_id("adr_panel" if panel != "all" else "adr", panel if panel != "all" else row.get("adr_id", row.get("adr_term", "unknown")), "unknown")
        for _, row in work.iterrows()
    ]
    work["_graph_source"] = work.apply(_source_text, axis=1)
    return work.drop_duplicates("_graph_edge_key", keep="first").reset_index(drop=True)


def _node_rows(edges: pd.DataFrame, panel: str) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for _, row in edges.iterrows():
        candidates = [
            (row["_graph_drug_node"], "drug", row.get("drug_id", "")),
            (row["_graph_target_node"], "target", row.get("target_id", row.get("target_gene", ""))),
            (row["_graph_adr_node"], "adr_panel" if panel != "all" else "adr", panel if panel != "all" else row.get("adr_term", "")),
            (row["_graph_event_node"], "mechanism_hyperedge", row["_graph_edge_key"]),
        ]
        for node_id, node_type, label in candidates:
            if node_id in seen:
                continue
            seen.add(node_id)
            rows.append({"node_id": node_id, "node_type": node_type, "label": _clean(label), "source": "mechanism_pu_table"})
    return pd.DataFrame(rows)


def _hyperedges(edges: pd.DataFrame, label_col: str) -> pd.DataFrame:
    keep = [
        "_graph_edge_key",
        "_graph_event_node",
        "_graph_drug_node",
        "_graph_target_node",
        "_graph_adr_node",
        label_col,
        "mechanism_pu_state",
        "mechanism_pu_training_role",
        "_sample_weight",
        "_graph_source",
        "drug_id",
        "target_id",
        "pdb_id",
        "adr_site_group",
        "protein_class",
        "target_family",
        "label_source",
        "source_family",
        "upstream_source",
    ]
    cols = [col for col in keep if col in edges.columns]
    out = edges[cols].copy()
    out = out.rename(
        columns={
            "_graph_edge_key": "hyperedge_id",
            "_graph_event_node": "event_node_id",
            "_graph_drug_node": "drug_node_id",
            "_graph_target_node": "target_node_id",
            "_graph_adr_node": "adr_node_id",
            label_col: "label",
            "_graph_source": "evidence_source",
        }
    )
    out["label"] = pd.to_numeric(out["label"], errors="coerce").astype("Int64")
    return out


def _triples(hyperedges: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for _, row in hyperedges.iterrows():
        split_payload = {
            "hyperedge_id": row["hyperedge_id"],
            "label": int(row["label"]) if pd.notna(row["label"]) else None,
            "evidence_source": row.get("evidence_source", ""),
            "mechanism_pu_state": row.get("mechanism_pu_state", ""),
        }
        rows.extend(
            [
                {
                    "head": row["drug_node_id"],
                    "relation": "drug_in_mechanism",
                    "tail": row["event_node_id"],
                    "hyperedge_id": row["hyperedge_id"],
                    **split_payload,
                },
                {
                    "head": row["target_node_id"],
                    "relation": "target_in_mechanism",
                    "tail": row["event_node_id"],
                    "hyperedge_id": row["hyperedge_id"],
                    **split_payload,
                },
                {
                    "head": row["event_node_id"],
                    "relation": "mechanism_has_adr_panel",
                    "tail": row["adr_node_id"],
                    "hyperedge_id": row["hyperedge_id"],
                    **split_payload,
                },
            ]
        )
    return pd.DataFrame(rows)


def _group_col_for_mode(split_mode: str) -> str | None:
    return {
        "drug_holdout": "drug_id",
        "target_holdout": "target_id",
        "source_holdout": "label_source",
        "random": None,
    }.get(split_mode)


def _assign_splits(
    hyperedges: pd.DataFrame,
    *,
    split_mode: str,
    seed: int,
    validation_fraction: float,
    test_fraction: float,
) -> pd.Series:
    rng = random.Random(seed)
    split = pd.Series("train", index=hyperedges.index, dtype="object")
    if hyperedges.empty:
        return split
    group_col = _group_col_for_mode(split_mode)
    if group_col is None or group_col not in hyperedges.columns:
        indices = list(hyperedges.index)
        rng.shuffle(indices)
        n_test = max(1, round(len(indices) * test_fraction))
        n_val = max(1, round(len(indices) * validation_fraction)) if len(indices) >= 10 else 0
        split.loc[indices[:n_test]] = "test"
        split.loc[indices[n_test : n_test + n_val]] = "validation"
        return split
    groups = sorted({_safe_token(value, "missing") for value in hyperedges[group_col].fillna("").astype(str)})
    rng.shuffle(groups)
    n_test_groups = max(1, round(len(groups) * test_fraction)) if groups else 0
    n_val_groups = max(1, round(len(groups) * validation_fraction)) if len(groups) >= 10 else 0
    test_groups = set(groups[:n_test_groups])
    val_groups = set(groups[n_test_groups : n_test_groups + n_val_groups])
    values = hyperedges[group_col].fillna("").astype(str).map(lambda value: _safe_token(value, "missing"))
    split.loc[values.isin(test_groups)] = "test"
    split.loc[values.isin(val_groups)] = "validation"
    return split


def _split_manifest(
    hyperedges: pd.DataFrame,
    *,
    split_mode: str,
    seed: int,
    validation_fraction: float,
    test_fraction: float,
) -> pd.DataFrame:
    table = hyperedges.copy()
    table["split"] = _assign_splits(
        table,
        split_mode=split_mode,
        seed=seed,
        validation_fraction=validation_fraction,
        test_fraction=test_fraction,
    )
    table["fold"] = split_mode
    hash_cols = [col for col in ["hyperedge_id", "drug_id", "target_id", "adr_site_group", "label"] if col in table.columns]
    table["row_identity_hash"] = table.apply(lambda row: row_identity_hash(row, hash_cols), axis=1)
    return table


def _counts(table: pd.DataFrame) -> dict[str, Any]:
    if table.empty:
        return {"rows": 0}
    label = pd.to_numeric(table.get("label", pd.Series(pd.NA, index=table.index)), errors="coerce")
    out: dict[str, Any] = {
        "rows": int(len(table)),
        "positive_edges": int(label.eq(1).sum()),
        "negative_edges": int(label.eq(0).sum()),
    }
    if "split" in table.columns:
        for split_name, group in table.groupby("split", dropna=False):
            group_label = pd.to_numeric(group.get("label", pd.Series(pd.NA, index=group.index)), errors="coerce")
            out[f"{split_name}_rows"] = int(len(group))
            out[f"{split_name}_positive_edges"] = int(group_label.eq(1).sum())
            out[f"{split_name}_negative_edges"] = int(group_label.eq(0).sum())
    return out


def _write_negative_policy(out: Path, *, panel: str) -> None:
    policy = {
        "panel": panel,
        "default_negative_source": "Only rows already labeled mechanism_pu_label=0 are emitted as production negative edges.",
        "unknown_policy": "Rows with missing mechanism_pu_label remain unknown/background and are not written as negative KGE or hypergraph labels.",
        "corruption_policy": "Corruption-based negative sampling is disabled by default. If enabled later, it must be marked benchmark_or_sensitivity_only and generated within each training fold.",
        "leakage_policy": "Do not use Atlas score, SCORCH score, heatmap rank, FDR, or mechanism_graph_score to choose negative edges when those are model features.",
        "conflict_policy": "Excluded/ambiguous/conflicting rows are omitted from production graph labels.",
        "cardiac_qt_focus": "The cardiac_qt panel is the first claim candidate because current labels are site-group matched and comparatively dense.",
    }
    (out / "negative_sampling_policy.json").write_text(json.dumps(policy, indent=2, sort_keys=True), encoding="utf-8")


def build_mechanism_graph_training_splits(
    mechanism_pu_table_path: str | Path,
    out_dir: str | Path,
    *,
    panel: str = "cardiac_qt",
    split_modes: tuple[str, ...] | list[str] = DEFAULT_SPLIT_MODES,
    label_col: str = LABEL_COL,
    seed: int = 42,
    validation_fraction: float = 0.1,
    test_fraction: float = 0.2,
) -> dict[str, Any]:
    source = pd.read_csv(mechanism_pu_table_path, low_memory=False)
    edges = _label_frame(source, panel, label_col)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    nodes = _node_rows(edges, panel)
    hyperedges = _hyperedges(edges, label_col) if not edges.empty else pd.DataFrame()
    triples = _triples(hyperedges) if not hyperedges.empty else pd.DataFrame()
    nodes.to_csv(out / "graph_nodes.csv", index=False)
    hyperedges.to_csv(out / "mechanism_hyperedges.csv", index=False)
    triples.to_csv(out / "kge_reified_triples.csv", index=False)
    negative_edges = hyperedges.loc[pd.to_numeric(hyperedges.get("label", pd.Series(dtype="float")), errors="coerce").eq(0)].copy()
    negative_edges.to_csv(out / "negative_edges.csv", index=False)
    split_summaries: dict[str, Any] = {}
    for mode in split_modes:
        split_table = _split_manifest(
            hyperedges,
            split_mode=mode,
            seed=seed,
            validation_fraction=validation_fraction,
            test_fraction=test_fraction,
        ) if not hyperedges.empty else pd.DataFrame()
        mode_dir = out / f"split_{mode}"
        mode_dir.mkdir(parents=True, exist_ok=True)
        split_table.to_csv(mode_dir / "edge_split_manifest.csv", index=False)
        if not split_table.empty:
            split_triples = triples.merge(
                split_table[["hyperedge_id", "split", "fold", "row_identity_hash"]],
                on="hyperedge_id",
                how="left",
            )
        else:
            split_triples = triples.copy()
        split_triples.to_csv(mode_dir / "kge_reified_triples.split.csv", index=False)
        split_summary = _counts(split_table)
        split_summary.update({"split_mode": mode, "manifest": str(mode_dir / "edge_split_manifest.csv")})
        (mode_dir / "edge_split_manifest.json").write_text(
            json.dumps(
                {
                    "dataset_hash": dataframe_content_hash(hyperedges) if not hyperedges.empty else None,
                    "split_mode": mode,
                    "seed": seed,
                    "validation_fraction": validation_fraction,
                    "test_fraction": test_fraction,
                    "summary": split_summary,
                },
                indent=2,
                sort_keys=True,
            ),
            encoding="utf-8",
        )
        split_summaries[mode] = split_summary
    _write_negative_policy(out, panel=panel)
    manifest = {
        "source_table": str(mechanism_pu_table_path),
        "out_dir": str(out),
        "panel": panel,
        "label_col": label_col,
        "nodes": str(out / "graph_nodes.csv"),
        "hyperedges": str(out / "mechanism_hyperedges.csv"),
        "triples": str(out / "kge_reified_triples.csv"),
        "negative_edges": str(out / "negative_edges.csv"),
        "negative_sampling_policy": str(out / "negative_sampling_policy.json"),
        "counts": {
            "source_rows": int(len(source)),
            "panel_labelable_edges": int(len(hyperedges)),
            "nodes": int(len(nodes)),
            "triples": int(len(triples)),
            "negative_edges": int(len(negative_edges)),
        },
        "split_summaries": split_summaries,
        "policy": "KGE/hypergraph training uses reified drug-target-ADR panel hyperedges. Unknowns are not production negatives; negative labels come only from explicit mechanism_pu_label=0 rows.",
    }
    (out / "graph_training_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest

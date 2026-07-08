from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.feature_metadata import enrich_ml_feature_metadata
from analysis.ml.mechanism_panel_audit import filter_mechanism_panel, write_panel_readiness_audit
from analysis.ml.mechanism_graph_splits import DEFAULT_SPLIT_MODES, build_mechanism_graph_training_splits
from analysis.ml.mechanism_pu_tables import build_mechanism_pu_table
from analysis.ml.mechanism_recovery import DEFAULT_SCORE_COLS, write_mechanism_topk_recovery
from analysis.ml.train_classifier import train_ml_model


def _source_values(df: pd.DataFrame, source_col: str) -> list[str]:
    if source_col not in df.columns:
        return []
    values: set[str] = set()
    for raw in df[source_col].dropna().astype(str):
        for token in raw.split(";"):
            token = token.strip()
            if token:
                values.add(token)
    return sorted(values)


def _source_mask(df: pd.DataFrame, source_col: str, source: str) -> pd.Series:
    values = df[source_col].fillna("").astype(str)
    return values.eq(source) | values.str.split(";").map(lambda parts: source in parts)


def _label_counts(df: pd.DataFrame, label_col: str) -> dict[str, int]:
    label = pd.to_numeric(df.get(label_col, pd.Series(pd.NA, index=df.index)), errors="coerce")
    return {
        "rows": int(len(df)),
        "labelable_rows": int(label.notna().sum()),
        "positive_rows": int(label.eq(1).sum()),
        "negative_rows": int(label.eq(0).sum()),
        "unknown_rows": int(label.isna().sum()),
    }


def run_leave_one_source_out(
    pu_table_path: str | Path,
    out_dir: str | Path,
    *,
    source_col: str = "label_source",
    label_col: str = "mechanism_pu_label",
    feature_set: str = "spd_mechanism_nonleaky",
    model_type: str = "logistic_regression",
    pu_mode: str = "stratified_bagging_pu",
    pu_bags: int = 20,
    seed: int = 42,
    min_train_positives: int = 5,
    min_train_negatives: int = 5,
    min_test_positives: int = 5,
    min_test_negatives: int = 5,
    max_sources: int = 0,
    strict_feature_set: bool = False,
) -> pd.DataFrame:
    df = pd.read_csv(pu_table_path, low_memory=False)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    sources = _source_values(df, source_col)
    if max_sources and max_sources > 0:
        sources = sources[:max_sources]
    rows: list[dict[str, Any]] = []
    temp_dir = out / "_source_holdout_tables"
    temp_dir.mkdir(parents=True, exist_ok=True)
    for source in sources:
        heldout = _source_mask(df, source_col, source)
        label = pd.to_numeric(df.get(label_col, pd.Series(pd.NA, index=df.index)), errors="coerce")
        # Prevent held-out-source unknowns from entering the fold-local PU background.
        temp = df[~(heldout & label.isna())].copy()
        train_counts = _label_counts(temp[~_source_mask(temp, source_col, source)], label_col)
        test_counts = _label_counts(temp[_source_mask(temp, source_col, source)], label_col)
        row: dict[str, Any] = {
            "heldout_source": source,
            "status": "pending",
            "train_labelable_rows": train_counts["labelable_rows"],
            "train_positive_rows": train_counts["positive_rows"],
            "train_negative_rows": train_counts["negative_rows"],
            "train_unknown_rows": train_counts["unknown_rows"],
            "test_labelable_rows": test_counts["labelable_rows"],
            "test_positive_rows": test_counts["positive_rows"],
            "test_negative_rows": test_counts["negative_rows"],
        }
        if test_counts["positive_rows"] < min_test_positives or test_counts["negative_rows"] < min_test_negatives:
            row["status"] = "skipped_insufficient_heldout_labels"
            rows.append(row)
            continue
        if train_counts["positive_rows"] < min_train_positives or train_counts["negative_rows"] < min_train_negatives:
            row["status"] = "skipped_insufficient_train_labels"
            rows.append(row)
            continue
        safe_name = "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in source)[:80]
        temp_path = temp_dir / f"mechanism_pu_without_heldout_unknowns_{safe_name}.csv"
        temp.to_csv(temp_path, index=False)
        model_dir = out / "leave_one_source_out" / safe_name
        try:
            trained = train_ml_model(
                temp_path,
                label_col=label_col,
                feature_set=feature_set,
                model_type=model_type,
                split_mode="source_holdout",
                split_column=source_col,
                test_values=[source],
                out_dir=model_dir,
                seed=seed,
                pu_mode=pu_mode,
                pu_bags=pu_bags,
                n_bootstraps=100,
                claim_mode="exploratory",
                dataset_provenance={
                    "stage": "mechanism_leave_one_source_out_pu",
                    "heldout_source": source,
                    "heldout_source_unknowns_excluded_from_pu_pool": True,
                },
                strict_feature_set=strict_feature_set,
            )
        except Exception as exc:  # noqa: BLE001 - surfaced in summary table.
            row["status"] = "failed"
            row["error"] = str(exc)
            rows.append(row)
            continue
        metrics = dict(trained.get("metrics", {}))
        metric_keys = [
            "AUROC",
            "AUPRC",
            "Brier",
            "ECE",
            "F1@0.5",
            "precision@0.5",
            "recall@0.5",
            "specificity@0.5",
            "accuracy@0.5",
            "balanced_accuracy@0.5",
        ]
        row.update({"status": "trained", "model_dir": str(model_dir)})
        row.update({key: metrics.get(key) for key in metric_keys})
        rows.append(row)
    result = pd.DataFrame(rows)
    result.to_csv(out / "mechanism_leave_one_source_out_summary.csv", index=False)
    return result



def _write_full_feature_table(
    table_path: Path,
    out_path: Path,
    *,
    run_dir: Path | None,
    repo_root: Path,
    chemical_cluster: str,
    target_family: str,
    source_lineage: str,
) -> dict[str, Any]:
    df = pd.read_csv(table_path, low_memory=False)
    enriched, summary = enrich_ml_feature_metadata(
        df,
        chemical_cluster=chemical_cluster,
        target_family=target_family,
        source_lineage=source_lineage,
        run_dir=run_dir,
        repo_root=repo_root,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    enriched.to_csv(out_path, index=False)
    summary.update({"input": str(table_path), "output": str(out_path), "columns": int(len(enriched.columns))})
    out_path.with_suffix(".manifest.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=str),
        encoding="utf-8",
    )
    return summary

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build mechanism PU labels and source-held/top-K recovery reports.")
    parser.add_argument("--mechanism-table", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--source-col", default="label_source")
    parser.add_argument("--label-col", default="mechanism_pu_label")
    parser.add_argument("--feature-set", default="spd_mechanism_nonleaky")
    parser.add_argument("--model", default="logistic_regression")
    parser.add_argument(
        "--models",
        nargs="*",
        default=None,
        help="Optional model-family comparison list; overrides --model when provided.",
    )
    parser.add_argument("--panel", default="all", help="ADR/site panel filter, e.g. cardiac_qt, cns_sedation, liver, gi.")
    parser.add_argument("--pu-mode", default="stratified_bagging_pu")
    parser.add_argument("--pu-bags", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--score-cols", nargs="*", default=None)
    parser.add_argument("--group-cols", nargs="*", default=["target_id", "pdb_id", "adr_site_group"])
    parser.add_argument("--top-k", nargs="*", type=int, default=[5, 10, 20, 50])
    parser.add_argument("--skip-source-held", action="store_true")
    parser.add_argument("--write-graph-splits", action="store_true")
    parser.add_argument("--graph-split-modes", nargs="*", default=list(DEFAULT_SPLIT_MODES))
    parser.add_argument("--max-sources", type=int, default=0)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--run-dir", type=Path, default=None, help="Atlas run data directory used to backfill full ML features.")
    parser.add_argument("--feature-metadata", choices=["auto", "always", "never"], default="auto")
    parser.add_argument("--chemical-cluster", choices=["auto", "scaffold", "smiles", "chemotype", "none"], default="auto")
    parser.add_argument("--target-family", choices=["auto", "protein_class", "gene_heuristic", "none"], default="auto")
    parser.add_argument("--source-lineage", choices=["auto", "none"], default="auto")
    parser.add_argument("--strict-feature-set", action="store_true")
    parser.add_argument("--min-train-positives", type=int, default=5)
    parser.add_argument("--min-train-negatives", type=int, default=5)
    parser.add_argument("--min-test-positives", type=int, default=5)
    parser.add_argument("--min-test-negatives", type=int, default=5)
    args = parser.parse_args(argv)

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    pu_table = out / "mechanism_pu_table.csv"
    build_mechanism_pu_table(args.mechanism_table, pu_table)
    repo_root = args.repo_root.resolve()
    feature_metadata_summary: dict[str, Any] = {}
    if args.feature_metadata in {"auto", "always"}:
        full_feature_table = out / "mechanism_pu_table.full_features.csv"
        feature_metadata_summary = _write_full_feature_table(
            pu_table,
            full_feature_table,
            run_dir=args.run_dir,
            repo_root=repo_root,
            chemical_cluster=args.chemical_cluster,
            target_family=args.target_family,
            source_lineage=args.source_lineage,
        )
        pu_table = full_feature_table
    write_panel_readiness_audit(
        pu_table,
        out / "mechanism_adr_panel_readiness.csv",
        label_col=args.label_col,
        source_balance_path=out / "mechanism_adr_panel_source_balance.csv",
    )
    eval_table = pu_table
    if args.panel != "all":
        panel_df = filter_mechanism_panel(pd.read_csv(pu_table, low_memory=False), args.panel)
        suffix = ".full_features" if args.feature_metadata in {"auto", "always"} else ""
        eval_table = out / f"mechanism_pu_table.panel_{args.panel}{suffix}.csv"
        panel_df.to_csv(eval_table, index=False)
    score_cols = args.score_cols if args.score_cols else DEFAULT_SCORE_COLS
    topk = write_mechanism_topk_recovery(
        eval_table,
        out / (f"mechanism_topk_recovery.panel_{args.panel}.csv" if args.panel != "all" else "mechanism_topk_recovery.csv"),
        label_col=args.label_col,
        score_cols=score_cols,
        group_cols=args.group_cols,
        ks=args.top_k,
    )
    loso = pd.DataFrame()
    model_list = args.models if args.models else [args.model]
    if not args.skip_source_held:
        model_frames: list[pd.DataFrame] = []
        for model_name in model_list:
            model_out = out if len(model_list) == 1 else out / "model_family_comparison" / model_name
            frame = run_leave_one_source_out(
                eval_table,
                model_out,
                source_col=args.source_col,
                label_col=args.label_col,
                feature_set=args.feature_set,
                model_type=model_name,
                pu_mode=args.pu_mode,
                pu_bags=args.pu_bags,
                seed=args.seed,
                min_train_positives=args.min_train_positives,
                min_train_negatives=args.min_train_negatives,
                min_test_positives=args.min_test_positives,
                min_test_negatives=args.min_test_negatives,
                max_sources=args.max_sources,
                strict_feature_set=bool(args.strict_feature_set),
            )
            frame.insert(0, "model_type", model_name)
            model_frames.append(frame)
        if model_frames:
            loso = pd.concat(model_frames, ignore_index=True)
            if len(model_list) > 1:
                loso.to_csv(out / "mechanism_model_family_comparison.csv", index=False)
                loso.to_csv(out / "mechanism_leave_one_source_out_summary.csv", index=False)
    graph_manifest: dict[str, Any] = {}
    if args.write_graph_splits:
        graph_manifest = build_mechanism_graph_training_splits(
            eval_table,
            out / (f"graph_splits.panel_{args.panel}" if args.panel != "all" else "graph_splits"),
            panel="all",
            split_modes=args.graph_split_modes,
            label_col=args.label_col,
            seed=args.seed,
        )
    manifest = {
        "mechanism_table": str(args.mechanism_table),
        "out_dir": str(out),
        "pu_table": str(pu_table),
        "eval_table": str(eval_table),
        "feature_metadata": feature_metadata_summary,
        "panel": args.panel,
        "models": model_list,
        "panel_readiness": str(out / "mechanism_adr_panel_readiness.csv"),
        "topk_rows": int(len(topk)),
        "leave_one_source_rows": int(len(loso)),
        "graph_splits": graph_manifest,
        "policy": "Mechanism PU/source-held evaluation keeps unknowns separate and excludes held-out-source unknowns from fold-local PU background.",
    }
    (out / "mechanism_pu_evaluation_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

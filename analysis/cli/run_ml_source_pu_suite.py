from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from analysis.io import load_config
from analysis.ml.audit_suite import run_ml_audit_suite
from analysis.ml.feature_metadata import refresh_tables
from analysis.ml.source_benchmark_tables import build_source_benchmark_tables
from analysis.ml.train_classifier import train_ml_model


SOURCE_EXPERTS: list[dict[str, Any]] = [
    {
        "name": "curated_bioactivity_within_source",
        "rel_path": "source_benchmark_comparison_v2/tables/model_ready/ml_curated_bioactivity_model_ready.csv",
        "label": "curated_bioactivity_label",
        "feature_set": "pilot_binding_only",
        "split": "drug_holdout",
    },
    {
        "name": "toxcast_hts_within_source",
        "rel_path": "source_benchmark_comparison_v2/tables/model_ready/ml_toxcast_hts_model_ready.csv",
        "label": "toxcast_hts_label",
        "feature_set": "pilot_binding_only",
        "split": "drug_holdout",
    },
    {
        "name": "spd_exposure_within_source",
        "rel_path": "source_benchmark_comparison_v2/tables/model_ready/ml_spd_exposure_model_ready.csv",
        "label": "spd_exposure_label",
        "feature_set": "spd_exposure_nonleaky",
        "split": "drug_holdout",
    },
    {
        "name": "curated_to_toxcast_hard_transfer",
        "rel_path": "source_benchmark_comparison_v2/tables/ml_source_transfer_table.csv",
        "label": "source_specific_activity_label",
        "feature_set": "pilot_binding_only",
        "split": "custom_source_holdout",
        "split_column": "source_objective",
        "train_values": ["curated_bioactivity"],
        "test_values": ["toxcast_hts"],
    },
    {
        "name": "curated_to_toxcast_matched_transfer",
        "rel_path": "source_benchmark_comparison_v2/tables/ml_matched_source_transfer_table.csv",
        "label": "source_specific_activity_label",
        "feature_set": "pilot_binding_only",
        "split": "custom_source_holdout",
        "split_column": "source_objective",
        "train_values": ["curated_bioactivity"],
        "test_values": ["toxcast_hts"],
    },
]

PU_MODES = [
    "standard_binary",
    "bagging_pu",
    "stratified_bagging_pu",
    "propensity_weighted_pu",
    "elkan_noto",
    "pulsnar_style",
]


def _run_train_and_audit(
    *,
    dataset: Path,
    label: str,
    feature_set: str,
    model_type: str,
    split: str,
    out_dir: Path,
    seed: int,
    bootstraps: int,
    pu_mode: str,
    class_weight: str | None,
    split_column: str | None = None,
    train_values: list[str] | None = None,
    test_values: list[str] | None = None,
    pu_bags: int = 50,
    pu_unlabeled_ratio: float = 1.0,
    claim_mode: str = "exploratory",
    group_reweight_cols: list[str] | None = None,
    group_reweight_include_label: bool = True,
    group_reweight_max_factor: float = 5.0,
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    result: dict[str, Any] = {
        "dataset": str(dataset),
        "label": label,
        "feature_set": feature_set,
        "split": split,
        "pu_mode": pu_mode,
        "status": "pending",
    }
    try:
        trained = train_ml_model(
            dataset,
            label,
            feature_set,
            model_type,
            split,
            out_dir / "model",
            seed,
            split_column=split_column,
            train_values=train_values,
            test_values=test_values,
            n_bootstraps=bootstraps,
            pu_mode=pu_mode,
            class_weight=class_weight,
            nested_model_selection=True,
            pu_bags=pu_bags,
            pu_unlabeled_ratio=pu_unlabeled_ratio,
            claim_mode=claim_mode,
            repo_root=Path.cwd(),
            group_reweight_cols=group_reweight_cols,
            group_reweight_include_label=group_reweight_include_label,
            group_reweight_max_factor=group_reweight_max_factor,
        )
        audit = run_ml_audit_suite(
            dataset,
            label,
            out_dir / "audit",
            feature_set=feature_set,
            source_col="label_source",
            predictions_path=out_dir / "model" / "model_predictions.csv",
            source_transfer_train=train_values,
            source_transfer_test=test_values,
            seed=seed,
            claim_mode=claim_mode,
        )
        result.update(
            {
                "status": "trained",
                "metrics": trained.get("metrics", {}),
                "n_train": trained.get("n_train"),
                "n_test": trained.get("n_test"),
                "features": trained.get("features", []),
                "model_dir": str(out_dir / "model"),
                "audit_dir": str(out_dir / "audit"),
                "audit_claim_readiness": audit.get("claim_readiness", {}),
            }
        )
    except Exception as exc:
        result["status"] = "failed"
        result["error"] = str(exc)
    return result


def _maybe_build_tables(args: argparse.Namespace, run_dir: Path) -> None:
    if not args.bioactivity_source:
        return
    tables_dir = run_dir / "source_benchmark_comparison_v2" / "tables"
    build_source_benchmark_tables(
        args.bioactivity_source,
        tables_dir,
        spd_table_path=args.spd_table,
    )


def _metadata_refreshed_dataset(
    dataset: Path,
    *,
    run_dir: Path,
    out_dir: Path,
    mode: str,
    chemical_cluster: str,
    target_family: str,
    source_lineage: str,
    drop_columns: list[str] | None,
) -> Path:
    if mode == "never":
        return dataset
    manifest = refresh_tables(
        [dataset],
        out_dir=out_dir / "_feature_metadata",
        in_place=False,
        run_dir=run_dir,
        chemical_cluster=chemical_cluster,
        target_family=target_family,
        source_lineage=source_lineage,
        drop_columns=drop_columns,
    )
    outputs = manifest.get("outputs") or []
    return Path(outputs[0]) if outputs else dataset


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run Atlas source-specific ML, hard source-transfer stress tests, "
            "and PU sensitivity models from one command."
        )
    )
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--run-dir", type=Path, default=Path("data/pilotstudy"))
    parser.add_argument("--bioactivity-source", type=Path, default=None)
    parser.add_argument("--spd-table", type=Path, default=None)
    parser.add_argument("--binding-dataset", type=Path, default=None)
    parser.add_argument("--binding-label", default="four_state_ml_label")
    parser.add_argument("--binding-feature-set", default="consensus_banana_binding")
    parser.add_argument("--model", default="logistic_regression")
    parser.add_argument("--bootstraps", type=int, default=100)
    parser.add_argument("--pu-bags", type=int, default=50)
    parser.add_argument("--pu-unlabeled-ratio", type=float, default=1.0)
    parser.add_argument("--class-weight", choices=["balanced", "none"], default="balanced")
    parser.add_argument(
        "--feature-metadata",
        choices=["auto", "always", "never"],
        default="auto",
        help="Refresh chemical_cluster, target_family, and source-lineage columns before training.",
    )
    parser.add_argument("--chemical-cluster", choices=["auto", "scaffold", "smiles", "chemotype", "none"], default="auto")
    parser.add_argument("--target-family", choices=["auto", "protein_class", "gene_heuristic", "none"], default="auto")
    parser.add_argument("--source-lineage", choices=["auto", "none"], default="auto")
    parser.add_argument(
        "--drop-column",
        nargs="*",
        default=None,
        help="Optional columns to remove from refreshed training tables for sensitivity passes.",
    )
    parser.add_argument("--skip-source-suite", action="store_true")
    parser.add_argument("--skip-pu-suite", action="store_true")
    parser.add_argument("--pu-modes", nargs="*", choices=PU_MODES, default=PU_MODES)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--claim-mode", choices=["exploratory", "publication"], default="exploratory")
    parser.add_argument(
        "--group-reweight-cols",
        nargs="*",
        default=None,
        help="Optional capped inverse-frequency reweighting over source/family strata.",
    )
    parser.add_argument(
        "--group-reweight-include-label",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Include the supervised label in group-reweight strata.",
    )
    parser.add_argument("--group-reweight-max-factor", type=float, default=5.0)
    args = parser.parse_args(argv)

    run_dir = args.run_dir.resolve()
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    config = load_config(args.config)
    seed = int(config.get("project", {}).get("random_seed", 42))
    class_weight = None if args.class_weight == "none" else args.class_weight
    _maybe_build_tables(args, run_dir)

    results: list[dict[str, Any]] = []
    if not args.skip_source_suite:
        for spec in SOURCE_EXPERTS:
            dataset = run_dir / spec["rel_path"]
            if not dataset.exists():
                results.append(
                    {
                        "suite": "source_specific",
                        "name": spec["name"],
                        "status": "skipped",
                        "reason": "dataset_missing",
                        "dataset": str(dataset),
                    }
                )
                continue
            if args.feature_metadata in {"auto", "always"}:
                dataset = _metadata_refreshed_dataset(
                    dataset,
                    run_dir=run_dir,
                    out_dir=out,
                    mode=args.feature_metadata,
                    chemical_cluster=args.chemical_cluster,
                    target_family=args.target_family,
                    source_lineage=args.source_lineage,
                    drop_columns=args.drop_column,
                )
            results.append(
                {
                    "suite": "source_specific",
                    "name": spec["name"],
                    **_run_train_and_audit(
                        dataset=dataset,
                        label=spec["label"],
                        feature_set=spec["feature_set"],
                        model_type=args.model,
                        split=spec["split"],
                        out_dir=out / "source_specific" / spec["name"],
                        seed=seed,
                        bootstraps=args.bootstraps,
                        pu_mode="standard_binary",
                        class_weight=class_weight,
                        split_column=spec.get("split_column"),
                        train_values=spec.get("train_values"),
                        test_values=spec.get("test_values"),
                        claim_mode=args.claim_mode,
                        group_reweight_cols=args.group_reweight_cols,
                        group_reweight_include_label=bool(args.group_reweight_include_label),
                        group_reweight_max_factor=float(args.group_reweight_max_factor),
                    ),
                }
            )

    if not args.skip_pu_suite:
        binding_dataset = args.binding_dataset or (
            run_dir / "moe_experts" / "banana_four_state_with_scorch_backfill_model_ready.csv"
        )
        if not binding_dataset.exists():
            results.append(
                {
                    "suite": "pu_sensitivity",
                    "status": "skipped",
                    "reason": "binding_dataset_missing",
                    "dataset": str(binding_dataset),
                }
            )
        else:
            if args.feature_metadata in {"auto", "always"}:
                binding_dataset = _metadata_refreshed_dataset(
                    binding_dataset,
                    run_dir=run_dir,
                    out_dir=out,
                    mode=args.feature_metadata,
                    chemical_cluster=args.chemical_cluster,
                    target_family=args.target_family,
                    source_lineage=args.source_lineage,
                    drop_columns=args.drop_column,
                )
            for mode in args.pu_modes:
                results.append(
                    {
                        "suite": "pu_sensitivity",
                        "name": mode,
                        **_run_train_and_audit(
                            dataset=binding_dataset,
                            label=args.binding_label,
                            feature_set=args.binding_feature_set,
                            model_type=args.model,
                            split="drug_holdout",
                            out_dir=out / "pu_sensitivity" / mode,
                            seed=seed,
                            bootstraps=args.bootstraps,
                            pu_mode=mode,
                            class_weight=class_weight,
                            pu_bags=args.pu_bags,
                            pu_unlabeled_ratio=args.pu_unlabeled_ratio,
                            claim_mode=args.claim_mode,
                            group_reweight_cols=args.group_reweight_cols,
                            group_reweight_include_label=bool(args.group_reweight_include_label),
                            group_reweight_max_factor=float(args.group_reweight_max_factor),
                        ),
                    }
                )

    manifest = {
        "run_dir": str(run_dir),
        "out_dir": str(out),
        "policy": [
            "Within-source models are reported separately.",
            "ChEMBL/Papyrus to ToxCast is a hard transfer stress test, not the main claim.",
            "PU modes preserve measured inactives separately from unlabeled rows when the input table retains unknown labels.",
        ],
        "group_reweighting": {
            "group_reweight_cols": args.group_reweight_cols or [],
            "include_label": bool(args.group_reweight_include_label),
            "max_factor": float(args.group_reweight_max_factor),
        },
        "feature_metadata": {
            "mode": args.feature_metadata,
            "chemical_cluster": args.chemical_cluster,
            "target_family": args.target_family,
            "source_lineage": args.source_lineage,
            "drop_columns": args.drop_column or [],
        },
        "results": results,
    }
    (out / "ml_source_pu_suite_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

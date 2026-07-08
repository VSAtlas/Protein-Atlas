from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from analysis.io import load_config
from analysis.ml.train_classifier import train_ml_model


def _alias(value: Any, aliases: dict[str, str]) -> str:
    raw = str(value)
    return aliases.get(raw, raw)


def _list_config(value: Any) -> list[str] | None:
    if value is None:
        return None
    if isinstance(value, list):
        return [str(item) for item in value]
    return [str(value)]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run one config-driven Atlas ML model update and write the model-run ledger.")
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dataset", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--notes", default="")
    parser.add_argument("--allow-missing-feature-set", action="store_true", help="Allow exploratory runs to use the available subset of a named feature set.")
    args = parser.parse_args(argv)
    config = load_config(args.config)
    dataset = args.dataset or Path(str(config.get("dataset", "")))
    out_dir = args.out_dir or Path(str(config.get("out_dir", "outputs/models/model_update")))
    if not str(dataset):
        raise ValueError("dataset must be supplied in the config or with --dataset")
    task = str(config.get("task", "model_update"))
    feature_set = _alias(config["feature_set"], {"clean_binding_v1": "spd_binding_nonleaky"})
    split = _alias(config.get("split", "drug_holdout"), {"scaffold_grouped_cv": "scaffold_holdout"})
    pu_strategy = _alias(config.get("pu_strategy", "standard_binary"), {"supervised": "standard_binary"})
    train_ml_model(
        dataset,
        str(config["label"]),
        feature_set,
        str(config.get("model", "logistic_regression")),
        split,
        out_dir,
        int(config.get("seed", 42)),
        pu_mode=pu_strategy,
        exclude_features=_list_config(config.get("excluded_columns")),
        calibration_method=str(config.get("calibration", "none")),
        claim_mode=str(config.get("claim_mode", "exploratory")),
        nested_model_selection=bool(config.get("nested_model_selection", True)),
        strict_feature_set=(not args.allow_missing_feature_set and bool(config.get("strict_feature_set", True))),
        repo_root=args.repo_root,
        dataset_provenance={
            "stage": "run_model_update",
            "task": task,
            "config_path": str(args.config),
            "metrics_requested": config.get("metrics", []),
            "notes": args.notes,
            "run_id": args.run_id,
        },
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

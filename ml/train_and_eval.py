#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime, timezone
from pathlib import Path
import sys

import joblib
import numpy as np
import pandas as pd

if __package__ in {None, ""}:
    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

from ml.config import config_to_dict, load_atlas_cfg, load_config
from ml.data.bigbind import load_train_holdout_from_bigbind
from ml.evaluate import evaluate_holdout_metrics
from ml.featurize import BigBindFeaturizer
from ml.fpocket_bigbind import merge_fpocket_metrics_on_pocket, precompute_fpocket_for_bigbind_df
from ml.modeling import train_logistic_regression


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "CPU-only BigBind training workflow: load activities, featurize, "
            "train logistic regression, and evaluate holdout metrics."
        )
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to ml config file (JSON or YAML).",
    )
    return parser


def _setup_logging() -> logging.Logger:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    return logging.getLogger("ml.train_and_eval")


def _resolve_run_id(config_run_id: str | None) -> str:
    if config_run_id:
        return config_run_id
    return datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")


def _save_outputs(
    *,
    output_dir: Path,
    model: object,
    config: dict[str, object],
    featurizer_metadata: dict[str, object],
    metrics: dict[str, float | int],
    holdout_predictions: pd.DataFrame,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(model, output_dir / "model.joblib")

    (output_dir / "config_snapshot.txt").write_text(
        json.dumps(config, indent=2),
        encoding="utf-8",
    )
    (output_dir / "featurizer_metadata.json").write_text(
        json.dumps(featurizer_metadata, indent=2),
        encoding="utf-8",
    )
    (output_dir / "metrics_report.json").write_text(
        json.dumps(metrics, indent=2),
        encoding="utf-8",
    )
    pd.DataFrame([metrics]).to_csv(output_dir / "metrics_report.csv", index=False)
    holdout_predictions.to_csv(output_dir / "holdout_predictions.csv", index=False)


def run_pipeline(config_path: str) -> tuple[Path, dict[str, float | int]]:
    logger = _setup_logging()
    config = load_config(config_path)
    atlas_cfg = load_atlas_cfg(config)
    np.random.seed(config.random_seed)
    run_id = _resolve_run_id(config.run_id)
    output_dir = Path(__file__).resolve().parent / "outputs" / run_id
    output_dir.mkdir(parents=True, exist_ok=True)

    dataset = load_train_holdout_from_bigbind(
        bigbind_dir=config.bigbind_dir,
        train_pdb=config.train_pdb,
        splits=config.splits,
        max_rows=config.max_rows,
        random_seed=config.random_seed,
        dataset_split_mode=config.dataset_split_mode,
        exclude_target_prefixes=config.exclude_target_prefixes,
        logger=logger,
    )
    train_df = dataset.train_df.copy()
    holdout_df = dataset.test_df.copy()
    if config.features.pocket_fpocket:
        combined = pd.concat([train_df, holdout_df], ignore_index=False)
        metrics_df = precompute_fpocket_for_bigbind_df(
            df=combined,
            bigbind_root=dataset.bigbind_root,
            atlas_cfg=atlas_cfg,
            run_dir=output_dir,
        )
        train_df = merge_fpocket_metrics_on_pocket(train_df, metrics_df)
        holdout_df = merge_fpocket_metrics_on_pocket(holdout_df, metrics_df)

    featurizer = BigBindFeaturizer(
        bigbind_root=dataset.bigbind_root,
        features=config.features,
        atlas_cfg=atlas_cfg,
        fpocket_center_columns=config.fpocket_center_columns,
        fpocket_variant_column=config.fpocket_variant_column,
        fpocket_ph_column=config.fpocket_ph_column,
        fpocket_centers_by_pdb=config.fpocket_centers_by_pdb,
        logger=logger,
    )
    train_features = featurizer.transform(
        train_df,
        audit_dir=output_dir,
        audit_tag="train",
    )
    holdout_features = featurizer.transform(
        holdout_df,
        audit_dir=output_dir,
        audit_tag="holdout",
    )

    logger.info(
        "[ml] Featurized train rows=%d dropped_invalid_smiles=%d",
        train_features.X.shape[0],
        train_features.dropped_invalid_smiles,
    )
    logger.info(
        "[ml] Featurized holdout rows=%d dropped_invalid_smiles=%d",
        holdout_features.X.shape[0],
        holdout_features.dropped_invalid_smiles,
    )
    logger.info(
        "[ml] fpocket loaded=%d missing_center=%d missing_info=%d",
        featurizer.loaded_fpocket_count,
        featurizer.missing_center_count,
        featurizer.missing_info_count,
    )

    model = train_logistic_regression(
        X_train=train_features.X,
        y_train=train_features.y,
        model_config=config.model,
        random_seed=config.random_seed,
    )
    holdout_prob_active = model.predict_proba(holdout_features.X)[:, 1]

    metrics = evaluate_holdout_metrics(
        holdout_features.y,
        holdout_prob_active,
        fractions=(0.01, 0.02, 0.05, 0.10),
    )

    holdout_used = holdout_df.loc[holdout_features.source_index].copy()
    holdout_predictions = holdout_used[
        ["split", "lig_smiles", "active", "pocket", "ex_rec_pdb"]
    ].copy()
    holdout_predictions["murcko_scaffold"] = holdout_features.murcko_scaffolds
    holdout_predictions["p_active"] = holdout_prob_active

    _save_outputs(
        output_dir=output_dir,
        model=model,
        config=config_to_dict(config),
        featurizer_metadata=featurizer.metadata(),
        metrics=metrics,
        holdout_predictions=holdout_predictions,
    )

    return output_dir, metrics


def main() -> None:
    parser = _build_arg_parser()
    args = parser.parse_args()
    output_dir, metrics = run_pipeline(args.config)
    print(
        "[ml] "
        f"saved={output_dir} "
        f"N={metrics['N']} "
        f"n_actives={metrics['n_actives']} "
        f"PR_AUC={metrics['PR_AUC']:.6f} "
        f"EF@1%={metrics['EF@1%']:.6f} "
        f"EF@5%={metrics['EF@5%']:.6f}"
    )


if __name__ == "__main__":
    main()

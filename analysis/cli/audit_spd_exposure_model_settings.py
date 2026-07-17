from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re
from typing import Any

import pandas as pd

from analysis.ml.model_run_ledger import LEDGER_COLUMNS
from analysis.ml.spd_exposure_grouped_oof import (
    run_spd_exposure_grouped_oof,
    save_spd_exposure_group_folds,
)


_DEFAULT_RIDGE_ALPHAS = (0.1, 1.0, 10.0)
_DEFAULT_LIGHTGBM_SEEDS = (42, 137, 314)
_DEFAULT_LIGHTGBM_SETTINGS: tuple[dict[str, Any], ...] = (
    {
        "id": "lightgbm_regularized_shallow",
        "params": {
            "n_estimators": 200,
            "learning_rate": 0.04,
            "max_depth": 3,
            "num_leaves": 7,
            "min_child_samples": 15,
            "min_child_weight": 0.001,
            "subsample": 0.8,
            "subsample_freq": 1,
            "colsample_bytree": 0.8,
            "reg_alpha": 0.1,
            "reg_lambda": 2.0,
        },
    },
    {
        "id": "lightgbm_regularized_strong",
        "params": {
            "n_estimators": 300,
            "learning_rate": 0.03,
            "max_depth": 3,
            "num_leaves": 7,
            "min_child_samples": 25,
            "min_child_weight": 0.001,
            "subsample": 0.8,
            "subsample_freq": 1,
            "colsample_bytree": 0.75,
            "reg_alpha": 0.5,
            "reg_lambda": 5.0,
        },
    },
)
_BASE_METRICS = (
    "potency_exact_r2",
    "potency_exact_median_fold_error",
    "free_cmax_r2",
    "free_cmax_median_fold_error",
    "prevalence",
    "AUPRC",
    "AUROC",
    "Brier",
    "null_Brier",
    "Brier_skill",
    "ECE",
)
_LOWER_IS_BETTER = frozenset(
    {
        "potency_exact_median_fold_error",
        "free_cmax_median_fold_error",
        "Brier",
        "ECE",
    }
)


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return parsed


def _thread_count(value: str) -> int:
    parsed = _positive_int(value)
    if parsed > 4:
        raise argparse.ArgumentTypeError("value must be at most 4")
    return parsed


def _positive_float(value: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed) or parsed <= 0:
        raise argparse.ArgumentTypeError("value must be a positive finite number")
    return parsed


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str),
        encoding="utf-8",
    )


def _setting_id(value: object) -> str:
    setting = str(value).strip()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", setting):
        raise ValueError(
            "setting IDs must start with an alphanumeric and contain only "
            "alphanumerics, dot, underscore, or hyphen"
        )
    return setting


def _deduplicated_numbers(
    values: Sequence[Any],
    *,
    name: str,
    integer: bool,
) -> list[int] | list[float]:
    if integer:
        parsed_ints: list[int] = []
        for value in values:
            if isinstance(value, bool):
                raise ValueError(f"{name} values must be integers")
            numeric = int(value)
            if float(value) != numeric:
                raise ValueError(f"{name} values must be integers")
            parsed_ints.append(numeric)
        if any(value <= 0 for value in parsed_ints):
            raise ValueError(f"{name} values must be positive")
        return list(dict.fromkeys(parsed_ints))

    parsed_floats = [float(value) for value in values]
    if any(not math.isfinite(value) or value <= 0 for value in parsed_floats):
        raise ValueError(f"{name} values must be positive and finite")
    return list(dict.fromkeys(parsed_floats))


def _resolved_settings(
    config_path: Path | None,
    *,
    ridge_alphas: Sequence[float] | None,
    lightgbm_seeds: Sequence[int] | None,
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "ridge_alphas": list(_DEFAULT_RIDGE_ALPHAS),
        "lightgbm_seeds": list(_DEFAULT_LIGHTGBM_SEEDS),
        "lightgbm_settings": [
            {"id": row["id"], "params": dict(row["params"])}
            for row in _DEFAULT_LIGHTGBM_SETTINGS
        ],
    }
    if config_path is not None:
        loaded = json.loads(config_path.read_text(encoding="utf-8"))
        if not isinstance(loaded, dict):
            raise ValueError("settings config must contain one JSON object")
        unknown = sorted(
            set(loaded) - {"ridge_alphas", "lightgbm_seeds", "lightgbm_settings"}
        )
        if unknown:
            raise ValueError("unknown settings config keys: " + ", ".join(unknown))
        payload.update(loaded)
    if ridge_alphas is not None:
        payload["ridge_alphas"] = list(ridge_alphas)
    if lightgbm_seeds is not None:
        payload["lightgbm_seeds"] = list(lightgbm_seeds)

    raw_alphas = payload.get("ridge_alphas")
    raw_seeds = payload.get("lightgbm_seeds")
    raw_lightgbm = payload.get("lightgbm_settings")
    if not isinstance(raw_alphas, list):
        raise ValueError("ridge_alphas must be a JSON list")
    if not isinstance(raw_seeds, list):
        raise ValueError("lightgbm_seeds must be a JSON list")
    if not isinstance(raw_lightgbm, list):
        raise ValueError("lightgbm_settings must be a JSON list")
    alphas = _deduplicated_numbers(
        raw_alphas,
        name="ridge_alphas",
        integer=False,
    )
    seeds = _deduplicated_numbers(
        raw_seeds,
        name="lightgbm_seeds",
        integer=True,
    )

    lightgbm_settings: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    for raw in raw_lightgbm:
        if not isinstance(raw, dict) or set(raw) != {"id", "params"}:
            raise ValueError(
                "each lightgbm_settings entry must contain exactly id and params"
            )
        setting = _setting_id(raw["id"])
        if setting in seen_ids:
            raise ValueError(f"duplicate LightGBM setting ID: {setting}")
        params = raw["params"]
        if not isinstance(params, dict):
            raise ValueError(f"LightGBM params for {setting} must be a JSON object")
        seen_ids.add(setting)
        lightgbm_settings.append({"id": setting, "params": dict(params)})
    if lightgbm_settings and not seeds:
        raise ValueError("at least one lightgbm_seeds value is required")
    if not alphas and not lightgbm_settings:
        raise ValueError(
            "settings audit requires at least one Ridge or LightGBM setting"
        )
    return {
        "ridge_alphas": alphas,
        "lightgbm_seeds": seeds,
        "lightgbm_settings": lightgbm_settings,
    }


def _alpha_token(alpha: float) -> str:
    return format(float(alpha), ".12g").replace("-", "m").replace(".", "p")


def _run_definitions(
    settings: Mapping[str, Any], *, ridge_seed: int
) -> list[dict[str, Any]]:
    definitions: list[dict[str, Any]] = []
    for alpha in settings["ridge_alphas"]:
        setting_id = f"ridge_alpha_{_alpha_token(float(alpha))}"
        definitions.append(
            {
                "setting_id": setting_id,
                "run_id": setting_id,
                "model_family": "ridge",
                "model_type": "ridge",
                "model_seed": int(ridge_seed),
                "model_params": {"alpha": float(alpha)},
            }
        )
    ridge_setting_ids = {str(row["setting_id"]) for row in definitions}
    for setting in settings["lightgbm_settings"]:
        if str(setting["id"]) in ridge_setting_ids:
            raise ValueError(
                f"LightGBM setting ID collides with Ridge: {setting['id']}"
            )
        for seed in settings["lightgbm_seeds"]:
            definitions.append(
                {
                    "setting_id": str(setting["id"]),
                    "run_id": f"{setting['id']}_seed_{int(seed)}",
                    "model_family": "lightgbm",
                    "model_type": "lightgbm",
                    "model_seed": int(seed),
                    "model_params": dict(setting["params"]),
                }
            )
    run_ids = [str(row["run_id"]) for row in definitions]
    if len(run_ids) != len(set(run_ids)):
        raise ValueError("resolved settings produce duplicate run IDs")
    return definitions


def _metric_columns(frame: pd.DataFrame, top_k: Sequence[int]) -> list[str]:
    requested = list(_BASE_METRICS)
    for value in top_k:
        requested.extend(
            [
                f"precision_at_{value}",
                f"recall_at_{value}",
                f"enrichment_at_{value}",
                f"positives_recovered_at_{value}",
            ]
        )
    return [column for column in requested if column in frame.columns]


def _reference_setting(summary: pd.DataFrame) -> str:
    preferred = "ridge_alpha_1"
    settings = summary["setting_id"].astype(str).tolist()
    return preferred if preferred in settings else settings[0]


def _same_fold_values(
    fold_metrics: pd.DataFrame,
    *,
    setting_id: str,
    metric: str,
) -> pd.Series:
    selected = fold_metrics.loc[
        fold_metrics["setting_id"].eq(setting_id),
        ["fold", metric],
    ].copy()
    selected[metric] = pd.to_numeric(selected[metric], errors="coerce")
    return selected.groupby("fold", sort=True)[metric].mean()


def _architecture_sensitivity(
    summary: pd.DataFrame,
    fold_metrics: pd.DataFrame,
    *,
    metrics: Sequence[str],
    reference_setting: str,
) -> pd.DataFrame:
    reference_summary = summary.loc[summary["setting_id"].eq(reference_setting)]
    rows: list[dict[str, Any]] = []
    for setting_id, group in summary.groupby("setting_id", sort=False):
        family = str(group["model_family"].iloc[0])
        for metric in metrics:
            values = pd.to_numeric(group[metric], errors="coerce").dropna()
            reference_values = pd.to_numeric(
                reference_summary[metric],
                errors="coerce",
            ).dropna()
            pooled_mean = float(values.mean()) if not values.empty else math.nan
            reference_mean = (
                float(reference_values.mean())
                if not reference_values.empty
                else math.nan
            )
            current_folds = _same_fold_values(
                fold_metrics,
                setting_id=str(setting_id),
                metric=metric,
            )
            reference_folds = _same_fold_values(
                fold_metrics,
                setting_id=reference_setting,
                metric=metric,
            )
            paired = pd.concat(
                [current_folds.rename("current"), reference_folds.rename("reference")],
                axis=1,
                join="inner",
            ).dropna()
            deltas = paired["current"] - paired["reference"]
            rows.append(
                {
                    "setting_id": setting_id,
                    "model_family": family,
                    "n_model_seeds": int(group["model_seed"].nunique()),
                    "model_seeds": json.dumps(
                        sorted(int(value) for value in group["model_seed"].unique())
                    ),
                    "metric": metric,
                    "direction": (
                        "lower_is_better"
                        if metric in _LOWER_IS_BETTER
                        else "higher_is_better_or_descriptive"
                    ),
                    "pooled_mean_across_seeds": pooled_mean,
                    "pooled_std_across_seeds": (
                        float(values.std(ddof=1)) if len(values) > 1 else 0.0
                    ),
                    "reference_setting_id": reference_setting,
                    "reference_pooled_mean": reference_mean,
                    "pooled_delta_from_reference": pooled_mean - reference_mean,
                    "n_paired_folds": int(len(deltas)),
                    "mean_same_fold_delta": (
                        float(deltas.mean()) if not deltas.empty else math.nan
                    ),
                    "median_same_fold_delta": (
                        float(deltas.median()) if not deltas.empty else math.nan
                    ),
                    "std_same_fold_delta": (
                        float(deltas.std(ddof=1)) if len(deltas) > 1 else 0.0
                    ),
                }
            )
    return pd.DataFrame(rows)


def _stochastic_variation(
    summary: pd.DataFrame,
    fold_metrics: pd.DataFrame,
    *,
    metrics: Sequence[str],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for setting_id, group in summary.groupby("setting_id", sort=False):
        if group["model_seed"].nunique() < 2:
            continue
        family = str(group["model_family"].iloc[0])
        fold_group = fold_metrics.loc[fold_metrics["setting_id"].eq(setting_id)]
        for metric in metrics:
            pooled = pd.to_numeric(group[metric], errors="coerce").dropna()
            per_fold = fold_group[["fold", "model_seed", metric]].copy()
            per_fold[metric] = pd.to_numeric(per_fold[metric], errors="coerce")
            by_fold = per_fold.groupby("fold", sort=True)[metric]
            fold_ranges = by_fold.max() - by_fold.min()
            fold_stds = by_fold.std(ddof=1)
            rows.append(
                {
                    "setting_id": setting_id,
                    "model_family": family,
                    "n_model_seeds": int(group["model_seed"].nunique()),
                    "model_seeds": json.dumps(
                        sorted(int(value) for value in group["model_seed"].unique())
                    ),
                    "metric": metric,
                    "pooled_seed_mean": (
                        float(pooled.mean()) if not pooled.empty else math.nan
                    ),
                    "pooled_seed_std": (
                        float(pooled.std(ddof=1)) if len(pooled) > 1 else 0.0
                    ),
                    "pooled_seed_min": (
                        float(pooled.min()) if not pooled.empty else math.nan
                    ),
                    "pooled_seed_max": (
                        float(pooled.max()) if not pooled.empty else math.nan
                    ),
                    "pooled_seed_range": (
                        float(pooled.max() - pooled.min())
                        if not pooled.empty
                        else math.nan
                    ),
                    "mean_same_fold_seed_std": (
                        float(fold_stds.mean()) if fold_stds.notna().any() else math.nan
                    ),
                    "mean_same_fold_seed_range": (
                        float(fold_ranges.mean())
                        if fold_ranges.notna().any()
                        else math.nan
                    ),
                    "max_same_fold_seed_range": (
                        float(fold_ranges.max())
                        if fold_ranges.notna().any()
                        else math.nan
                    ),
                }
            )
    return pd.DataFrame(rows)


def run_spd_exposure_model_settings_audit(
    dataset_path: str | Path,
    out_dir: str | Path,
    *,
    potency_features: Sequence[str],
    pk_features: Sequence[str],
    settings: Mapping[str, Any],
    label_col: str = "spd_exposure_label",
    group_col: str = "drug_id",
    n_splits: int = 5,
    fold_seed: int = 42,
    max_threads: int = 4,
    censored_policy: str = "exclude",
    top_k: Sequence[int] = (10, 20, 100),
) -> dict[str, Any]:
    """Run a compact, same-fold parameter/seed sensitivity audit sequentially."""
    if not 1 <= int(max_threads) <= 4:
        raise ValueError("max_threads must be between 1 and 4")
    top_k_values = list(dict.fromkeys(int(value) for value in top_k))
    if not top_k_values or any(value <= 0 for value in top_k_values):
        raise ValueError("top_k values must be positive integers")
    definitions = _run_definitions(settings, ridge_seed=fold_seed)
    dataset = Path(dataset_path)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    shared_folds = out / "shared_fold_assignments.csv"
    fold_manifest = save_spd_exposure_group_folds(
        dataset,
        shared_folds,
        label_col=label_col,
        group_col=group_col,
        n_splits=n_splits,
        seed=fold_seed,
    )
    effective_splits = int(fold_manifest["effective_n_splits"])
    fold_signature = str(fold_manifest["fold_signature_sha256"])
    resolved_config = {
        "schema_version": 1,
        "dataset": str(dataset),
        "label_col": label_col,
        "group_col": group_col,
        "potency_features": list(potency_features),
        "pk_features": list(pk_features),
        "fold_seed": int(fold_seed),
        "requested_n_splits": int(n_splits),
        "effective_n_splits": effective_splits,
        "fold_signature_sha256": fold_signature,
        "max_threads": int(max_threads),
        "censored_policy": censored_policy,
        "top_k": top_k_values,
        "permutation_repeats": 0,
        "settings": dict(settings),
        "resolved_runs": definitions,
    }
    resolved_config_path = out / "resolved_settings_config.json"
    _write_json(resolved_config_path, resolved_config)

    summary_rows: list[dict[str, Any]] = []
    fold_frames: list[pd.DataFrame] = []
    ledger_frames: list[pd.DataFrame] = []
    model_setting_frames: list[pd.DataFrame] = []
    run_manifests: list[str] = []
    for definition in definitions:
        run_dir = out / "runs" / str(definition["run_id"])
        manifest = run_spd_exposure_grouped_oof(
            dataset,
            run_dir,
            label_col=label_col,
            potency_features=potency_features,
            pk_features=pk_features,
            group_col=group_col,
            n_splits=effective_splits,
            model_type=str(definition["model_type"]),
            model_params=definition["model_params"],
            seed=int(definition["model_seed"]),
            max_threads=int(max_threads),
            fold_assignments_path=shared_folds,
            fold_seed=fold_seed,
            censored_policy=censored_policy,
            permutation_repeats=0,
            top_k=top_k_values,
        )
        observed_signature = str(manifest["split"]["fold_signature_sha256"])
        if observed_signature != fold_signature:
            raise AssertionError(
                f"run {definition['run_id']} did not reuse the shared folds"
            )
        pooled = dict(manifest["pooled_metrics"])
        exact_settings = dict(manifest["model_configuration"])
        summary_rows.append(
            {
                "setting_id": definition["setting_id"],
                "run_id": definition["run_id"],
                "model_family": definition["model_family"],
                "model_type": definition["model_type"],
                "model_seed": int(definition["model_seed"]),
                "fold_signature_sha256": fold_signature,
                "model_params_json": json.dumps(
                    definition["model_params"],
                    sort_keys=True,
                    default=str,
                ),
                "exact_settings_json": json.dumps(
                    exact_settings,
                    sort_keys=True,
                    default=str,
                ),
                "run_dir": str(run_dir),
                **pooled,
            }
        )
        fold_frame = pd.read_csv(run_dir / "per_fold_metrics.csv")
        fold_frame["setting_id"] = definition["setting_id"]
        fold_frame["run_id"] = definition["run_id"]
        fold_frame["model_family"] = definition["model_family"]
        fold_frame["model_seed"] = int(definition["model_seed"])
        fold_frames.append(fold_frame)
        ledger_frames.append(pd.read_csv(run_dir / "model_run_record.csv"))
        model_settings = pd.read_csv(run_dir / "model_settings.csv")
        model_settings["setting_id"] = definition["setting_id"]
        model_settings["run_id"] = definition["run_id"]
        model_settings["model_family"] = definition["model_family"]
        model_setting_frames.append(model_settings)
        run_manifests.append(str(manifest["outputs"]["manifest"]))

    summary = pd.DataFrame(summary_rows)
    fold_metrics = pd.concat(fold_frames, ignore_index=True)
    ledger = pd.concat(ledger_frames, ignore_index=True).reindex(columns=LEDGER_COLUMNS)
    resolved_model_settings = pd.concat(model_setting_frames, ignore_index=True)
    metrics = _metric_columns(summary, top_k_values)
    reference_setting = _reference_setting(summary)
    architecture = _architecture_sensitivity(
        summary,
        fold_metrics,
        metrics=metrics,
        reference_setting=reference_setting,
    )
    stochastic = _stochastic_variation(
        summary,
        fold_metrics,
        metrics=metrics,
    )

    outputs = {
        "resolved_settings_config": resolved_config_path,
        "shared_fold_assignments": shared_folds,
        "shared_fold_manifest": Path(str(fold_manifest["manifest_path"])),
        "settings_summary": out / "settings_summary.csv",
        "settings_fold_metrics": out / "settings_fold_metrics.csv",
        "resolved_model_settings": out / "resolved_model_settings.csv",
        "architecture_sensitivity": out / "architecture_sensitivity.csv",
        "stochastic_variation": out / "stochastic_variation.csv",
        "model_ledger": out / "ml_model_run_ledger.csv",
        "manifest": out / "spd_exposure_model_settings_audit_manifest.json",
    }
    summary.to_csv(outputs["settings_summary"], index=False)
    fold_metrics.to_csv(outputs["settings_fold_metrics"], index=False)
    resolved_model_settings.to_csv(outputs["resolved_model_settings"], index=False)
    architecture.to_csv(outputs["architecture_sensitivity"], index=False)
    stochastic.to_csv(outputs["stochastic_variation"], index=False)
    ledger.to_csv(outputs["model_ledger"], index=False)

    audit_manifest: dict[str, Any] = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "audit": "spd_exposure_same_fold_model_settings",
        "scope": "compact_settings_sensitivity_only",
        "full_ml_audit_run": False,
        "dataset": str(dataset),
        "out_dir": str(out),
        "fold_signature_sha256": fold_signature,
        "fold_manifest": fold_manifest,
        "reference_setting_id": reference_setting,
        "n_settings": int(summary["setting_id"].nunique()),
        "n_runs": int(len(summary)),
        "max_threads": int(max_threads),
        "execution": "sequential",
        "metrics": metrics,
        "settings": dict(settings),
        "run_manifests": run_manifests,
        "interpretation": {
            "architecture_sensitivity": (
                "Compare setting means and same-fold deltas after averaging repeated "
                "LightGBM seeds within each setting."
            ),
            "stochastic_variation": (
                "Compare seed standard deviations and ranges only within the same "
                "LightGBM setting; folds, preprocessing contract, and hyperparameters "
                "remain fixed."
            ),
        },
        "outputs": {key: str(value) for key, value in outputs.items()},
    }
    _write_json(outputs["manifest"], audit_manifest)
    return audit_manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compare Ridge alphas and a small regularized LightGBM grid/seeds for "
            "the two-stage SPD exposure model on one saved grouped-OOF split."
        )
    )
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--label", default="spd_exposure_label")
    parser.add_argument("--potency-feature", action="append", required=True)
    parser.add_argument("--pk-feature", action="append", required=True)
    parser.add_argument("--group-col", default="drug_id")
    parser.add_argument("--n-splits", type=_positive_int, default=5)
    parser.add_argument("--fold-seed", type=int, default=42)
    parser.add_argument(
        "--settings-config",
        type=Path,
        default=None,
        help=(
            "Optional JSON object with ridge_alphas, lightgbm_seeds, and "
            "lightgbm_settings; resolved settings are always written to the output."
        ),
    )
    parser.add_argument(
        "--ridge-alpha",
        action="append",
        type=_positive_float,
        default=None,
        help="Ridge alpha; repeat to override the default 0.1, 1, 10 grid.",
    )
    parser.add_argument(
        "--lightgbm-seed",
        action="append",
        type=int,
        default=None,
        help="LightGBM seed; repeat to override the default 42, 137, 314 seeds.",
    )
    parser.add_argument("--max-threads", type=_thread_count, default=4)
    parser.add_argument(
        "--censored-policy",
        choices=["exclude", "bound"],
        default="exclude",
    )
    parser.add_argument(
        "--top-k",
        action="append",
        type=_positive_int,
        default=None,
        help="Top-K cutoff; repeat for multiple cutoffs (default: 10, 20, 100).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    settings = _resolved_settings(
        args.settings_config,
        ridge_alphas=args.ridge_alpha,
        lightgbm_seeds=args.lightgbm_seed,
    )
    manifest = run_spd_exposure_model_settings_audit(
        args.dataset,
        args.out_dir,
        potency_features=args.potency_feature,
        pk_features=args.pk_feature,
        settings=settings,
        label_col=args.label,
        group_col=args.group_col,
        n_splits=args.n_splits,
        fold_seed=args.fold_seed,
        max_threads=args.max_threads,
        censored_policy=args.censored_policy,
        top_k=args.top_k or (10, 20, 100),
    )
    print(json.dumps(manifest, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

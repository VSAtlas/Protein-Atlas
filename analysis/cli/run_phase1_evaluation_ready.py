from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.feature_metadata import refresh_tables


DEFAULT_DATASET = Path(
    "data/spd90_fda_fda_dud_spr10_40h_patchdist_20260525_161503/"
    "ml_training_pass_phase1_v1_20260706/merged_spd_external_four_state/"
    "spd90_model_ready_plus_external_four_state.csv"
)


EVALUATION_SPLITS = [
    "drug_holdout",
    "chemical_cluster_holdout",
    "target_holdout",
    "target_family_holdout",
]


CORE_MODELS = ["logistic_regression", "elastic_net", "random_forest", "lightgbm", "catboost", "ebm"]


def _run(command: list[str], *, log_path: Path, timeout: int | None = None) -> dict[str, Any]:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc)
    env = os.environ.copy()
    thread_limit = str(max(1, min(int(env.get("ATLAS_EVAL_THREADS", "8")), 32)))
    for name in [
        "OPENBLAS_NUM_THREADS",
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "BLIS_NUM_THREADS",
    ]:
        env[name] = thread_limit
    with log_path.open("w", encoding="utf-8") as handle:
        handle.write(" ".join(command) + "\n\n")
        handle.write(f"[thread_limit] {thread_limit}\n\n")
        try:
            completed = subprocess.run(
                command,
                check=False,
                env=env,
                stdout=handle,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=timeout,
            )
            returncode = int(completed.returncode)
            timed_out = False
        except subprocess.TimeoutExpired as exc:
            handle.write(f"\n[TIMEOUT] command exceeded timeout={timeout}: {exc}\n")
            returncode = 124
            timed_out = True
    finished = datetime.now(timezone.utc)
    return {
        "command": command,
        "log": str(log_path),
        "returncode": returncode,
        "timed_out": timed_out,
        "started_at": started.isoformat(),
        "finished_at": finished.isoformat(),
        "elapsed_seconds": (finished - started).total_seconds(),
    }


def _copy_csv(src: Path, dst: Path) -> bool:
    if not src.exists():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    pd.read_csv(src, low_memory=False).to_csv(dst, index=False)
    return True


def _combine_summaries(paths: list[tuple[str, Path]], out_path: Path) -> int:
    frames: list[pd.DataFrame] = []
    for split, path in paths:
        if not path.exists():
            continue
        frame = pd.read_csv(path, low_memory=False)
        if "split_method" not in frame.columns:
            frame.insert(0, "split_method", split)
        else:
            frame["split_method"] = frame["split_method"].fillna(split)
        frames.append(frame)
    if not frames:
        return 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.concat(frames, ignore_index=True, sort=False).to_csv(out_path, index=False)
    return int(sum(len(frame) for frame in frames))


def _write_label_status(dataset: Path, out_dir: Path) -> dict[str, Any]:
    frame = pd.read_csv(dataset, low_memory=False)
    rows: list[dict[str, Any]] = []
    for label in [
        "spd_binding_label",
        "spd_exposure_relevant",
        "tissue_site_label",
        "external_four_state_ml_label",
        "combined_activity_ml_label",
    ]:
        if label not in frame.columns:
            rows.append({"label": label, "status": "missing"})
            continue
        counts = frame[label].value_counts(dropna=False).to_dict()
        rows.append(
            {
                "label": label,
                "status": "present",
                "n_rows": int(len(frame)),
                "n_positive": int((pd.to_numeric(frame[label], errors="coerce") == 1).sum()),
                "n_negative": int((pd.to_numeric(frame[label], errors="coerce") == 0).sum()),
                "n_unknown": int(pd.to_numeric(frame[label], errors="coerce").isna().sum()),
                "raw_counts_json": json.dumps({str(k): int(v) for k, v in counts.items()}, sort_keys=True),
            }
        )
    status_path = out_dir / "evaluation_ready_label_status.csv"
    pd.DataFrame(rows).to_csv(status_path, index=False)
    return {"path": str(status_path), "rows": rows}


def _best_rows(summary_path: Path, out_path: Path) -> int:
    if not summary_path.exists():
        return 0
    frame = pd.read_csv(summary_path, low_memory=False)
    if frame.empty:
        return 0
    score_col = "PR_AUC" if "PR_AUC" in frame.columns else "AUPRC"
    groups = [col for col in ["split_method", "expert"] if col in frame.columns]
    if not groups:
        groups = ["split_method"] if "split_method" in frame.columns else []
    rows = []
    if groups:
        for _, sub in frame[frame.get("status", "ok").fillna("ok").isin(["ok", "trained", "skipped_existing"])].groupby(groups, dropna=False):
            usable = sub.dropna(subset=[score_col]) if score_col in sub else sub
            if usable.empty:
                usable = sub.dropna(subset=["enrichment_at_K"]) if "enrichment_at_K" in sub else sub
            if not usable.empty:
                key = score_col if score_col in usable and usable[score_col].notna().any() else "enrichment_at_K"
                rows.append(usable.sort_values(key, ascending=False).iloc[0])
    else:
        rows = [frame.iloc[0]]
    if not rows:
        return 0
    pd.DataFrame(rows).to_csv(out_path, index=False)
    return len(rows)



def _summary_has_split(summary_path: Path, split: str) -> bool:
    if not summary_path.exists():
        return False
    frame = pd.read_csv(summary_path, low_memory=False)
    return "split_method" in frame.columns and frame["split_method"].astype(str).eq(split).any()


def _summary_has_column(summary_path: Path, column_candidates: list[str]) -> bool:
    if not summary_path.exists():
        return False
    frame = pd.read_csv(summary_path, low_memory=False, nrows=5)
    return any(column in frame.columns for column in column_candidates)


def _glob_present(root: Path, pattern: str) -> bool:
    return any(root.glob(pattern))


def _existing_output_path(configured: str | None, *fallbacks: Path) -> Path:
    if configured:
        path = Path(configured)
        if path.exists():
            return path
    for path in fallbacks:
        if path.exists():
            return path
    return fallbacks[0] if fallbacks else Path("")



def _safe_metrics(label: pd.Series, score: pd.Series, *, top_k: int) -> dict[str, Any]:
    from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

    y = pd.to_numeric(label, errors="coerce")
    s = pd.to_numeric(score, errors="coerce")
    mask = y.isin([0, 1]) & s.notna()
    y = y[mask].astype(int)
    s = s[mask].astype(float)
    n = int(len(y))
    positives = int(y.sum())
    negatives = int(n - positives)
    if n == 0:
        return {
            "n_test": 0,
            "n_positive": 0,
            "n_negative": 0,
            "AUROC": None,
            "PR_AUC": None,
            "precision_at_K": None,
            "enrichment_at_K": None,
            "Brier_score": None,
            "prevalence": None,
        }
    auroc = float(roc_auc_score(y, s)) if positives > 0 and negatives > 0 else None
    pr_auc = float(average_precision_score(y, s)) if positives > 0 else None
    k = min(int(top_k), n)
    order = s.sort_values(ascending=False).index[:k]
    precision = float(y.loc[order].mean()) if k > 0 else None
    prevalence = float(positives / n) if n else None
    enrichment = float(precision / prevalence) if precision is not None and prevalence else None
    brier = None
    if float(s.min()) >= 0.0 and float(s.max()) <= 1.0:
        brier = float(brier_score_loss(y, s.clip(0, 1)))
    return {
        "n_test": n,
        "n_positive": positives,
        "n_negative": negatives,
        "AUROC": auroc,
        "PR_AUC": pr_auc,
        "precision_at_K": precision,
        "enrichment_at_K": enrichment,
        "Brier_score": brier,
        "prevalence": prevalence,
    }


def _write_exposure_prior_baseline(dataset: Path, out_dir: Path, splits: list[str], *, seed: int, top_k: int) -> str | None:
    from analysis.ml.splits import make_split

    frame = pd.read_csv(dataset, low_memory=False)
    label = "spd_exposure_relevant"
    score_candidates = [
        "binding_expert_score",
        "banana_atlas_blend_score",
        "banana_binding_probability",
        "consensus_score_normalized",
        "atlas_score_normalized",
        "banana_score_normalized",
        "consensus_score",
    ]
    if label not in frame.columns:
        return None
    score_col = next(
        (col for col in score_candidates if col in frame.columns and pd.to_numeric(frame[col], errors="coerce").notna().any()),
        None,
    )
    if score_col is None:
        return None
    work = frame[pd.to_numeric(frame[label], errors="coerce").isin([0, 1])].copy()
    rows: list[dict[str, Any]] = []
    for split in splits:
        try:
            _, test_idx = make_split(work, split_mode=split, seed=seed)
            test = work.loc[test_idx]
            status = "ok"
            error = ""
        except Exception as exc:
            test = work
            status = "failed_split"
            error = str(exc)
        metric_row = _safe_metrics(test[label], test[score_col], top_k=top_k)
        rows.append(
            {
                "expert": "exposure",
                "model_type": "binding_prior_only",
                "status": status,
                "label_used": label,
                "feature_set": "frozen_binding_prior_only",
                "split_method": split,
                "PU_strategy": "none",
                "calibration_method": "none",
                "score_col": score_col,
                "top_k": int(top_k),
                "interpretation": (
                    "Exposure baseline that ranks SPD exposure relevance using only the frozen binding prior. "
                    "It excludes AC50, free Cmax, Cmax, exposure margin, and observed SPD assay fields."
                ),
                "error": error,
                **metric_row,
            }
        )
    out_path = out_dir / "evaluation_ready_exposure_binding_prior_baseline.csv"
    pd.DataFrame(rows).to_csv(out_path, index=False)
    return str(out_path)


def _write_metric_coverage(out: Path, *, outputs: dict[str, Any]) -> str:
    ready_outputs = outputs.get("evaluation_ready", {})
    binding_ablation = _existing_output_path(
        ready_outputs.get("binding_ablation_comparison"),
        out / "evaluation_ready_binding_ablation_comparison.csv",
        out / "binding_ablation" / "binding_ablation_comparison.csv",
    )
    binding_family = _existing_output_path(
        ready_outputs.get("binding_model_family_summary"),
        out / "evaluation_ready_binding_model_family_summary.csv",
    )
    exposure_family = _existing_output_path(
        ready_outputs.get("exposure_model_family_summary"),
        out / "evaluation_ready_exposure_model_family_summary.csv",
    )
    grouped_binding = _existing_output_path(
        ready_outputs.get("binding_grouped_cv_stability_summary.csv"),
        out / "evaluation_ready_binding_grouped_cv_stability_summary.csv",
    )
    grouped_exposure = _existing_output_path(
        ready_outputs.get("exposure_grouped_cv_stability_summary.csv"),
        out / "evaluation_ready_exposure_grouped_cv_stability_summary.csv",
    )

    def both_have(columns: list[str]) -> bool:
        binding_ok = _summary_has_column(binding_ablation, columns) or _summary_has_column(binding_family, columns)
        exposure_ok = _summary_has_column(exposure_family, columns)
        return binding_ok and exposure_ok

    rows = [
        {
            "evaluation": "PR-AUC / average precision",
            "status": "present" if both_have(["PR_AUC", "AUPRC"]) else "partial_or_missing",
            "ready_for_phase1_claim": bool(both_have(["PR_AUC", "AUPRC"])),
            "evidence_file": ";".join(str(p) for p in [binding_ablation, binding_family, exposure_family] if p.exists()),
            "notes": "Primary imbalanced-class ranking metric for SPD binding and exposure.",
        },
        {
            "evaluation": "Precision@K",
            "status": "present" if both_have(["precision_at_K"]) else "partial_or_missing",
            "ready_for_phase1_claim": bool(both_have(["precision_at_K"])),
            "evidence_file": ";".join(str(p) for p in [binding_ablation, binding_family, exposure_family] if p.exists()),
            "notes": "Decision-focused top-ranked candidate precision.",
        },
        {
            "evaluation": "Enrichment@K",
            "status": "present" if both_have(["enrichment_at_K", "EF@1%", "EF@5%", "EF@10%"]) else "partial_or_missing",
            "ready_for_phase1_claim": bool(both_have(["enrichment_at_K", "EF@1%", "EF@5%", "EF@10%"])),
            "evidence_file": ";".join(str(p) for p in [binding_ablation, binding_family, exposure_family] if p.exists()),
            "notes": "Shows top-ranked enrichment over label prevalence.",
        },
        {
            "evaluation": "Calibration curve / reliability table",
            "status": "present"
            if _glob_present(out, "binding_model_family/**/ml_calibration_curve.png")
            and _glob_present(out, "exposure_model_family/**/ml_calibration_curve.png")
            else "partial_or_missing",
            "ready_for_phase1_claim": bool(
                _glob_present(out, "binding_model_family/**/model_reliability_table.csv")
                and _glob_present(out, "exposure_model_family/**/model_reliability_table.csv")
            ),
            "evidence_file": "binding_model_family/**/ml_calibration_curve.png;exposure_model_family/**/ml_calibration_curve.png",
            "notes": "Binding ablation is metrics-only; calibrated probabilities come from model-family runs.",
        },
        {
            "evaluation": "Brier score",
            "status": "present" if both_have(["Brier_score", "Brier"]) else "partial_or_missing",
            "ready_for_phase1_claim": bool(both_have(["Brier_score", "Brier"])),
            "evidence_file": ";".join(str(p) for p in [binding_ablation, binding_family, exposure_family] if p.exists()),
            "notes": "Probability-quality metric; only meaningful for probability-like scores.",
        },
        {
            "evaluation": "Ablation tests",
            "status": "present" if binding_ablation.exists() else "missing",
            "ready_for_phase1_claim": bool(binding_ablation.exists()),
            "evidence_file": str(binding_ablation) if binding_ablation.exists() else "",
            "notes": "BANANA-only, consensus-only, RDKit-only, Atlas/BANANA, and full clean binding comparisons.",
        },
        {
            "evaluation": "Held-out scaffold / chemotype split",
            "status": "present"
            if any(_summary_has_split(path, "chemical_cluster_holdout") for path in [binding_ablation, binding_family, exposure_family])
            else "missing",
            "ready_for_phase1_claim": bool(
                _summary_has_split(exposure_family, "chemical_cluster_holdout")
                and (_summary_has_split(binding_ablation, "chemical_cluster_holdout") or _summary_has_split(binding_family, "chemical_cluster_holdout"))
            ),
            "evidence_file": ";".join(str(p) for p in [binding_ablation, binding_family, exposure_family] if p.exists()),
            "notes": "Uses chemical_cluster/scaffold fallback for chemistry generalization.",
        },
        {
            "evaluation": "Held-out drug split",
            "status": "present"
            if any(_summary_has_split(path, "drug_holdout") for path in [binding_ablation, binding_family, exposure_family])
            else "missing",
            "ready_for_phase1_claim": bool(
                _summary_has_split(exposure_family, "drug_holdout")
                and (_summary_has_split(binding_ablation, "drug_holdout") or _summary_has_split(binding_family, "drug_holdout"))
            ),
            "evidence_file": ";".join(str(p) for p in [binding_ablation, binding_family, exposure_family] if p.exists()),
            "notes": "Tests new-drug generalization.",
        },
        {
            "evaluation": "Held-out target/PDB split",
            "status": "present"
            if any(_summary_has_split(path, "target_holdout") for path in [binding_ablation, binding_family, exposure_family])
            else "missing",
            "ready_for_phase1_claim": bool(
                _summary_has_split(exposure_family, "target_holdout")
                and (_summary_has_split(binding_ablation, "target_holdout") or _summary_has_split(binding_family, "target_holdout"))
            ),
            "evidence_file": ";".join(str(p) for p in [binding_ablation, binding_family, exposure_family] if p.exists()),
            "notes": "Tests new-target/PDB generalization.",
        },
        {
            "evaluation": "Family-only diagnostic",
            "status": "present"
            if any(_summary_has_split(path, "target_family_holdout") for path in [binding_ablation, binding_family, exposure_family])
            else "missing",
            "ready_for_phase1_claim": False,
            "evidence_file": ";".join(str(p) for p in [binding_ablation, binding_family, exposure_family] if p.exists()),
            "notes": "Diagnostic for target-family memorization; not a primary Phase 1 claim if positive/negative families are imbalanced.",
        },
        {
            "evaluation": "Repeated grouped CV / confidence intervals",
            "status": "present" if grouped_binding.exists() or grouped_exposure.exists() else "missing",
            "ready_for_phase1_claim": bool(grouped_binding.exists() and grouped_exposure.exists()),
            "evidence_file": ";".join(str(p) for p in [grouped_binding, grouped_exposure] if p.exists()),
            "notes": "Pooled out-of-fold predictions, bootstrap CIs, and permutation checks when grouped CV completes.",
        },
        {
            "evaluation": "External validation",
            "status": "deferred_to_source_transfer_suite",
            "ready_for_phase1_claim": False,
            "evidence_file": "python -m analysis.cli.run_ml_source_pu_suite --out-dir <run_dir>/ml_source_pu_suite",
            "notes": "ChEMBL/BindingDB/ToxCast labels are separate activity objectives; report source-transfer as a stress test, not as direct SPD exposure validation.",
        },
    ]
    path = out / "evaluation_ready_metric_coverage.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run Atlas Phase 1 evaluation-ready validation outputs for binding and direct SPD exposure."
    )
    parser.add_argument("--dataset", type=Path, default=DEFAULT_DATASET)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--models", nargs="+", default=CORE_MODELS)
    parser.add_argument("--splits", nargs="+", default=EVALUATION_SPLITS)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--top-k", type=int, default=20)
    parser.add_argument("--bootstraps", type=int, default=50)
    parser.add_argument("--grouped-cv-repeats", type=int, default=10)
    parser.add_argument("--grouped-cv-bootstraps", type=int, default=100)
    parser.add_argument("--grouped-cv-permutations", type=int, default=200)
    parser.add_argument("--model-n-jobs", type=int, default=8)
    parser.add_argument("--per-run-timeout", type=int, default=1800)
    parser.add_argument("--skip-binding-ablation", action="store_true")
    parser.add_argument("--skip-binding-model-family", action="store_true")
    parser.add_argument("--skip-exposure-splits", action="store_true")
    parser.add_argument(
        "--run-grouped-cv",
        action="store_true",
        help="Run expensive grouped/repeated CV with bootstrap and permutation CIs. Intended for final/offline validation only.",
    )
    parser.add_argument("--skip-grouped-cv", action="store_true", help="Deprecated compatibility flag; grouped CV is skipped unless --run-grouped-cv is set.")
    parser.add_argument("--skip-existing-models", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--grouped-cv-stratified-target-holdout-by-family", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--grouped-cv-stratified-target-holdout-repeats", type=int, default=None)
    parser.add_argument("--grouped-cv-stratified-target-holdout-fraction", type=float, default=0.2)
    parser.add_argument("--metrics-only-binding", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--refresh-feature-metadata", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--chemical-cluster", choices=["auto", "scaffold", "smiles", "chemotype", "butina", "ecfp", "none"], default="auto")
    args = parser.parse_args(argv)
    if not args.run_grouped_cv:
        args.skip_grouped_cv = True

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    dataset = args.dataset
    feature_metadata_manifest: dict[str, Any] | None = None
    if args.refresh_feature_metadata:
        feature_metadata_manifest = refresh_tables(
            [args.dataset],
            out_dir=out / "_feature_metadata",
            run_dir=None,
            chemical_cluster=args.chemical_cluster,
            target_family="auto",
            source_lineage="auto",
        )
        if feature_metadata_manifest.get("outputs"):
            dataset = Path(feature_metadata_manifest["outputs"][0])
    commands: list[dict[str, Any]] = []
    outputs: dict[str, Any] = {
        "evaluation_ready": {},
        "exploratory": {},
    }

    label_status = _write_label_status(dataset, out)
    outputs["evaluation_ready"]["label_status"] = label_status["path"]
    exposure_prior = _write_exposure_prior_baseline(dataset, out, args.splits, seed=args.seed, top_k=args.top_k)
    if exposure_prior:
        outputs["evaluation_ready"]["exposure_binding_prior_baseline"] = exposure_prior

    if not args.skip_binding_ablation:
        binding_out = out / "binding_ablation"
        command = [
            sys.executable,
            "-m",
            "analysis.cli.run_binding_ablation_comparison",
            "--dataset",
            str(dataset),
            "--label",
            "spd_binding_label",
            "--out-dir",
            str(binding_out),
            "--splits",
            *args.splits,
            "--models",
            *args.models,
            "--seed",
            str(args.seed),
            "--bootstraps",
            str(args.bootstraps),
            "--top-k",
            str(args.top_k),
            "--calibration",
            "sigmoid",
            "--model-n-jobs",
            str(args.model_n_jobs),
            "--claim-mode",
            "publication",
        ]
        if args.metrics_only_binding:
            command.append("--metrics-only")
        commands.append(_run(command, log_path=out / "logs" / "binding_ablation.log", timeout=args.per_run_timeout))
        src = binding_out / "binding_ablation_comparison.csv"
        dst = out / "evaluation_ready_binding_ablation_comparison.csv"
        if _copy_csv(src, dst):
            outputs["evaluation_ready"]["binding_ablation_comparison"] = str(dst)
            best = out / "evaluation_ready_binding_best_by_split.csv"
            _best_rows(dst, best)
            outputs["evaluation_ready"]["binding_best_by_split"] = str(best)

    binding_summaries: list[tuple[str, Path]] = []
    if not args.skip_binding_model_family:
        for split in args.splits:
            split_out = out / "binding_model_family" / split
            command = [
                sys.executable,
                "-m",
                "analysis.cli.run_ml_model_family_sweep",
                "--binding-dataset",
                str(dataset),
                "--experts",
                "binding",
                "--models",
                *args.models,
                "--split",
                split,
                "--calibration",
                "sigmoid",
                "--bootstraps",
                str(args.bootstraps),
                "--model-n-jobs",
                str(args.model_n_jobs),
                "--per-run-timeout",
                str(args.per_run_timeout),
                "--group-reweight-cols",
                "source_family",
                "target_family",
                "--group-reweight-max-factor",
                "3",
                "--out-dir",
                str(split_out),
            ]
            if args.skip_existing_models:
                command.append("--skip-existing")
            commands.append(_run(command, log_path=out / "logs" / f"binding_{split}.log", timeout=args.per_run_timeout * max(1, len(args.models))))
            binding_summaries.append((split, split_out / "model_family_sweep_summary.csv"))
        combined = out / "evaluation_ready_binding_model_family_summary.csv"
        _combine_summaries(binding_summaries, combined)
        outputs["evaluation_ready"]["binding_model_family_summary"] = str(combined)
        best = out / "evaluation_ready_binding_model_family_best_by_split.csv"
        _best_rows(combined, best)
        outputs["evaluation_ready"]["binding_model_family_best_by_split"] = str(best)

    exposure_summaries: list[tuple[str, Path]] = []
    if not args.skip_exposure_splits:
        for split in args.splits:
            split_out = out / "exposure_model_family" / split
            command = [
                sys.executable,
                "-m",
                "analysis.cli.run_ml_model_family_sweep",
                "--exposure-dataset",
                str(dataset),
                "--experts",
                "exposure",
                "--models",
                *args.models,
                "--split",
                split,
                "--calibration",
                "sigmoid",
                "--bootstraps",
                str(args.bootstraps),
                "--model-n-jobs",
                str(args.model_n_jobs),
                "--per-run-timeout",
                str(args.per_run_timeout),
                "--group-reweight-cols",
                "source_family",
                "target_family",
                "--group-reweight-max-factor",
                "3",
                "--out-dir",
                str(split_out),
            ]
            if args.skip_existing_models:
                command.append("--skip-existing")
            commands.append(_run(command, log_path=out / "logs" / f"exposure_{split}.log", timeout=args.per_run_timeout * max(1, len(args.models))))
            exposure_summaries.append((split, split_out / "model_family_sweep_summary.csv"))
        combined = out / "evaluation_ready_exposure_model_family_summary.csv"
        _combine_summaries(exposure_summaries, combined)
        outputs["evaluation_ready"]["exposure_model_family_summary"] = str(combined)
        best = out / "evaluation_ready_exposure_best_by_split.csv"
        _best_rows(combined, best)
        outputs["evaluation_ready"]["exposure_best_by_split"] = str(best)

    if not args.skip_grouped_cv:
        grouped_specs = [
            ("binding", "spd_binding_label", "spd_binding_nonleaky"),
            ("exposure", "spd_exposure_relevant", "spd_exposure_nonleaky"),
        ]
        for name, label, feature_set in grouped_specs:
            cv_out = out / f"{name}_grouped_cv"
            command = [
                sys.executable,
                "-m",
                "analysis.cli.run_grouped_cv_stability",
                "--dataset",
                str(dataset),
                "--label",
                label,
                "--feature-set",
                feature_set,
                "--out-dir",
                str(cv_out),
                "--models",
                *args.models,
                "--group-cols",
                "target_id",
                "chemical_cluster",
                "target_family",
                "--leave-one-group-cols",
                "target_family",
                "--n-splits",
                "5",
                "--repeats",
                str(args.grouped_cv_repeats),
                "--min-test-positives",
                "10",
                "--min-test-negatives",
                "10",
                "--top-k",
                str(args.top_k),
                "--bootstraps",
                str(args.grouped_cv_bootstraps),
                "--permutations",
                str(args.grouped_cv_permutations),
                "--seed",
                str(args.seed),
            ]
            if args.grouped_cv_stratified_target_holdout_by_family:
                command.append("--stratified-target-holdout-by-family")
                if args.grouped_cv_stratified_target_holdout_repeats is not None:
                    command.extend(["--stratified-target-holdout-repeats", str(args.grouped_cv_stratified_target_holdout_repeats)])
                command.extend([
                    "--stratified-target-holdout-fraction",
                    str(args.grouped_cv_stratified_target_holdout_fraction),
                ])
            commands.append(_run(command, log_path=out / "logs" / f"{name}_grouped_cv.log", timeout=args.per_run_timeout * 2))
            for filename in [
                "grouped_cv_stability_summary.csv",
                "pooled_oof_metrics.csv",
                "pooled_oof_bootstrap_ci.csv",
                "pooled_oof_permutation_p.csv",
                "group_feasibility.csv",
            ]:
                src = cv_out / filename
                dst = out / f"evaluation_ready_{name}_{filename}"
                if _copy_csv(src, dst):
                    outputs["evaluation_ready"][f"{name}_{filename}"] = str(dst)

    tissue_note = out / "exploratory_tissue_evidence_status.csv"
    pd.DataFrame(
        [
            {
                "component": "tissue",
                "status": "exploratory_evidence_layer",
                "reason": "tissue_site_label currently has positives plus unknowns; use PU top-K recovery, not supervised AUROC/AUPRC, until independent negatives/controls are staged.",
            },
            {
                "component": "mechanism",
                "status": "diagnostic_only",
                "reason": "global mechanism labels are source/family sensitive; use graph/top-K/case studies and panel-specific controls before publication-grade supervised mechanism claims.",
            },
        ]
    ).to_csv(tissue_note, index=False)
    outputs["exploratory"]["tissue_mechanism_status"] = str(tissue_note)
    outputs["evaluation_ready"]["metric_coverage"] = _write_metric_coverage(out, outputs=outputs)

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": str(dataset),
        "out_dir": str(out),
        "models": args.models,
        "splits": args.splits,
        "top_k": int(args.top_k),
        "feature_metadata_manifest": feature_metadata_manifest,
        "chemical_cluster_mode": args.chemical_cluster,
        "evaluation_ready_definition": [
            "Binding/exposure outputs use measured SPD labels and leakage-controlled splits.",
            "Files prefixed evaluation_ready_ are intended for validation tables/figures after manual review of caveats.",
            "Grouped/repeated CV with bootstrap/permutation CIs is offline final-validation mode and runs only when --run-grouped-cv is set.",
            "Tissue and mechanism outputs are explicitly exploratory/diagnostic unless independent negatives or panel controls are added.",
        ],
        "outputs": outputs,
        "commands": commands,
        "failed_commands": [row for row in commands if row.get("returncode") != 0],
    }
    (out / "evaluation_ready_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    command_columns = ["command", "log", "returncode", "timed_out", "started_at", "finished_at", "elapsed_seconds"]
    pd.DataFrame(commands, columns=command_columns).to_csv(out / "evaluation_ready_command_log.csv", index=False)
    if manifest["failed_commands"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

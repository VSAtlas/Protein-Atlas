#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional


@dataclass(frozen=True)
class DeltaRow:
    key: str
    delta_roc_auc: float
    delta_pr_auc: float
    source_path: str


def _read_json(path: Path) -> Mapping[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return payload if isinstance(payload, dict) else {}


def _safe_float(value: Any) -> Optional[float]:
    try:
        if value is None:
            return None
        token = str(value).strip()
        if not token:
            return None
        return float(token)
    except Exception:
        return None


def _quantile(values_sorted: List[float], q: float) -> float:
    if not values_sorted:
        return 0.0
    if len(values_sorted) == 1:
        return float(values_sorted[0])
    q_clamped = min(1.0, max(0.0, float(q)))
    idx = q_clamped * float(len(values_sorted) - 1)
    lo = int(math.floor(idx))
    hi = int(math.ceil(idx))
    if lo == hi:
        return float(values_sorted[lo])
    frac = idx - float(lo)
    return float(values_sorted[lo] * (1.0 - frac) + values_sorted[hi] * frac)


def _mean(values: List[float]) -> float:
    return float(sum(values) / float(len(values))) if values else 0.0


def _stddev(values: List[float]) -> float:
    if len(values) <= 1:
        return 0.0
    mu = _mean(values)
    var = sum((v - mu) ** 2 for v in values) / float(len(values) - 1)
    return float(math.sqrt(max(0.0, var)))


def _load_delta_rows(comparison_json: Path) -> List[DeltaRow]:
    payload = _read_json(comparison_json)
    rows_raw = payload.get("quality_vs_baseline")
    if not isinstance(rows_raw, list):
        return []
    out: List[DeltaRow] = []
    for row in rows_raw:
        if not isinstance(row, Mapping):
            continue
        key = str(row.get("key", "")).strip()
        if not key:
            continue
        delta_roc = _safe_float(row.get("delta_roc_auc"))
        delta_pr = _safe_float(row.get("delta_pr_auc"))
        if delta_roc is None or delta_pr is None:
            continue
        out.append(
            DeltaRow(
                key=key,
                delta_roc_auc=float(delta_roc),
                delta_pr_auc=float(delta_pr),
                source_path=str(comparison_json),
            )
        )
    return out


def _resolve_comparison_path(token: str, repo_root: Path) -> Path:
    raw = Path(token)
    if raw.suffix.lower() == ".json" and raw.exists():
        return raw.resolve()
    candidate = repo_root / "analysis" / "benchmarks" / token / "comparison.json"
    return candidate.resolve()


def _build_profile(
    rows: Iterable[DeltaRow],
    *,
    sigma_multiplier: float,
    threshold_floor: float,
) -> Dict[str, Any]:
    by_key: Dict[str, Dict[str, List[float]]] = {}
    for row in rows:
        slot = by_key.setdefault(row.key, {"roc": [], "pr": []})
        slot["roc"].append(float(row.delta_roc_auc))
        slot["pr"].append(float(row.delta_pr_auc))

    per_key: Dict[str, Dict[str, Any]] = {}
    roc_thresholds: List[float] = []
    pr_thresholds: List[float] = []
    for key in sorted(by_key.keys()):
        roc_vals = by_key[key]["roc"]
        pr_vals = by_key[key]["pr"]
        roc_abs_sorted = sorted(abs(v) for v in roc_vals)
        pr_abs_sorted = sorted(abs(v) for v in pr_vals)
        mean_roc = _mean(roc_vals)
        mean_pr = _mean(pr_vals)
        std_roc = _stddev(roc_vals)
        std_pr = _stddev(pr_vals)
        p95_abs_roc = _quantile(roc_abs_sorted, 0.95)
        p95_abs_pr = _quantile(pr_abs_sorted, 0.95)
        thresh_roc = max(
            float(threshold_floor),
            abs(mean_roc) + float(sigma_multiplier) * std_roc,
            p95_abs_roc,
        )
        thresh_pr = max(
            float(threshold_floor),
            abs(mean_pr) + float(sigma_multiplier) * std_pr,
            p95_abs_pr,
        )
        roc_thresholds.append(float(thresh_roc))
        pr_thresholds.append(float(thresh_pr))
        per_key[key] = {
            "n": int(len(roc_vals)),
            "mean_delta_roc_auc": mean_roc,
            "mean_delta_pr_auc": mean_pr,
            "std_delta_roc_auc": std_roc,
            "std_delta_pr_auc": std_pr,
            "p95_abs_delta_roc_auc": p95_abs_roc,
            "p95_abs_delta_pr_auc": p95_abs_pr,
            "threshold_roc_auc": float(thresh_roc),
            "threshold_pr_auc": float(thresh_pr),
        }

    global_profile = {
        "threshold_roc_auc": max(roc_thresholds) if roc_thresholds else float(threshold_floor),
        "threshold_pr_auc": max(pr_thresholds) if pr_thresholds else float(threshold_floor),
    }
    return {
        "schema_version": 1,
        "sigma_multiplier": float(sigma_multiplier),
        "threshold_floor": float(threshold_floor),
        "keys": per_key,
        "global": global_profile,
    }


def _evaluate(profile: Mapping[str, Any], rows: Iterable[DeltaRow]) -> Dict[str, Any]:
    keys = profile.get("keys")
    keys_map = keys if isinstance(keys, Mapping) else {}
    global_cfg = profile.get("global")
    global_map = global_cfg if isinstance(global_cfg, Mapping) else {}
    global_roc = _safe_float(global_map.get("threshold_roc_auc")) or 0.05
    global_pr = _safe_float(global_map.get("threshold_pr_auc")) or 0.05

    out_rows: List[Dict[str, Any]] = []
    flagged = 0
    for row in rows:
        key_cfg_raw = keys_map.get(row.key)
        key_cfg = key_cfg_raw if isinstance(key_cfg_raw, Mapping) else {}
        thr_roc = _safe_float(key_cfg.get("threshold_roc_auc")) or float(global_roc)
        thr_pr = _safe_float(key_cfg.get("threshold_pr_auc")) or float(global_pr)
        abs_roc = abs(float(row.delta_roc_auc))
        abs_pr = abs(float(row.delta_pr_auc))
        pass_roc = abs_roc <= thr_roc
        pass_pr = abs_pr <= thr_pr
        pass_both = pass_roc and pass_pr
        if not pass_both:
            flagged += 1
        out_rows.append(
            {
                "key": row.key,
                "delta_roc_auc": float(row.delta_roc_auc),
                "delta_pr_auc": float(row.delta_pr_auc),
                "abs_delta_roc_auc": abs_roc,
                "abs_delta_pr_auc": abs_pr,
                "threshold_roc_auc": float(thr_roc),
                "threshold_pr_auc": float(thr_pr),
                "pass_roc": bool(pass_roc),
                "pass_pr": bool(pass_pr),
                "pass_both": bool(pass_both),
                "source_path": row.source_path,
            }
        )

    return {
        "rows": out_rows,
        "summary": {
            "checked_rows": int(len(out_rows)),
            "flagged_rows": int(flagged),
            "likely_noise_only": bool(flagged == 0),
        },
    }


def _default_profile_path(repo_root: Path) -> Path:
    return repo_root / "analysis" / "benchmarks" / "bench2_aa_variance_profile.json"


def _default_eval_path(repo_root: Path) -> Path:
    return repo_root / "analysis" / "benchmarks" / "bench2_aa_variance_eval.json"


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build and apply an A/A variance profile for bench2 quality deltas "
            "(ROC_AUC/PR_AUC) from benchmark_bench2_compare comparison artifacts."
        ),
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    fit = sub.add_parser("fit", help="Build an A/A variance profile.")
    fit.add_argument(
        "--aa-comparison",
        action="append",
        required=True,
        help=(
            "Comparison JSON path or run id (resolved to "
            "analysis/benchmarks/<run_id>/comparison.json). Repeat for multiple A/A pairs."
        ),
    )
    fit.add_argument("--sigma-multiplier", type=float, default=3.0)
    fit.add_argument("--threshold-floor", type=float, default=0.02)
    fit.add_argument("--out", help="Output profile JSON path.")

    assess = sub.add_parser("assess", help="Assess candidate compare artifact vs A/A profile.")
    assess.add_argument("--profile", required=True, help="Path to variance profile JSON.")
    assess.add_argument(
        "--candidate-comparison",
        required=True,
        help="Comparison JSON path or run id for candidate.",
    )
    assess.add_argument("--out", help="Output assessment JSON path.")
    return parser


def main() -> int:
    args = _build_parser().parse_args()
    repo_root = Path(__file__).resolve().parents[1]

    if args.cmd == "fit":
        rows: List[DeltaRow] = []
        for token in list(args.aa_comparison):
            path = _resolve_comparison_path(str(token), repo_root)
            rows.extend(_load_delta_rows(path))
        profile = _build_profile(
            rows,
            sigma_multiplier=float(args.sigma_multiplier),
            threshold_floor=float(args.threshold_floor),
        )
        out_path = Path(args.out).resolve() if args.out else _default_profile_path(repo_root)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(profile, indent=2, sort_keys=True), encoding="utf-8")
        print(
            json.dumps(
                {
                    "path": str(out_path),
                    "aa_rows": int(len(rows)),
                    "keys": int(len(profile.get("keys", {}))),
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0

    profile_path = Path(args.profile).resolve()
    profile_payload = _read_json(profile_path)
    candidate_path = _resolve_comparison_path(str(args.candidate_comparison), repo_root)
    candidate_rows = _load_delta_rows(candidate_path)
    report = _evaluate(profile_payload, candidate_rows)
    out_path = Path(args.out).resolve() if args.out else _default_eval_path(repo_root)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    summary = report.get("summary", {})
    print(
        json.dumps(
            {
                "path": str(out_path),
                "checked_rows": int(summary.get("checked_rows", 0)),
                "flagged_rows": int(summary.get("flagged_rows", 0)),
                "likely_noise_only": bool(summary.get("likely_noise_only", False)),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

import csv
import json
import logging
import math
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from analysis import dud_eval

_LOG = logging.getLogger(__name__)
_ROC_AUC_SCORE = getattr(dud_eval, "roc_auc_score")
_EF_AT_FRACTIONS = getattr(dud_eval, "ef_at_fractions")


def _read_scores_csv(path: Path) -> List[Dict[str, Any]]:
    if not path.exists():
        return []
    try:
        if path.stat().st_size == 0:
            return []
    except Exception:
        return []

    rows: List[Dict[str, Any]] = []
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                ligand_id = row.get("ligand_id")
                label = row.get("label")
                score_raw = row.get("score")
                if not ligand_id or label not in ("strong", "non"):
                    continue
                if score_raw is None or score_raw == "":
                    continue
                try:
                    score = float(score_raw)
                except Exception:
                    continue
                rows.append(
                    {"ligand_id": str(ligand_id), "label": label, "score": score}
                )
    except Exception:
        return []
    return rows


def _metric_or_default(value: Any, default: float) -> float:
    if value is None:
        return default
    try:
        val = float(value)
    except Exception:
        return default
    if math.isnan(val):
        return default
    return val


def _compute_basic_metrics(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    y_true = np.array([1 if r["label"] == "strong" else 0 for r in rows], dtype=int)
    scores = np.array([float(r["score"]) for r in rows], dtype=float)
    n_pos = int(y_true.sum())
    n_neg = int(len(y_true) - n_pos)

    metrics: Dict[str, Any] = {
        "n_pos": n_pos,
        "n_neg": n_neg,
        "n_total": int(len(y_true)),
    }

    if n_pos == 0 or n_neg == 0:
        metrics["error"] = "missing_class"
        metrics["auc"] = float("nan")
        metrics["ef1"] = float("nan")
        return metrics

    try:
        metrics["auc"] = float(_ROC_AUC_SCORE(y_true, -scores))
    except Exception:
        metrics["auc"] = float("nan")

    ef = _EF_AT_FRACTIONS(y_true, scores, fractions=(0.01,))
    metrics["ef1"] = float(ef.get("EF@1%", float("nan")))
    return metrics


def _ensure_ranking_metrics(
    metrics: Optional[Dict[str, Any]], rows: List[Dict[str, Any]]
) -> Dict[str, Any]:
    out = dict(metrics or {})
    if "auc_mean" not in out or "ef1_mean" not in out:
        computed = _compute_basic_metrics(rows)
        out.setdefault("auc_mean", computed.get("auc", float("nan")))
        out.setdefault("ef1_mean", computed.get("ef1", float("nan")))
        out.setdefault("auc_std", 0.0)
    return out


def _split_scores(
    rows: List[Dict[str, Any]]
) -> Tuple[Dict[str, float], Dict[str, float]]:
    strong_scores: Dict[str, float] = {}
    non_scores: Dict[str, float] = {}
    for row in rows:
        ligand_id = str(row["ligand_id"])
        score = float(row["score"])
        if row["label"] == "strong":
            strong_scores[ligand_id] = score
        elif row["label"] == "non":
            non_scores[ligand_id] = score
    return strong_scores, non_scores


def _rank_pockets(pockets: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    def _sort_key(entry: Dict[str, Any]) -> Tuple[float, float, float, str]:
        metrics = entry.get("metrics") or {}
        auc = _metric_or_default(metrics.get("auc_mean"), float("-inf"))
        ef = _metric_or_default(metrics.get("ef1_mean"), float("-inf"))
        std = _metric_or_default(metrics.get("auc_std"), float("inf"))
        pocket_id = str(entry.get("pocket_id") or "")
        return (-auc, -ef, std, pocket_id)

    return sorted(pockets, key=_sort_key)


def _assign_strong_ligands(
    best: Dict[str, Any],
    second: Dict[str, Any],
    score_margin_threshold: float,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], Dict[str, int]]:
    best_id = str(best["pocket_id"])
    second_id = str(second["pocket_id"])
    best_scores = best["strong_scores"]
    second_scores = second["strong_scores"]

    strong_ids = sorted(set(best_scores) & set(second_scores))
    assignments: List[Dict[str, Any]] = []
    ambiguous: List[Dict[str, Any]] = []
    count_best = 0
    count_second = 0

    for ligand_id in strong_ids:
        score_best = float(best_scores[ligand_id])
        score_second = float(second_scores[ligand_id])
        margin = abs(score_best - score_second)
        payload = {
            "ligand_id": ligand_id,
            "label": "strong",
            "scores": {best_id: score_best, second_id: score_second},
            "margin": margin,
        }
        if margin < score_margin_threshold:
            payload["assigned_pocket_id"] = None
            ambiguous.append(payload)
            continue

        assigned = best_id if score_best <= score_second else second_id
        payload["assigned_pocket_id"] = assigned
        assignments.append(payload)
        if assigned == best_id:
            count_best += 1
        else:
            count_second += 1

    return (
        assignments,
        ambiguous,
        {
            "count_best": count_best,
            "count_second": count_second,
            "total_non_ambiguous": count_best + count_second,
        },
    )


def _compute_sorted_metrics(
    pocket: Dict[str, Any],
    assigned_ids: List[str],
    ambiguous_count: int,
) -> Dict[str, Any]:
    strong_scores = pocket["strong_scores"]
    non_scores = pocket["non_scores"]

    assigned_scores = [
        float(strong_scores[lid])
        for lid in sorted(assigned_ids)
        if lid in strong_scores
    ]
    n_pos = len(strong_scores)
    n_neg = len(non_scores)
    n_pos_assigned = len(assigned_scores)

    metrics: Dict[str, Any] = {
        "n_pos": n_pos,
        "n_neg": n_neg,
        "n_pos_assigned": n_pos_assigned,
        "n_ambiguous_strong": int(ambiguous_count),
    }

    if n_pos_assigned == 0:
        metrics["error"] = "no_assigned_positives"
        metrics["auc"] = float("nan")
        metrics["ef1"] = float("nan")
        return metrics
    if n_neg == 0:
        metrics["error"] = "missing_class"
        metrics["auc"] = float("nan")
        metrics["ef1"] = float("nan")
        return metrics

    y_true = np.array([1] * n_pos_assigned + [0] * n_neg, dtype=int)
    scores = np.array(
        assigned_scores + [float(non_scores[lid]) for lid in sorted(non_scores)],
        dtype=float,
    )

    try:
        metrics["auc"] = float(_ROC_AUC_SCORE(y_true, -scores))
    except Exception:
        metrics["auc"] = float("nan")
    ef = _EF_AT_FRACTIONS(y_true, scores, fractions=(0.01,))
    metrics["ef1"] = float(ef.get("EF@1%", float("nan")))
    return metrics


def run_second_pass(
    pocket_performance_path: Path,
    *,
    auc_gap_threshold: float = 0.03,
    strong_split_min_fraction: float = 0.30,
    score_margin_threshold: float = 0.50,
    max_pockets: int = 2,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, Any]:
    """
    Returns a dict containing at least:
      - selected_pockets: List[str]
      - multi_pocket_detected: bool
      - paths: {"assignments": "...", "performance_sorted": "..."}
    Writes pocket_assignments.json and pocket_performance_sorted.json into pocket_eval dir.
    """
    log = logger or _LOG
    perf_path = Path(pocket_performance_path)
    if not perf_path.exists():
        raise FileNotFoundError(f"pocket_performance.json missing: {perf_path}")

    try:
        bytes_on_disk = None
        try:
            bytes_on_disk = perf_path.stat().st_size
        except Exception:
            bytes_on_disk = None
        if bytes_on_disk is not None:
            log.info(
                "[pocket_eval.json.read.start] path=%s bytes=%d",
                perf_path,
                bytes_on_disk,
            )
        start = time.perf_counter()
        raw = perf_path.read_text(encoding="utf-8")
        payload = json.loads(raw)
        elapsed_s = time.perf_counter() - start
        keys = list(payload.keys()) if isinstance(payload, dict) else []
        log.info(
            "[pocket_eval.json.read.done] path=%s elapsed_s=%.3f keys=%s",
            perf_path,
            elapsed_s,
            keys[:8],
        )
    except Exception as exc:
        raise ValueError(f"Failed to read pocket performance: {perf_path}") from exc

    pocket_eval_dir = perf_path.parent
    pockets_payload = payload.get("pockets") if isinstance(payload, dict) else None
    if not isinstance(pockets_payload, list):
        pockets_payload = []

    pockets: List[Dict[str, Any]] = []
    for entry in pockets_payload:
        if not isinstance(entry, dict):
            continue
        pocket_id = entry.get("pocket_id")
        if not pocket_id:
            continue
        artifacts_raw = entry.get("artifacts")
        artifacts = artifacts_raw if isinstance(artifacts_raw, dict) else {}
        scores_path = artifacts.get("scores_path")
        path = (
            Path(scores_path)
            if scores_path
            else pocket_eval_dir / str(pocket_id) / "calibration_scores.csv"
        )
        if not path.is_absolute():
            path = (pocket_eval_dir / path).resolve()
        if not path.exists():
            fallback = pocket_eval_dir / str(pocket_id) / "calibration_scores.csv"
            if path != fallback:
                path = fallback

        rows = _read_scores_csv(path)
        if not rows:
            log.warning(
                "[pocket-second-pass] scores missing pocket=%s path=%s",
                pocket_id,
                path,
            )
        metrics = _ensure_ranking_metrics(entry.get("metrics"), rows)
        strong_scores, non_scores = _split_scores(rows)
        pockets.append(
            {
                "pocket_id": str(pocket_id),
                "metrics": metrics,
                "scores_path": path,
                "rows": rows,
                "strong_scores": strong_scores,
                "non_scores": non_scores,
            }
        )

    ranked = _rank_pockets(pockets)
    if not ranked:
        raise ValueError("No pockets available for second pass")

    best = ranked[0]
    second = ranked[1] if len(ranked) > 1 else None

    selected_pockets = [str(best["pocket_id"])]
    multi_pocket_detected = False
    decision_reason = "single_pocket"

    assignments: List[Dict[str, Any]] = []
    ambiguous: List[Dict[str, Any]] = []
    assigned_ids_by_pocket: Dict[str, List[str]] = {}
    ambiguous_count = 0

    if second and max_pockets >= 2:
        best_auc = _metric_or_default(best["metrics"].get("auc_mean"), float("-inf"))
        second_auc = _metric_or_default(
            second["metrics"].get("auc_mean"), float("-inf")
        )
        gap = best_auc - second_auc
        if math.isnan(gap):
            gap = 0.0

        if gap > auc_gap_threshold:
            multi_pocket_detected = False
            decision_reason = "clear_winner"
        else:
            (
                assignments,
                ambiguous,
                split_stats,
            ) = _assign_strong_ligands(best, second, score_margin_threshold)
            ambiguous_count = len(ambiguous)

            total_non_ambiguous = split_stats.get("total_non_ambiguous", 0)
            split_present = False
            if total_non_ambiguous >= 4:
                frac_best = split_stats["count_best"] / total_non_ambiguous
                frac_second = split_stats["count_second"] / total_non_ambiguous
                split_present = (
                    frac_best >= strong_split_min_fraction
                    and frac_second >= strong_split_min_fraction
                )

            multi_pocket_detected = True
            selected_pockets = [str(best["pocket_id"]), str(second["pocket_id"])]
            decision_reason = "strong_split" if split_present else "similar_performance"

            assigned_ids_by_pocket = {selected_pockets[0]: [], selected_pockets[1]: []}
            for entry in assignments:
                assigned_id = entry.get("assigned_pocket_id")
                ligand_id = entry.get("ligand_id")
                if assigned_id in assigned_ids_by_pocket and ligand_id:
                    assigned_ids_by_pocket[assigned_id].append(str(ligand_id))

    if not multi_pocket_detected:
        assigned_ids_by_pocket = {
            selected_pockets[0]: sorted(best["strong_scores"].keys())
        }
        ambiguous_count = 0

    assignments_path = pocket_eval_dir / "pocket_assignments.json"
    assignments_payload: Dict[str, Any] = {
        "selected_pockets": selected_pockets,
        "multi_pocket_detected": multi_pocket_detected,
        "assignments": assignments if multi_pocket_detected else [],
        "ambiguous": ambiguous if multi_pocket_detected else [],
        "assignment_mode": (
            "sorted_strong_split" if multi_pocket_detected else "skipped_single_pocket"
        ),
    }
    assignments_text = json.dumps(assignments_payload, indent=2)
    start = time.perf_counter()
    assignments_path.write_text(assignments_text, encoding="utf-8")
    elapsed_s = time.perf_counter() - start
    log.info(
        "[pocket_eval.json.write.done] path=%s bytes=%d elapsed_s=%.3f keys=%s",
        assignments_path,
        len(assignments_text.encode("utf-8")),
        elapsed_s,
        list(assignments_payload.keys())[:8],
    )

    pocket_by_id = {p["pocket_id"]: p for p in pockets}
    pocket_entries: List[Dict[str, Any]] = []
    for pocket_id in selected_pockets:
        pocket = pocket_by_id.get(pocket_id)
        if not pocket:
            continue
        unsorted_metrics = _compute_basic_metrics(pocket["rows"])
        sorted_metrics = _compute_sorted_metrics(
            pocket,
            assigned_ids_by_pocket.get(pocket_id, []),
            ambiguous_count,
        )
        pocket_entries.append(
            {
                "pocket_id": pocket_id,
                "unsorted_metrics": unsorted_metrics,
                "sorted_metrics": sorted_metrics,
            }
        )

    performance_sorted_path = pocket_eval_dir / "pocket_performance_sorted.json"
    performance_payload = {
        "multi_pocket_detected": multi_pocket_detected,
        "selected_pockets": selected_pockets,
        "decision_reason": decision_reason,
        "pockets": pocket_entries,
    }
    performance_text = json.dumps(performance_payload, indent=2)
    start = time.perf_counter()
    performance_sorted_path.write_text(performance_text, encoding="utf-8")
    elapsed_s = time.perf_counter() - start
    log.info(
        "[pocket_eval.json.write.done] path=%s bytes=%d elapsed_s=%.3f keys=%s",
        performance_sorted_path,
        len(performance_text.encode("utf-8")),
        elapsed_s,
        list(performance_payload.keys())[:8],
    )

    return {
        "selected_pockets": selected_pockets,
        "multi_pocket_detected": multi_pocket_detected,
        "decision_reason": decision_reason,
        "paths": {
            "assignments": str(assignments_path),
            "performance_sorted": str(performance_sorted_path),
        },
    }

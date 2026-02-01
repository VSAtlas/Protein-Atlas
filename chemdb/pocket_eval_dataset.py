from __future__ import annotations

import csv
import hashlib
import json
import logging
import math
import os
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np

from analysis import dud_eval

try:  # optional dependency
    from sklearn.model_selection import StratifiedGroupKFold  # type: ignore[import-untyped]
except Exception:  # pragma: no cover - optional dependency
    StratifiedGroupKFold = None

DATASET_SCHEMA_VERSION = 1

_POCKET_META_KEYS = {
    "rank": "pocket_rank",
    "pocket_rank": "pocket_rank",
    "score": "pocket_score",
    "pocket_score": "pocket_score",
    "residue_summary": "pocket_residue_summary",
    "pocket_residue_summary": "pocket_residue_summary",
    "residue_count": "pocket_residue_count",
    "pocket_residue_count": "pocket_residue_count",
    "num_residues": "pocket_residue_count",
    "volume": "pocket_volume",
    "pocket_volume": "pocket_volume",
    "name": "pocket_name",
    "pocket_name": "pocket_name",
    "probability": "pocket_probability",
    "druggability": "pocket_druggability",
}


def _now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _json_dumps(value: Any) -> str:
    return json.dumps(
        value, ensure_ascii=True, sort_keys=True, separators=(",", ":"), default=str
    )


def _as_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    try:
        return float(value)
    except Exception:
        return None


def _stable_hash(text: str) -> int:
    return int(hashlib.sha1(text.encode("utf-8")).hexdigest()[:16], 16)


def _extract_pocket_meta(pocket: Dict[str, Any]) -> Tuple[Dict[str, Any], str]:
    meta: Dict[str, Any] = {}
    for key, out_key in _POCKET_META_KEYS.items():
        if key in pocket and out_key not in meta:
            val = pocket.get(key)
            if isinstance(val, (list, dict)):
                val = _json_dumps(val)
            meta[out_key] = val

    pocket_meta = dict(pocket)
    coords = pocket_meta.get("coords")
    if isinstance(coords, list) and len(coords) > 100:
        pocket_meta.pop("coords", None)
        pocket_meta["coords_omitted"] = True
    pocket_meta_json = _json_dumps(pocket_meta)
    return meta, pocket_meta_json


def _prep_report_to_sets(
    prep_report: Optional[Dict[str, Any]],
) -> Tuple[set[str], Dict[str, Dict[str, Any]]]:
    if not prep_report:
        return set(), {}
    missing = set(prep_report.get("missing_ligand_ids") or [])
    invalid_list = prep_report.get("invalid_ligands") or []
    invalid_map: Dict[str, Dict[str, Any]] = {}
    for entry in invalid_list:
        ligand_id = entry.get("ligand_id")
        if ligand_id:
            invalid_map[str(ligand_id)] = entry
    return missing, invalid_map


def build_pocket_eval_dataset_long(
    *,
    run_id: str,
    pdb_id: str,
    variant: str,
    calibrators_used: List[Dict[str, Any]],
    pockets_plan: List[Dict[str, Any]],
    scores_by_pocket: Dict[str, List[Dict[str, Any]]],
    dock_events_by_pocket: Optional[Dict[str, List[Dict[str, Any]]]] = None,
    prep_report: Optional[Dict[str, Any]] = None,
) -> List[Dict[str, Any]]:
    """Return long-format rows with one row per (ligand_id, pocket_id)."""
    calibrator_by_id: Dict[str, Dict[str, Any]] = {}
    for row in calibrators_used:
        ligand_id = row.get("ligand_id")
        if ligand_id:
            calibrator_by_id[str(ligand_id)] = row

    missing_ids, invalid_map = _prep_report_to_sets(prep_report)

    score_index: Dict[str, Dict[str, Dict[str, Any]]] = {}
    for pocket_id, rows in scores_by_pocket.items():
        if not rows:
            continue
        score_index[pocket_id] = {
            str(row.get("ligand_id")): row for row in rows if row.get("ligand_id")
        }

    event_index: Dict[str, Dict[str, Dict[str, Any]]] = {}
    if dock_events_by_pocket:
        for pocket_id, events in dock_events_by_pocket.items():
            if not events:
                continue
            event_index[pocket_id] = {
                str(ev.get("ligand_id")): ev for ev in events if ev.get("ligand_id")
            }

    rows_out: List[Dict[str, Any]] = []
    for plan in pockets_plan:
        pocket_id = str(plan.get("pocket_id") or "")
        if not pocket_id:
            continue
        pocket = plan.get("pocket") or {}
        center = plan.get("center") or (None, None, None)
        box_size = plan.get("box_size") or (None, None, None)
        pocket_meta_fields, pocket_meta_json = _extract_pocket_meta(pocket)

        pocket_scores = score_index.get(pocket_id, {})
        pocket_events = event_index.get(pocket_id, {})

        for ligand_id, cal_row in calibrator_by_id.items():
            label = cal_row.get("label")
            label_used = label in {"strong", "non"}

            dock_status = "missing_score"
            score_source = "missing"
            score_val: Optional[float] = None
            dock_ms: Optional[float] = None
            out_path: Optional[str] = None
            dock_error_type: Optional[str] = None
            dock_error: Optional[str] = None

            if ligand_id in missing_ids:
                dock_status = "prep_missing"
            elif ligand_id in invalid_map:
                dock_status = "prep_invalid"
                bad = invalid_map.get(ligand_id, {})
                dock_error_type = bad.get("bad_type")
                dock_error = bad.get("bad_type")
            else:
                event = pocket_events.get(ligand_id)
                if event:
                    dock_status = str(event.get("status") or "error")
                    score_source = str(event.get("score_source") or "run")
                    score_val = _as_float(event.get("score"))
                    dock_ms = _as_float(event.get("dock_ms"))
                    out_path = event.get("out_path")
                    dock_error_type = event.get("error_type")
                    dock_error = event.get("error_str")
                    if dock_status != "scored":
                        score_val = None
                else:
                    cached = pocket_scores.get(ligand_id)
                    if cached and cached.get("score") is not None:
                        score_val = _as_float(cached.get("score"))
                        dock_status = "scored"
                        score_source = "cache"

            row_out = {
                "run_id": run_id,
                "pdb_id": pdb_id,
                "variant": variant,
                "ligand_id": ligand_id,
                "label": label,
                "label_used_for_metrics": bool(label_used),
                "ligand_smiles": cal_row.get("smiles"),
                "label_source": cal_row.get("label_source"),
                "pocket_id": pocket_id,
                "score": score_val,
                "score_source": score_source,
                "dock_status": dock_status,
                "dock_ms": dock_ms,
                "out_path": out_path,
                "dock_error_type": dock_error_type,
                "dock_error": dock_error,
                "pocket_center_x": _as_float(center[0]) if center else None,
                "pocket_center_y": _as_float(center[1]) if center else None,
                "pocket_center_z": _as_float(center[2]) if center else None,
                "pocket_box_x": _as_float(box_size[0]) if box_size else None,
                "pocket_box_y": _as_float(box_size[1]) if box_size else None,
                "pocket_box_z": _as_float(box_size[2]) if box_size else None,
                "pocket_meta_json": pocket_meta_json,
            }
            row_out.update(pocket_meta_fields)
            rows_out.append(row_out)

    return rows_out


def make_pocket_eval_folds(
    rows: List[Dict[str, Any]],
    n_splits: int,
    seed: int,
    group_key: str = "ligand_id",
    stratify_labels: bool = True,
    strategy: str = "auto",
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Assign fold_id and return split metadata."""
    n_splits = max(1, int(n_splits or 1))
    seed = int(seed or 0)
    strategy_requested = str(strategy or "auto").lower()
    group_labels: Dict[str, Optional[str]] = {}

    for row in rows:
        group_id = str(row.get(group_key) or row.get("ligand_id") or "")
        row["group_id"] = group_id
        row["split_id"] = group_id
        label = row.get("label")
        if group_id:
            if label in {"strong", "non"}:
                group_labels[group_id] = label
            else:
                group_labels.setdefault(group_id, None)

    group_ids = sorted(group_labels.keys())
    if not group_ids:
        return rows, {
            "strategy_requested": strategy_requested,
            "strategy_used": "none",
            "seed": seed,
            "n_splits": 0,
            "group_key": group_key,
            "group_to_fold": {},
            "folds": [],
        }

    group_to_fold: Dict[str, int] = {}
    if strategy_requested == "disabled":
        group_to_fold = {gid: 0 for gid in group_ids}
        for row in rows:
            row["fold_id"] = 0
        return rows, _build_splits_payload(
            rows,
            group_to_fold,
            group_labels,
            seed,
            group_key,
            strategy_requested,
            "disabled",
        )

    strategy_used = "hash_group"

    can_stratify = False
    if stratify_labels:
        labels = [group_labels.get(gid) for gid in group_ids]
        if all(label in {"strong", "non"} for label in labels):
            counts = Counter(labels)
            can_stratify = min(counts.values()) >= n_splits

    if (
        StratifiedGroupKFold is not None
        and strategy_requested in {"auto", "sklearn"}
        and can_stratify
    ):
        try:
            y = np.array(
                [1 if group_labels[gid] == "strong" else 0 for gid in group_ids]
            )
            groups = np.array(group_ids)
            splitter = StratifiedGroupKFold(
                n_splits=n_splits, shuffle=True, random_state=seed
            )
            for fold_id, (_, test_idx) in enumerate(
                splitter.split(np.zeros(len(groups)), y, groups=groups)
            ):
                for idx in test_idx:
                    group_to_fold[str(groups[idx])] = int(fold_id)
            strategy_used = "sklearn_stratified_group"
        except Exception:
            group_to_fold.clear()

    if not group_to_fold:
        for gid in group_ids:
            fold_id = int(_stable_hash(f"{gid}:{seed}") % n_splits)
            group_to_fold[gid] = fold_id
        strategy_used = "hash_group"

    for row in rows:
        group_id = row.get("group_id") or row.get(group_key) or ""
        row["fold_id"] = int(group_to_fold.get(str(group_id), 0))

    splits_payload = _build_splits_payload(
        rows,
        group_to_fold,
        group_labels,
        seed,
        group_key,
        strategy_requested,
        strategy_used,
    )
    splits_payload["n_splits"] = n_splits
    return rows, splits_payload


def _build_splits_payload(
    rows: List[Dict[str, Any]],
    group_to_fold: Dict[str, int],
    group_labels: Dict[str, Optional[str]],
    seed: int,
    group_key: str,
    strategy_requested: str,
    strategy_used: str,
) -> Dict[str, Any]:
    fold_stats: Dict[int, Dict[str, Any]] = {}
    for gid, fold_id in group_to_fold.items():
        stats = fold_stats.setdefault(
            fold_id,
            {
                "fold_id": fold_id,
                "n_groups": 0,
                "n_strong": 0,
                "n_non": 0,
                "n_unknown": 0,
            },
        )
        stats["n_groups"] += 1
        label = group_labels.get(gid)
        if label == "strong":
            stats["n_strong"] += 1
        elif label == "non":
            stats["n_non"] += 1
        else:
            stats["n_unknown"] += 1

    for row in rows:
        fold_id_raw = row.get("fold_id")
        if fold_id_raw is None:
            continue
        try:
            fold_id = int(fold_id_raw)
        except Exception:
            continue
        stats = fold_stats.setdefault(
            fold_id,
            {
                "fold_id": fold_id,
                "n_groups": 0,
                "n_strong": 0,
                "n_non": 0,
                "n_unknown": 0,
            },
        )
        stats["n_rows"] = int(stats.get("n_rows", 0)) + 1

    return {
        "strategy_requested": strategy_requested,
        "strategy_used": strategy_used,
        "seed": seed,
        "group_key": group_key,
        "group_to_fold": group_to_fold,
        "folds": sorted(fold_stats.values(), key=lambda d: d["fold_id"]),
    }


def _atomic_path(path: Path) -> Path:
    token = f".tmp_{os.getpid()}_{int(time.time() * 1000)}"
    return path.with_name(path.name + token)


def _write_text_atomic(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = _atomic_path(path)
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)


def _write_json_atomic(path: Path, payload: Dict[str, Any]) -> None:
    text = json.dumps(payload, indent=2, ensure_ascii=True)
    _write_text_atomic(path, text)


def _write_csv_atomic(
    path: Path, rows: Iterable[Dict[str, Any]], fieldnames: List[str]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = _atomic_path(path)
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    tmp.replace(path)


def _build_fieldnames(rows: List[Dict[str, Any]]) -> List[str]:
    base_fields = [
        "run_id",
        "pdb_id",
        "variant",
        "ligand_id",
        "ligand_smiles",
        "label",
        "label_used_for_metrics",
        "label_source",
        "group_id",
        "split_id",
        "fold_id",
        "pocket_id",
        "score",
        "score_source",
        "dock_status",
        "dock_ms",
        "out_path",
        "dock_error_type",
        "dock_error",
        "pocket_center_x",
        "pocket_center_y",
        "pocket_center_z",
        "pocket_box_x",
        "pocket_box_y",
        "pocket_box_z",
        "pocket_rank",
        "pocket_score",
        "pocket_volume",
        "pocket_residue_summary",
        "pocket_residue_count",
        "pocket_name",
        "pocket_probability",
        "pocket_druggability",
        "pocket_meta_json",
    ]
    extra = sorted({k for row in rows for k in row.keys() if k not in base_fields})
    return base_fields + extra


def _dataset_counts(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    status_counts: Counter[str] = Counter()
    source_counts: Counter[str] = Counter()
    for row in rows:
        status_counts[str(row.get("dock_status"))] += 1
        source_counts[str(row.get("score_source"))] += 1
    return {
        "n_rows": len(rows),
        "n_ligands": len({row.get("ligand_id") for row in rows}),
        "n_pockets": len({row.get("pocket_id") for row in rows}),
        "dock_status_counts": dict(status_counts),
        "score_source_counts": dict(source_counts),
    }


def _pocket_rows(rows: List[Dict[str, Any]]) -> Dict[str, List[Dict[str, Any]]]:
    by_pocket: Dict[str, List[Dict[str, Any]]] = {}
    for row in rows:
        pocket_id = row.get("pocket_id")
        if not pocket_id:
            continue
        by_pocket.setdefault(str(pocket_id), []).append(row)
    return by_pocket


def _ligand_rows(rows: List[Dict[str, Any]]) -> Dict[str, List[Dict[str, Any]]]:
    by_ligand: Dict[str, List[Dict[str, Any]]] = {}
    for row in rows:
        ligand_id = row.get("ligand_id")
        if not ligand_id:
            continue
        by_ligand.setdefault(str(ligand_id), []).append(row)
    return by_ligand


def _auc_and_ef(y_true: List[int], scores_low: List[float]) -> Tuple[float, float]:
    if not y_true or len(y_true) != len(scores_low):
        return float("nan"), float("nan")
    if min(sum(y_true), len(y_true) - sum(y_true)) <= 0:
        return float("nan"), float("nan")
    y_arr = np.array(y_true, dtype=int)
    score_arr = np.array(scores_low, dtype=float)
    try:
        auc_val = float(dud_eval.roc_auc_score(y_arr, -score_arr))  # type: ignore[attr-defined]
    except Exception:
        auc_val = float("nan")
    ef = dud_eval.ef_at_fractions(y_arr, score_arr, fractions=(0.01,))  # type: ignore[attr-defined]
    ef1 = float(ef.get("EF@1%", float("nan")))
    return auc_val, ef1


def _compute_ensemble_metrics(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    by_ligand = _ligand_rows(rows)
    y_true: List[int] = []
    scores_low: List[float] = []
    scored_ligands = 0

    for lig_id, lig_rows in by_ligand.items():
        del lig_id
        label = None
        scores = []
        for row in lig_rows:
            if row.get("label") in {"strong", "non"}:
                label = row.get("label")
            score = row.get("score")
            if score is not None:
                scores.append(float(score))
        if scores:
            scored_ligands += 1
            if label in {"strong", "non"}:
                y_true.append(1 if label == "strong" else 0)
                scores_low.append(min(scores))

    auc_val, ef1 = _auc_and_ef(y_true, scores_low)
    return {
        "ensemble_min_score_auc": auc_val,
        "ensemble_min_score_ef1": ef1,
        "ligands_with_any_score": scored_ligands,
        "ligands_with_labels": len(y_true),
    }


def _compute_cv_assigned(
    rows: List[Dict[str, Any]],
    *,
    top_m: int,
    require_scores: bool,
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    by_pocket = _pocket_rows(rows)
    by_ligand = _ligand_rows(rows)
    fold_ids: List[int] = []
    for row in rows:
        fold_val = row.get("fold_id")
        if fold_val is None:
            continue
        try:
            fold_ids.append(int(fold_val))
        except Exception:
            continue
    folds = sorted(set(fold_ids))
    per_ligand_rows: List[Dict[str, Any]] = []
    fold_metrics: Dict[int, Dict[str, Any]] = {}

    for fold_id in folds:
        pocket_auc: List[Tuple[str, float]] = []
        for pocket_id, pocket_rows in by_pocket.items():
            y_true = []
            scores_low = []
            for row in pocket_rows:
                if row.get("fold_id") == fold_id:
                    continue
                if not row.get("label_used_for_metrics"):
                    continue
                score = row.get("score")
                if score is None:
                    continue
                y_true.append(1 if row.get("label") == "strong" else 0)
                scores_low.append(float(score))
            auc_val, _ = _auc_and_ef(y_true, scores_low)
            if math.isfinite(auc_val):
                pocket_auc.append((pocket_id, auc_val))
        pocket_auc.sort(key=lambda item: item[1], reverse=True)
        candidates = [pid for pid, _ in pocket_auc[: max(1, int(top_m))]]

        fold_y: List[int] = []
        fold_scores: List[float] = []
        for lig_id, lig_rows in by_ligand.items():
            if not lig_rows:
                continue
            if lig_rows[0].get("fold_id") != fold_id:
                continue
            label = lig_rows[0].get("label")
            if label not in {"strong", "non"}:
                continue
            best_score = None
            best_pocket = None
            for row in lig_rows:
                if row.get("pocket_id") not in candidates:
                    continue
                score = row.get("score")
                if score is None:
                    continue
                if best_score is None or float(score) < best_score:
                    best_score = float(score)
                    best_pocket = row.get("pocket_id")
            if best_score is None and require_scores:
                continue

            per_ligand_rows.append(
                {
                    "ligand_id": lig_id,
                    "label": label,
                    "fold_id": fold_id,
                    "cv_candidate_pockets": _json_dumps(candidates),
                    "cv_assigned_pocket": best_pocket,
                    "cv_assigned_score": best_score,
                }
            )
            if best_score is not None:
                fold_y.append(1 if label == "strong" else 0)
                fold_scores.append(float(best_score))

        auc_val, ef1 = _auc_and_ef(fold_y, fold_scores)
        fold_metrics[fold_id] = {
            "fold_id": fold_id,
            "n_assigned": len(fold_scores),
            "cv_assigned_auc": auc_val,
            "cv_assigned_ef1": ef1,
            "cv_candidate_pockets": candidates,
        }

    auc_vals = [m["cv_assigned_auc"] for m in fold_metrics.values()]
    ef_vals = [m["cv_assigned_ef1"] for m in fold_metrics.values()]
    return (
        {
            "cv_assigned_auc_mean": float(np.nanmean(auc_vals))
            if auc_vals
            else float("nan"),
            "cv_assigned_ef1_mean": float(np.nanmean(ef_vals))
            if ef_vals
            else float("nan"),
            "per_fold": list(fold_metrics.values()),
        },
        per_ligand_rows,
    )


def _compute_oracle(
    rows: List[Dict[str, Any]],
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    by_ligand = _ligand_rows(rows)
    fold_ids: List[int] = []
    for row in rows:
        fold_val = row.get("fold_id")
        if fold_val is None:
            continue
        try:
            fold_ids.append(int(fold_val))
        except Exception:
            continue
    folds = sorted(set(fold_ids))
    per_ligand_rows: List[Dict[str, Any]] = []
    fold_metrics: Dict[int, Dict[str, Any]] = {}

    for fold_id in folds:
        fold_y: List[int] = []
        fold_scores: List[float] = []
        for lig_id, lig_rows in by_ligand.items():
            if not lig_rows:
                continue
            if lig_rows[0].get("fold_id") != fold_id:
                continue
            label = lig_rows[0].get("label")
            if label not in {"strong", "non"}:
                continue
            scored = []
            for row in lig_rows:
                score = row.get("score")
                if score is not None:
                    scored.append((float(score), row.get("pocket_id")))
            if not scored:
                continue
            if label == "strong":
                best_score, best_pocket = min(scored, key=lambda item: item[0])
            else:
                best_score, best_pocket = max(scored, key=lambda item: item[0])
            per_ligand_rows.append(
                {
                    "ligand_id": lig_id,
                    "label": label,
                    "fold_id": fold_id,
                    "oracle_score": best_score,
                    "oracle_pocket": best_pocket,
                }
            )
            fold_y.append(1 if label == "strong" else 0)
            fold_scores.append(best_score)

        auc_val, ef1 = _auc_and_ef(fold_y, fold_scores)
        fold_metrics[fold_id] = {
            "fold_id": fold_id,
            "n_scored": len(fold_scores),
            "oracle_auc": auc_val,
            "oracle_ef1": ef1,
        }

    auc_vals = [m["oracle_auc"] for m in fold_metrics.values()]
    ef_vals = [m["oracle_ef1"] for m in fold_metrics.values()]
    return (
        {
            "oracle_auc_mean": float(np.nanmean(auc_vals))
            if auc_vals
            else float("nan"),
            "oracle_ef1_mean": float(np.nanmean(ef_vals)) if ef_vals else float("nan"),
            "upper_bound_note": "Oracle uses labels to pick best pocket per ligand; not a screening estimate.",
            "per_fold": list(fold_metrics.values()),
        },
        per_ligand_rows,
    )


def _maybe_write_parquet(
    path: Path, rows: List[Dict[str, Any]], logger: logging.Logger
) -> bool:
    try:
        import pandas as pd  # type: ignore
    except Exception as exc:
        logger.info("[pocket-eval] parquet.skip reason=import_error err=%s", exc)
        return False

    try:
        df = pd.DataFrame(rows)
        tmp = _atomic_path(path)
        df.to_parquet(tmp, index=False)
        tmp.replace(path)
        return True
    except Exception as exc:
        logger.info("[pocket-eval] parquet.skip reason=write_error err=%s", exc)
        return False


def _maybe_write_wide_csv(
    path: Path, rows: List[Dict[str, Any]], logger: logging.Logger
) -> bool:
    try:
        import pandas as pd  # type: ignore
    except Exception as exc:
        logger.info("[pocket-eval] wide_csv.skip reason=import_error err=%s", exc)
        return False

    try:
        df = pd.DataFrame(rows)
        if "score" not in df.columns or "pocket_id" not in df.columns:
            return False
        wide = df.pivot_table(
            index=["ligand_id", "label"],
            columns="pocket_id",
            values="score",
            aggfunc="min",
        )
        wide = wide.reset_index()
        tmp = _atomic_path(path)
        wide.to_csv(tmp, index=False)
        tmp.replace(path)
        return True
    except Exception as exc:
        logger.info("[pocket-eval] wide_csv.skip reason=write_error err=%s", exc)
        return False


def build_pockets_used(pockets_plan: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    pockets_used: List[Dict[str, Any]] = []
    for plan in pockets_plan:
        pocket = plan.get("pocket") or {}
        pocket_id = plan.get("pocket_id")
        center = plan.get("center")
        box_size = plan.get("box_size")
        pocket_meta_fields, pocket_meta_json = _extract_pocket_meta(pocket)
        pockets_used.append(
            {
                "pocket_id": pocket_id,
                "center": center,
                "box_size": box_size,
                "pocket_meta_json": pocket_meta_json,
                **pocket_meta_fields,
            }
        )
    return pockets_used


def write_pocket_eval_artifacts(
    dock_root: Path,
    rows: List[Dict[str, Any]],
    *,
    config_snapshot: Dict[str, Any],
    pockets_used: List[Dict[str, Any]],
    splits_payload: Dict[str, Any],
    logger: Optional[logging.Logger] = None,
    formats: str = "csv",
    enable_wide: bool = False,
    ensemble_enable: bool = False,
    cv_enable: bool = False,
    cv_top_m: int = 2,
    cv_require_scores: bool = True,
    oracle_enable: bool = False,
) -> Dict[str, Path]:
    log = logger or logging.getLogger(__name__)
    dock_root.mkdir(parents=True, exist_ok=True)

    dataset_path = dock_root / "pocket_eval_dataset_long.csv"
    fieldnames = _build_fieldnames(rows)
    _write_csv_atomic(dataset_path, rows, fieldnames)

    pockets_used_path = dock_root / "pocket_eval_pockets_used.json"
    _write_json_atomic(pockets_used_path, {"pockets": pockets_used})

    splits_path = dock_root / "pocket_eval_splits.json"
    _write_json_atomic(splits_path, splits_payload)

    meta_payload = {
        "schema_version": DATASET_SCHEMA_VERSION,
        "created_at": _now_iso(),
        "config_snapshot": config_snapshot,
        "counts": _dataset_counts(rows),
        "split_strategy": {
            "strategy_requested": splits_payload.get("strategy_requested"),
            "strategy_used": splits_payload.get("strategy_used"),
            "n_splits": splits_payload.get("n_splits"),
            "seed": splits_payload.get("seed"),
            "group_key": splits_payload.get("group_key"),
        },
        "label_definitions": {
            "positive": "strong",
            "negative": "non",
            "ignored": ["weak"],
        },
        "score_convention": "lower_is_better",
    }

    meta_path = dock_root / "pocket_eval_dataset_meta.json"
    _write_json_atomic(meta_path, meta_payload)

    artifacts: Dict[str, Path] = {
        "pocket_eval_dataset_long.csv": dataset_path,
        "pocket_eval_dataset_meta.json": meta_path,
        "pocket_eval_pockets_used.json": pockets_used_path,
        "pocket_eval_splits.json": splits_path,
    }

    formats_set = {
        f.strip().lower() for f in str(formats or "").split(",") if f.strip()
    }
    if "parquet" in formats_set:
        parquet_path = dock_root / "pocket_eval_dataset_long.parquet"
        if _maybe_write_parquet(parquet_path, rows, log):
            artifacts["pocket_eval_dataset_long.parquet"] = parquet_path

    if enable_wide:
        wide_path = dock_root / "pocket_eval_dataset_wide.csv"
        if _maybe_write_wide_csv(wide_path, rows, log):
            artifacts["pocket_eval_dataset_wide.csv"] = wide_path

    if ensemble_enable:
        ensemble_path = dock_root / "pocket_eval_ensemble_metrics.json"
        ensemble_payload = _compute_ensemble_metrics(rows)
        _write_json_atomic(ensemble_path, ensemble_payload)
        artifacts["pocket_eval_ensemble_metrics.json"] = ensemble_path

    if cv_enable:
        cv_payload, cv_rows = _compute_cv_assigned(
            rows, top_m=cv_top_m, require_scores=cv_require_scores
        )
        cv_path = dock_root / "pocket_eval_cv_assigned.json"
        _write_json_atomic(cv_path, cv_payload)
        artifacts["pocket_eval_cv_assigned.json"] = cv_path
        cv_rows_path = dock_root / "pocket_eval_cv_assigned_per_ligand.csv"
        if cv_rows:
            cv_fieldnames = _build_fieldnames(cv_rows)
            _write_csv_atomic(cv_rows_path, cv_rows, cv_fieldnames)
            artifacts["pocket_eval_cv_assigned_per_ligand.csv"] = cv_rows_path

    if oracle_enable:
        oracle_payload, oracle_rows = _compute_oracle(rows)
        oracle_path = dock_root / "pocket_eval_oracle_upper_bound.json"
        _write_json_atomic(oracle_path, oracle_payload)
        artifacts["pocket_eval_oracle_upper_bound.json"] = oracle_path
        if oracle_rows:
            oracle_rows_path = dock_root / "pocket_eval_oracle_per_ligand.csv"
            oracle_fieldnames = _build_fieldnames(oracle_rows)
            _write_csv_atomic(oracle_rows_path, oracle_rows, oracle_fieldnames)
            artifacts["pocket_eval_oracle_per_ligand.csv"] = oracle_rows_path

    return artifacts

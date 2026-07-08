# -*- coding: utf-8 -*-
import bisect
import logging
import math
from pathlib import Path
from typing import Any, Dict, List, Tuple

from analysis.reporting.fdr_methods import compute_score_fdr
from analysis.reporting.value_utils import as_float
from analysis.reporting.master_schema.constants import (
    COMPONENT,
    MIN_DECOYS_FOR_FDR,
    MIN_UNIQUE_DECOY_SCORES,
)
from analysis.reporting.master_schema.decoy_helpers import _row_is_decoy
from analysis.reporting.master_schema.fdr_scores import (
    _get_decoy_scores_for_combo,
    _get_docking_decoy_scores_for_combo,
)


def _apply_fdr_to_master_rows(
    *,
    master_rows: List[Dict[str, Any]],
    post_root: Path,
    docked_root: Path,
    run_id: str,
    decoy_prefix: str,
    logger: logging.Logger,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    # Initialize FDR fields
    for r in master_rows:
        r.setdefault("fdr_score_field", "")
        r.setdefault("fdr_null_source", "")
        r.setdefault("fdr_n_decoys", "")
        r.setdefault("fdr_unique_decoy_scores", "")
        r.setdefault("fdr_reliable", "")
        r.setdefault("fdr_p_empirical", "")
        r.setdefault("fdr_q_target", "")
        r.setdefault("fdr_hit_q05", "")
        r.setdefault("fdr_hit_q10", "")
        r.setdefault("scorch_fdr_scope", "")
        r.setdefault("scorch_fdr_n_decoys", "")
        r.setdefault("scorch_fdr_n_tested", "")
        r.setdefault("scorch_fdr_p_empirical", "")
        r.setdefault("scorch_fdr_q_bh", "")
        r.setdefault("scorch_fdr_decoy_competition_q", "")
        r.setdefault("scorch_fdr_decoy_competition_q_plus1", "")
        r.setdefault("scorch_fdr_decoy_competition_hit_q05", "")
        r.setdefault("scorch_fdr_decoy_competition_hit_q10", "")
        r.setdefault("scorch_fdr_decoy_competition_hit_q05_plus1", "")
        r.setdefault("scorch_fdr_decoy_competition_hit_q10_plus1", "")

    combo_indices: Dict[Tuple[str, str, str], List[int]] = {}
    for idx, row in enumerate(master_rows):
        key = (row["pdb_id"], row["variant"], row["ph_label"])
        combo_indices.setdefault(key, []).append(idx)

    fdr_summary_rows: List[Dict[str, Any]] = []
    scorch_fdr_summary_rows: List[Dict[str, Any]] = []
    for (pdb_id, variant, ph), indices in combo_indices.items():
        combo_dir = post_root / run_id / pdb_id / variant / ph
        primary_decoy_path = combo_dir / f"{decoy_prefix}_consensus_reranked_scorch.csv"
        legacy_decoy_path = combo_dir / "dud_consensus_reranked_scorch.csv"
        fallback_rows = [master_rows[i] for i in indices]
        (
            field,
            decoy_scores,
            used_fallback,
            _decoy_scores_by_field,
            fdr_null_source,
        ) = _get_decoy_scores_for_combo(
            combo_dir,
            fallback_rows,
            MIN_DECOYS_FOR_FDR,
            decoy_prefix=decoy_prefix,
        )
        if len(decoy_scores) < MIN_DECOYS_FOR_FDR:
            docking_combo_dir = docked_root / run_id / pdb_id / variant / ph
            (
                docking_field,
                docking_scores,
                _docking_scores_by_field,
                docking_source_path,
            ) = _get_docking_decoy_scores_for_combo(
                docking_combo_dir,
                MIN_DECOYS_FOR_FDR,
                decoy_prefix=decoy_prefix,
            )
            if len(docking_scores) > len(decoy_scores):
                field = docking_field
                decoy_scores = docking_scores
                fdr_null_source = "docked_consensus_decoy"
                logger.info(
                    "%s action=fdr_null_source_fallback pdb=%s variant=%s ph=%s source=docked_consensus_decoy path=%s n_decoys=%d",
                    COMPONENT,
                    pdb_id,
                    variant,
                    ph,
                    docking_source_path,
                    len(docking_scores),
                )
        n_decoys = len(decoy_scores)
        decoys_sorted = sorted(decoy_scores)
        unique_scores = len(set(decoy_scores))
        reliable = (
            n_decoys >= MIN_DECOYS_FOR_FDR and unique_scores >= MIN_UNIQUE_DECOY_SCORES
        )

        if used_fallback:
            log_path = (
                primary_decoy_path if primary_decoy_path.exists() else legacy_decoy_path
            )
            logger.warning(
                "%s action=fdr_decoy_fallback pdb=%s variant=%s ph=%s path=%s reason=%s",
                COMPONENT,
                pdb_id,
                variant,
                ph,
                log_path,
                "missing_or_empty" if not log_path.exists() else "no_decoys_found",
            )

        field_for_rows = field or "consensus_score"

        non_decoy_indices = [
            idx for idx in indices if not _row_is_decoy(master_rows[idx], decoy_prefix)
        ]
        for idx in non_decoy_indices:
            master_rows[idx]["fdr_score_field"] = field_for_rows
            master_rows[idx]["fdr_null_source"] = fdr_null_source
            master_rows[idx]["fdr_n_decoys"] = str(n_decoys)
            master_rows[idx]["fdr_unique_decoy_scores"] = str(unique_scores)
            master_rows[idx]["fdr_reliable"] = "1" if reliable else "0"

        tested: List[Tuple[int, float]] = []
        tested_scores: List[Tuple[int, float]] = []
        for idx in non_decoy_indices:
            row = master_rows[idx]
            score_val = as_float(row.get(field_for_rows))
            if score_val is None or not math.isfinite(score_val):
                continue
            n_ge = n_decoys - bisect.bisect_left(decoys_sorted, score_val)
            p_val = (1 + n_ge) / (1 + n_decoys)
            tested.append((idx, p_val))
            tested_scores.append((idx, float(score_val)))
            row["fdr_p_empirical"] = f"{p_val:.6g}"

        hits05 = 0
        hits10 = 0
        decoy_competition_hits05 = 0
        decoy_competition_hits10 = 0
        best_q = 1.0
        best_decoy_competition_q = 1.0
        tested_sorted: List[Tuple[int, float]] = []
        n_tested_rows = 0
        min_possible_bh_q = 1.0
        if tested:
            tested_sorted = sorted(tested, key=lambda x: x[1])
            m = len(tested_sorted)
            n_tested_rows = m
            min_possible_bh_q = min(1.0, (1.0 / (1.0 + n_decoys)) * m)
            raw_q: List[Tuple[int, float]] = []
            for rank, (idx, p_val) in enumerate(tested_sorted, start=1):
                raw_q.append((idx, p_val * m / rank))

            q_values: Dict[int, float] = {}
            prev = 1.0
            for idx, q_val in reversed(raw_q):
                adj = min(q_val, prev, 1.0)
                prev = adj
                q_values[idx] = adj

            for idx, q_val in q_values.items():
                master_rows[idx]["fdr_q_target"] = f"{q_val:.6g}"
                hit05 = q_val <= 0.05
                hit10 = q_val <= 0.10
                master_rows[idx]["fdr_hit_q05"] = "1" if hit05 else "0"
                master_rows[idx]["fdr_hit_q10"] = "1" if hit10 else "0"
                if hit05:
                    hits05 += 1
                if hit10:
                    hits10 += 1
                best_q = min(best_q, q_val)

            competition_raw: List[Tuple[int, float]] = []
            for rank, (idx, score_val) in enumerate(
                sorted(tested_scores, key=lambda item: item[1], reverse=True),
                start=1,
            ):
                n_decoys_ge = n_decoys - bisect.bisect_left(decoys_sorted, score_val)
                competition_raw.append((idx, n_decoys_ge / max(1, rank)))

            competition_q: Dict[int, float] = {}
            running = 1.0
            for idx, value in reversed(competition_raw):
                running = min(running, value, 1.0)
                competition_q[idx] = running
            for idx, q_val in competition_q.items():
                master_rows[idx]["fdr_decoy_competition_q"] = f"{q_val:.6g}"
                hit05 = q_val <= 0.05
                hit10 = q_val <= 0.10
                master_rows[idx]["fdr_decoy_competition_hit_q05"] = "1" if hit05 else "0"
                master_rows[idx]["fdr_decoy_competition_hit_q10"] = "1" if hit10 else "0"
                if hit05:
                    decoy_competition_hits05 += 1
                if hit10:
                    decoy_competition_hits10 += 1
                best_decoy_competition_q = min(best_decoy_competition_q, q_val)
        else:
            n_tested_rows = len(non_decoy_indices)

        for idx in non_decoy_indices:
            if not master_rows[idx].get("fdr_p_empirical"):
                master_rows[idx]["fdr_p_empirical"] = "1"
            if not master_rows[idx].get("fdr_q_target"):
                master_rows[idx]["fdr_q_target"] = "1"
            if not master_rows[idx].get("fdr_hit_q05"):
                master_rows[idx]["fdr_hit_q05"] = "0"
            if not master_rows[idx].get("fdr_hit_q10"):
                master_rows[idx]["fdr_hit_q10"] = "0"
            if not master_rows[idx].get("fdr_decoy_competition_q"):
                master_rows[idx]["fdr_decoy_competition_q"] = "1"
            if not master_rows[idx].get("fdr_decoy_competition_hit_q05"):
                master_rows[idx]["fdr_decoy_competition_hit_q05"] = "0"
            if not master_rows[idx].get("fdr_decoy_competition_hit_q10"):
                master_rows[idx]["fdr_decoy_competition_hit_q10"] = "0"

        fdr_summary_rows.append(
            {
                "pdb_id": pdb_id,
                "variant": variant,
                "ph_label": ph,
                "fdr_score_field": field_for_rows,
                "fdr_null_source": fdr_null_source,
                "fdr_n_decoys": str(n_decoys),
                "fdr_unique_decoy_scores": str(unique_scores),
                "fdr_reliable": "1" if reliable else "0",
                "n_tested": str(n_tested_rows),
                "n_hits_q05": str(hits05),
                "n_hits_q10": str(hits10),
                "best_q": f"{best_q:.6g}" if math.isfinite(best_q) else "",
                "min_possible_bh_q": f"{min_possible_bh_q:.6g}",
                "bh_resolvable_q05": "1" if min_possible_bh_q <= 0.05 else "0",
                "bh_resolvable_q10": "1" if min_possible_bh_q <= 0.10 else "0",
                "n_decoy_competition_hits_q05": str(decoy_competition_hits05),
                "n_decoy_competition_hits_q10": str(decoy_competition_hits10),
                "best_decoy_competition_q": f"{best_decoy_competition_q:.6g}"
                if math.isfinite(best_decoy_competition_q)
                else "",
            }
        )

        scorch_decoy_scores: List[float] = []
        scorch_tested_scores: List[Tuple[int, float]] = []
        for idx in indices:
            row = master_rows[idx]
            score_val = as_float(row.get("SCORCH_score_used") or row.get("scorch_composite"))
            if score_val is None or not math.isfinite(score_val):
                continue
            if _row_is_decoy(row, decoy_prefix):
                scorch_decoy_scores.append(float(score_val))
            else:
                scorch_tested_scores.append((idx, float(score_val)))

        scorch_fdr = compute_score_fdr(
            scorch_tested_scores,
            scorch_decoy_scores,
            higher_is_better=True,
        )
        if scorch_tested_scores or scorch_decoy_scores:
            s_summary = scorch_fdr.summary
            for idx, score_metrics in scorch_fdr.row_metrics.items():
                row = master_rows[idx]
                row["scorch_fdr_scope"] = "conditional_post_scorch"
                row["scorch_fdr_n_decoys"] = str(int(s_summary["n_decoys"]))
                row["scorch_fdr_n_tested"] = str(int(s_summary["n_tested"]))
                row["scorch_fdr_p_empirical"] = f"{score_metrics['p_empirical']:.6g}"
                row["scorch_fdr_q_bh"] = f"{score_metrics['q_bh']:.6g}"
                row["scorch_fdr_decoy_competition_q"] = (
                    f"{score_metrics['decoy_competition_q']:.6g}"
                )
                row["scorch_fdr_decoy_competition_q_plus1"] = (
                    f"{score_metrics['decoy_competition_q_plus1']:.6g}"
                )
                row["scorch_fdr_decoy_competition_hit_q05"] = (
                    "1" if score_metrics["decoy_competition_q"] <= 0.05 else "0"
                )
                row["scorch_fdr_decoy_competition_hit_q10"] = (
                    "1" if score_metrics["decoy_competition_q"] <= 0.10 else "0"
                )
                row["scorch_fdr_decoy_competition_hit_q05_plus1"] = (
                    "1" if score_metrics["decoy_competition_q_plus1"] <= 0.05 else "0"
                )
                row["scorch_fdr_decoy_competition_hit_q10_plus1"] = (
                    "1" if score_metrics["decoy_competition_q_plus1"] <= 0.10 else "0"
                )
            scorch_fdr_summary_rows.append(
                {
                    "pdb_id": pdb_id,
                    "variant": variant,
                    "ph_label": ph,
                    "score_field": "SCORCH_score_used",
                    "scope": "conditional_post_scorch",
                    "n_decoys": str(int(s_summary["n_decoys"])),
                    "n_tested": str(int(s_summary["n_tested"])),
                    "best_bh_q": f"{s_summary['best_bh_q']:.6g}",
                    "min_possible_bh_q": f"{s_summary['min_possible_bh_q']:.6g}",
                    "best_decoy_competition_q": (
                        f"{s_summary['best_decoy_competition_q']:.6g}"
                    ),
                    "best_decoy_competition_q_plus1": (
                        f"{s_summary['best_decoy_competition_q_plus1']:.6g}"
                    ),
                    "n_decoy_competition_hits_q05": str(
                        int(s_summary.get("n_decoy_competition_hits_q05", 0))
                    ),
                    "n_decoy_competition_hits_q10": str(
                        int(s_summary.get("n_decoy_competition_hits_q10", 0))
                    ),
                    "n_decoy_competition_hits_q05_plus1": str(
                        int(s_summary.get("n_decoy_competition_hits_q05_plus1", 0))
                    ),
                    "n_decoy_competition_hits_q10_plus1": str(
                        int(s_summary.get("n_decoy_competition_hits_q10_plus1", 0))
                    ),
                }
            )

        logger.info(
            "%s action=fdr_compute pdb=%s variant=%s ph=%s field=%s n_decoys=%d n_tested=%d n_hits_q05=%d",
            COMPONENT,
            pdb_id,
            variant,
            ph,
            field_for_rows,
            n_decoys,
            n_tested_rows,
            hits05,
        )
    return master_rows, fdr_summary_rows, scorch_fdr_summary_rows

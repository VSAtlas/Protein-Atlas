from __future__ import annotations

import os
from pathlib import Path
from typing import List, Optional, Set

import numpy as np
import pandas as pd

from analysis.dud_eval_core.labels import _collapse_best_scores, parse_name_and_label
from analysis.dud_eval_core.log import dbg
from analysis.dud_eval_core.reporting import _compute_metrics_and_plots, _write_pretty_summary
from analysis.dud_eval_core.schema import (
    _filter_consensus_no_data,
    filter_df_by_run_id,
    guess_consensus_ligfile_col,
    guess_consensus_score_col,
    guess_ligfile_col,
    guess_reranked_scorch_score_col,
    guess_score_col,
    parse_valid_mask,
    read_reranked_scorch_csv,
    resolve_valid_col,
)
from analysis.dud_eval_core.types import TargetEvaluation

def evaluate_target(
    pdb_id: str,
    csv_path: Path,
    out_dir: Path,
    lig_col_cli: Optional[str],
    score_col_cli: Optional[str],
    bedroc_alpha: float,
    logauc_lambda: float,
    run_id: Optional[str],
    valid_only: bool = False,
    valid_col_cli: Optional[str] = None,
    new_layout_active: bool = False,
    layout_label: str = "legacy",
) -> Optional[TargetEvaluation]:
    status_reason = "ok"
    try:
        df = pd.read_csv(csv_path)
        dbg("INFO", "csv", f"pdb={pdb_id} path={csv_path} rows={len(df)}")
    except Exception as e:
        dbg("ERROR", "csv", f"pdb={pdb_id} path={csv_path} err={e}")
        return TargetEvaluation(
            metrics=None,
            ligand_basenames=set(),
            has_run_id_column=False,
            run_ids=set(),
            status_reason="read_error",
        )

    has_run_id_col = "run_id" in df.columns
    df, filter_reason = filter_df_by_run_id(
        df,
        run_id,
        pdb_id=pdb_id,
        csv_path=csv_path,
        strict=new_layout_active,
        layout_label=layout_label,
    )
    status_reason = filter_reason if filter_reason != "ok" else "ok"
    if df.empty:
        final_reason = status_reason if status_reason != "ok" else "empty_after_filter"
        return TargetEvaluation(
            metrics=None,
            ligand_basenames=set(),
            has_run_id_column=has_run_id_col,
            run_ids=set(),
            status_reason=final_reason,
        )

    dbg(
        "DEBUG",
        "schema",
        f"pdb={pdb_id} cols={list(df.columns)} has_variant={'variant' in df.columns} has_run_id={'run_id' in df.columns}",
    )

    lig_col = guess_ligfile_col(df, lig_col_cli)
    score_col = guess_score_col(df, score_col_cli)
    dbg(
        "DEBUG",
        "schema",
        f"pdb={pdb_id} ligand_col={lig_col} source={'CLI' if lig_col_cli else 'auto'} score_col={score_col} source={'CLI' if score_col_cli else 'auto'}",
    )
    # Coerce scores to numeric; drop NaN/±inf early to avoid NaNs in metrics
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    before_nf = len(df)
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df = df.dropna(subset=[score_col])
    drop_nonfinite = before_nf - len(df)
    if df.empty:
        final_reason = status_reason if status_reason != "ok" else "empty_after_filter"
        return TargetEvaluation(
            metrics=None,
            ligand_basenames=set(),
            has_run_id_column=has_run_id_col,
            run_ids=set(),
            status_reason=final_reason,
        )

    # Apply pose validity filtering before best-pose aggregation.
    if valid_only:
        valid_col = resolve_valid_col(df, valid_col_cli, enabled=True)
        before_valid = len(df)
        valid_mask = parse_valid_mask(df[valid_col])
        df = df.loc[valid_mask].copy()
        dropped_valid = before_valid - len(df)
        dbg(
            "INFO",
            "filter",
            f"pdb={pdb_id} valid_only=ON valid_col={valid_col} rows_kept={len(df)} rows_dropped={dropped_valid}",
        )

    raw_lig = df[lig_col]
    missing_names = int(raw_lig.isna().sum())
    empty_names = int(raw_lig.astype(str).str.strip().eq("").sum())
    # derive ligand_id + label from filename
    lig_paths = raw_lig.astype(str)
    parsed = lig_paths.apply(parse_name_and_label)
    df["lig_id"] = parsed.apply(lambda t: t[0])
    df["is_active"] = parsed.apply(lambda t: t[1])
    used_ligand_basenames: Set[str] = {
        os.path.basename(p) for p in lig_paths.tolist() if p
    }

    run_ids_present: Set[str] = set()
    if has_run_id_col and "run_id" in df.columns:
        try:
            run_ids_present = set(df["run_id"].dropna().astype(str).unique().tolist())
        except Exception:
            run_ids_present = set()

    # drop rows where label couldn't be inferred
    before = len(df)
    df = df.dropna(subset=["is_active"])
    dropped_token = before - len(df)
    if df.empty:
        final_reason = status_reason if status_reason != "ok" else "empty_after_filter"
        return TargetEvaluation(
            metrics=None,
            ligand_basenames=used_ligand_basenames,
            has_run_id_column=has_run_id_col,
            run_ids=run_ids_present,
            status_reason=final_reason,
        )

    best, drop_nan_best = _collapse_best_scores(df, score_col, best_is_min=True)
    missing_name_detected = missing_names + empty_names

    variant_label: Optional[str] = None
    if "variant" in df.columns:
        variant_series = df["variant"].dropna()
        if not variant_series.empty:
            normalized = variant_series.astype(str).str.strip()
            normalized = normalized[normalized != ""]
            if not normalized.empty:
                upper = normalized.str.upper()
                has_apo = bool((upper == "APO").any())
                has_holo = bool((upper == "HOLO").any())
                if has_apo and not has_holo:
                    variant_label = "APO"
                elif has_holo and not has_apo:
                    variant_label = "HOLO"
                elif has_apo or has_holo:
                    variant_label = "MIXED"
                unique_vals = sorted(set(upper.tolist()))
                counts = {val: int((upper == val).sum()) for val in unique_vals}
                dbg(
                    "DEBUG",
                    "variant",
                    f"pdb={pdb_id} unique={unique_vals} counts={counts}",
                )

    return _compute_metrics_and_plots(
        pdb_id=pdb_id,
        best=best,
        out_dir=out_dir,
        bedroc_alpha=bedroc_alpha,
        logauc_lambda=logauc_lambda,
        score_high_is_better=False,
        missing_name_detected=missing_name_detected,
        drop_nonfinite=drop_nonfinite,
        dropped_token=dropped_token,
        drop_nan_best=drop_nan_best,
        token_regex="active|decoy",
        variant_label=variant_label,
        run_id=run_id,
        ligand_basenames=used_ligand_basenames,
        has_run_id_col=has_run_id_col,
        run_ids_present=run_ids_present,
        hist_xlabel="Best docking score (lower = better)",
        title_suffix="",
        mode_label="docking",
        status_reason=status_reason,
    )

def evaluate_target_consensus(
    pdb_id: str,
    csv_path: Path,
    out_dir: Path,
    lig_col_cli: Optional[str],
    score_col_cli: Optional[str],
    bedroc_alpha: float,
    logauc_lambda: float,
    run_id: Optional[str],
    variant_hint: Optional[str] = None,
    ph_tag: Optional[str] = None,  # kept for symmetry/future use
    new_layout_active: bool = False,
    layout_label: str = "legacy",
) -> Optional[TargetEvaluation]:
    status_reason = "ok"
    try:
        df = pd.read_csv(csv_path)
        dbg("INFO", "consensus.csv", f"pdb={pdb_id} path={csv_path} rows={len(df)}")
    except Exception as e:
        dbg("ERROR", "consensus.csv", f"pdb={pdb_id} path={csv_path} err={e}")
        return TargetEvaluation(
            metrics=None,
            ligand_basenames=set(),
            has_run_id_column=False,
            run_ids=set(),
            status_reason="read_error",
        )

    has_run_id_col = "run_id" in df.columns
    df, filter_reason = filter_df_by_run_id(
        df,
        run_id,
        pdb_id=pdb_id,
        csv_path=csv_path,
        strict=new_layout_active,
        layout_label=layout_label,
    )
    status_reason = filter_reason if filter_reason != "ok" else "ok"
    if df.empty:
        final_reason = status_reason if status_reason != "ok" else "empty_after_filter"
        return TargetEvaluation(
            metrics=None,
            ligand_basenames=set(),
            has_run_id_column=has_run_id_col,
            run_ids=set(),
            status_reason=final_reason,
        )

    lig_col = guess_consensus_ligfile_col(df, lig_col_cli)
    score_col = guess_consensus_score_col(df, score_col_cli)
    dbg(
        "DEBUG",
        "consensus.schema",
        f"pdb={pdb_id} ligand_col={lig_col} score_col={score_col} cols={list(df.columns)}",
    )

    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df, dropped_no_data = _filter_consensus_no_data(df, score_col)

    before_nf = len(df)
    df = df.dropna(subset=[score_col])
    drop_nonfinite = before_nf - len(df)
    if df.empty:
        final_reason = status_reason if status_reason != "ok" else "empty_after_filter"
        return TargetEvaluation(
            metrics=None,
            ligand_basenames=set(),
            has_run_id_column=has_run_id_col,
            run_ids=set(),
            status_reason=final_reason,
        )

    raw_lig = df[lig_col]
    missing_names = int(raw_lig.isna().sum())
    empty_names = int(raw_lig.astype(str).str.strip().eq("").sum())
    lig_paths = raw_lig.astype(str)
    parsed = lig_paths.apply(parse_name_and_label)
    df["lig_id"] = parsed.apply(lambda t: t[0])
    df["is_active"] = parsed.apply(lambda t: t[1])
    used_ligand_basenames: Set[str] = {
        os.path.basename(p) for p in lig_paths.tolist() if p
    }

    run_ids_present: Set[str] = set()
    if has_run_id_col and "run_id" in df.columns:
        try:
            run_ids_present = set(df["run_id"].dropna().astype(str).unique().tolist())
        except Exception:
            run_ids_present = set()

    before = len(df)
    df = df.dropna(subset=["is_active"])
    dropped_token = before - len(df)
    if df.empty:
        final_reason = status_reason if status_reason != "ok" else "empty_after_filter"
        return TargetEvaluation(
            metrics=None,
            ligand_basenames=used_ligand_basenames,
            has_run_id_column=has_run_id_col,
            run_ids=run_ids_present,
            status_reason=final_reason,
        )

    best, drop_nan_best = _collapse_best_scores(df, score_col, best_is_min=False)
    missing_name_detected = missing_names + empty_names

    variant_label: Optional[str] = variant_hint.upper() if variant_hint else None
    if variant_label is None and "variant" in df.columns:
        variant_series = df["variant"].dropna()
        if not variant_series.empty:
            normalized = variant_series.astype(str).str.strip()
            normalized = normalized[normalized != ""]
            if not normalized.empty:
                upper = normalized.str.upper()
                has_apo = bool((upper == "APO").any())
                has_holo = bool((upper == "HOLO").any())
                if has_apo and not has_holo:
                    variant_label = "APO"
                elif has_holo and not has_apo:
                    variant_label = "HOLO"
                elif has_apo or has_holo:
                    variant_label = "MIXED"
                unique_vals = sorted(set(upper.tolist()))
                counts = {val: int((upper == val).sum()) for val in unique_vals}
                dbg(
                    "DEBUG",
                    "consensus.variant",
                    f"pdb={pdb_id} unique={unique_vals} counts={counts}",
                )

    dbg(
        "DEBUG",
        "consensus.filter",
        f"pdb={pdb_id} drop_no_data={dropped_no_data} drop_nonfinite={drop_nonfinite} drop_missing_token={dropped_token} drop_nan_best={drop_nan_best}",
    )

    return _compute_metrics_and_plots(
        pdb_id=pdb_id,
        best=best,
        out_dir=out_dir,
        bedroc_alpha=bedroc_alpha,
        logauc_lambda=logauc_lambda,
        score_high_is_better=True,
        missing_name_detected=missing_name_detected,
        drop_nonfinite=drop_nonfinite,
        dropped_token=dropped_token,
        drop_nan_best=drop_nan_best,
        token_regex="active|decoy",
        variant_label=variant_label,
        run_id=run_id,
        ligand_basenames=used_ligand_basenames,
        has_run_id_col=has_run_id_col,
        run_ids_present=run_ids_present,
        hist_xlabel="Consensus score (higher = better; negated for plot)",
        title_suffix=" (consensus)",
        mode_label="consensus",
        status_reason=status_reason,
    )

def evaluate_target_post_docked_reranked_scorch(
    pdb_id: str,
    csv_path: Path,
    out_dir: Path,
    lig_col_cli: Optional[str],
    score_col_cli: Optional[str],
    bedroc_alpha: float,
    logauc_lambda: float,
    run_id: Optional[str],
    variant_hint: Optional[str] = None,
    ph_tag: Optional[str] = None,
    new_layout_active: bool = False,
    layout_label: str = "post_docked",
) -> Optional[TargetEvaluation]:
    status_reason = "ok"
    try:
        df = read_reranked_scorch_csv(csv_path)
        dbg("INFO", "reranked.csv", f"pdb={pdb_id} path={csv_path} rows={len(df)}")
    except Exception as e:
        dbg("ERROR", "reranked.csv", f"pdb={pdb_id} path={csv_path} err={e}")
        return TargetEvaluation(
            metrics=None,
            ligand_basenames=set(),
            has_run_id_column=False,
            run_ids=set(),
            status_reason="read_error",
        )

    has_run_id_col = "run_id" in df.columns
    df, filter_reason = filter_df_by_run_id(
        df,
        run_id,
        pdb_id=pdb_id,
        csv_path=csv_path,
        strict=new_layout_active,
        layout_label=layout_label,
    )
    status_reason = filter_reason if filter_reason != "ok" else "ok"
    if df.empty:
        final_reason = status_reason if status_reason != "ok" else "empty_after_filter"
        return TargetEvaluation(
            metrics=None,
            ligand_basenames=set(),
            has_run_id_column=has_run_id_col,
            run_ids=set(),
            status_reason=final_reason,
        )

    lig_col = guess_consensus_ligfile_col(df, lig_col_cli)
    score_col = guess_reranked_scorch_score_col(df, score_col_cli)
    dbg(
        "DEBUG",
        "reranked.schema",
        f"pdb={pdb_id} ligand_col={lig_col} score_col={score_col} cols={list(df.columns)}",
    )

    if score_col in ("scorch_composite", "SCORCH_score_used"):
        dbg(
            "WARN",
            "post_scorch.score_col",
            f"pdb={pdb_id} score_col={score_col} may evaluate rescored subset only if rescoring was partial; prefer final_score for full-library metrics",
        )

    rescored_frac = None
    if "rescored_flag" in df.columns:
        try:
            rescored_series = pd.to_numeric(df["rescored_flag"], errors="coerce")
            rescored_frac = rescored_series.fillna(0).mean()
        except Exception:
            rescored_frac = None

    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    score_eval_col = "_score_eval"
    if score_col.lower() == "final_rank":
        df[score_eval_col] = -pd.to_numeric(df[score_col], errors="coerce")
    else:
        df[score_eval_col] = df[score_col]

    before_nf = len(df)
    df = df.dropna(subset=[score_eval_col])
    drop_nonfinite = before_nf - len(df)
    if rescored_frac is not None and rescored_frac < 0.999:
        dbg(
            "INFO",
            "post_scorch.partial",
            f"pdb={pdb_id} rescored_frac={rescored_frac:.4f} score_col={score_col}",
        )
        if score_col != "final_score":
            dbg(
                "WARN",
                "post_scorch.subset",
                f"pdb={pdb_id} score_col={score_col} may drop unrescored ligands; consider --score-col final_score",
            )
    elif score_col != "final_score" and drop_nonfinite > 0:
        dbg(
            "WARN",
            "post_scorch.subset",
            f"pdb={pdb_id} score_col={score_col} dropped={drop_nonfinite} rows without scores",
        )

    if df.empty:
        final_reason = status_reason if status_reason != "ok" else "empty_after_filter"
        return TargetEvaluation(
            metrics=None,
            ligand_basenames=set(),
            has_run_id_column=has_run_id_col,
            run_ids=set(),
            status_reason=final_reason,
        )

    raw_lig = df[lig_col]
    missing_names = int(raw_lig.isna().sum())
    empty_names = int(raw_lig.astype(str).str.strip().eq("").sum())
    lig_paths = raw_lig.astype(str)
    parsed = lig_paths.apply(parse_name_and_label)
    df["lig_id"] = parsed.apply(lambda t: t[0])
    df["is_active"] = parsed.apply(lambda t: t[1])
    used_ligand_basenames: Set[str] = {
        os.path.basename(p) for p in lig_paths.tolist() if p
    }

    run_ids_present: Set[str] = set()
    if has_run_id_col and "run_id" in df.columns:
        try:
            run_ids_present = set(df["run_id"].dropna().astype(str).unique().tolist())
        except Exception:
            run_ids_present = set()

    before = len(df)
    df = df.dropna(subset=["is_active"])
    dropped_token = before - len(df)
    if df.empty:
        final_reason = status_reason if status_reason != "ok" else "empty_after_filter"
        return TargetEvaluation(
            metrics=None,
            ligand_basenames=used_ligand_basenames,
            has_run_id_column=has_run_id_col,
            run_ids=run_ids_present,
            status_reason=final_reason,
        )

    best, drop_nan_best = _collapse_best_scores(df, score_eval_col, best_is_min=False)
    missing_name_detected = missing_names + empty_names

    variant_label: Optional[str] = variant_hint.upper() if variant_hint else None
    if variant_label is None and "variant" in df.columns:
        variant_series = df["variant"].dropna()
        if not variant_series.empty:
            normalized = variant_series.astype(str).str.strip()
            normalized = normalized[normalized != ""]
            if not normalized.empty:
                upper = normalized.str.upper()
                has_apo = bool((upper == "APO").any())
                has_holo = bool((upper == "HOLO").any())
                if has_apo and not has_holo:
                    variant_label = "APO"
                elif has_holo and not has_apo:
                    variant_label = "HOLO"
                elif has_apo or has_holo:
                    variant_label = "MIXED"
                unique_vals = sorted(set(upper.tolist()))
                counts = {val: int((upper == val).sum()) for val in unique_vals}
                dbg(
                    "DEBUG",
                    "reranked.variant",
                    f"pdb={pdb_id} unique={unique_vals} counts={counts}",
                )

    return _compute_metrics_and_plots(
        pdb_id=pdb_id,
        best=best,
        out_dir=out_dir,
        bedroc_alpha=bedroc_alpha,
        logauc_lambda=logauc_lambda,
        score_high_is_better=True,
        missing_name_detected=missing_name_detected,
        drop_nonfinite=drop_nonfinite,
        dropped_token=dropped_token,
        drop_nan_best=drop_nan_best,
        token_regex="active|decoy",
        variant_label=variant_label,
        run_id=run_id,
        ligand_basenames=used_ligand_basenames,
        has_run_id_col=has_run_id_col,
        run_ids_present=run_ids_present,
        hist_xlabel="Reranked SCORCH score (higher = better)",
        title_suffix=" (reranked_scorch)",
        mode_label="post_docked_scorch",
        status_reason=status_reason,
    )

def emit_consensus_summary(
    df: pd.DataFrame, analysis_root: Path, run_label: str, pretty_enabled: bool
) -> None:
    if df is None or df.empty:
        dbg("INFO", "consensus.summary", "no consensus rows to write")
        return

    preferred_cols = (
        "variant",
        "pH",
        "run_id",
        "target_name",
        "library_name",
        "pdb_id",
    )
    summary_path = analysis_root / f"consensus_summary_{run_label}.tsv"
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    variant_col_present = "variant" in df.columns
    variant_upper = (
        df["variant"].astype(str).str.upper() if variant_col_present else None
    )
    apo_mask: pd.Series | None = (variant_upper == "APO") if variant_col_present else None
    holo_mask: pd.Series | None = (variant_upper == "HOLO") if variant_col_present else None
    has_sections = bool(
        variant_col_present
        and (
            (apo_mask is not None and bool(apo_mask.any()))
            or (holo_mask is not None and bool(holo_mask.any()))
        )
    )

    with open(summary_path, "w", newline="") as fh:
        sections: List[tuple[str | None, pd.DataFrame]] = []
        if has_sections:
            _df = df.copy()
            cols = [c for c in preferred_cols if c in _df.columns] + [
                c for c in _df.columns if c not in preferred_cols
            ]
            apo_rows = (
                df.loc[apo_mask].sort_values("pdb_id")
                if apo_mask is not None
                else pd.DataFrame()
            )
            holo_rows = (
                df.loc[holo_mask].sort_values("pdb_id")
                if holo_mask is not None
                else pd.DataFrame()
            )
            n_apo = len(apo_rows)
            n_holo = len(holo_rows)

            if n_apo:
                apo_out = apo_rows
                a_cols = [c for c in preferred_cols if c in apo_out.columns] + [
                    c for c in apo_out.columns if c not in preferred_cols
                ]
                fh.write("Apo\n")
                apo_out[a_cols].to_csv(fh, sep="\t", index=False, float_format="%.3f")
                sections.append(("CONSENSUS - Apo", apo_out[a_cols]))

            if n_holo:
                if n_apo:
                    fh.write("\n")
                holo_out = holo_rows
                h_cols = [c for c in preferred_cols if c in holo_out.columns] + [
                    c for c in holo_out.columns if c not in preferred_cols
                ]
                fh.write("Holo\n")
                holo_out[h_cols].to_csv(fh, sep="\t", index=False, float_format="%.3f")
                sections.append(("CONSENSUS - Holo", holo_out[h_cols]))

            leftover_mask = (
                ~(apo_mask | holo_mask)
                if (apo_mask is not None and holo_mask is not None)
                else pd.Series(False, index=df.index)
            )
            leftover_rows = (
                df.loc[leftover_mask].sort_values("pdb_id")
                if not df.empty
                else pd.DataFrame()
            )
            if not leftover_rows.empty:
                if n_apo or n_holo:
                    fh.write("\n")
                lo_base = leftover_rows
                l_cols = [c for c in preferred_cols if c in lo_base.columns] + [
                    c for c in lo_base.columns if c not in preferred_cols
                ]
                lo_base[l_cols].to_csv(fh, sep="\t", index=False)
                sections.append(("CONSENSUS - Unlabeled", lo_base[l_cols]))

            dbg(
                "INFO",
                "consensus.summary",
                f"sections=Apo:{n_apo} Holo:{n_holo} other={len(leftover_rows)} out={summary_path}",
            )
            if pretty_enabled and sections:
                pretty_path = summary_path.with_name(summary_path.stem + "_pretty.txt")
                _write_pretty_summary(sections, pretty_path)
                print(f"[consensus] pretty_summary_out={pretty_path}")
        else:
            _df = df.copy()
            cols = [c for c in preferred_cols if c in _df.columns] + [
                c for c in _df.columns if c not in preferred_cols
            ]
            _df[cols].to_csv(fh, sep="\t", index=False)
            dbg("INFO", "consensus.summary", f"out={summary_path} rows={len(_df)}")
            if pretty_enabled:
                pretty_path = summary_path.with_name(summary_path.stem + "_pretty.txt")
                _write_pretty_summary([("CONSENSUS", _df[cols])], pretty_path)
                print(f"[consensus] pretty_summary_out={pretty_path}")

def emit_reranked_scorch_summary(
    df: pd.DataFrame, analysis_root: Path, run_label: str, pretty_enabled: bool
) -> None:
    if df is None or df.empty:
        dbg("INFO", "reranked.summary", "no reranked rows to write")
        return

    preferred_cols = (
        "variant",
        "pH",
        "run_id",
        "target_name",
        "library_name",
        "pdb_id",
    )
    summary_root = analysis_root / "post_docked"
    summary_path = summary_root / f"consensus_reranked_scorch_summary_{run_label}.tsv"
    summary_root.mkdir(parents=True, exist_ok=True)

    variant_col_present = "variant" in df.columns
    variant_upper = (
        df["variant"].astype(str).str.upper() if variant_col_present else None
    )
    apo_mask: pd.Series | None = (variant_upper == "APO") if variant_col_present else None
    holo_mask: pd.Series | None = (variant_upper == "HOLO") if variant_col_present else None
    has_sections = bool(
        variant_col_present
        and (
            (apo_mask is not None and bool(apo_mask.any()))
            or (holo_mask is not None and bool(holo_mask.any()))
        )
    )

    with open(summary_path, "w", newline="") as fh:
        sections: List[tuple[str | None, pd.DataFrame]] = []
        if has_sections:
            apo_rows = (
                df.loc[apo_mask].sort_values("pdb_id")
                if apo_mask is not None
                else pd.DataFrame()
            )
            holo_rows = (
                df.loc[holo_mask].sort_values("pdb_id")
                if holo_mask is not None
                else pd.DataFrame()
            )
            n_apo = len(apo_rows)
            n_holo = len(holo_rows)

            if n_apo:
                apo_out = apo_rows
                a_cols = [c for c in preferred_cols if c in apo_out.columns] + [
                    c for c in apo_out.columns if c not in preferred_cols
                ]
                fh.write("Apo\n")
                apo_out[a_cols].to_csv(fh, sep="\t", index=False, float_format="%.3f")
                sections.append(("POST_DOCKED_SCORCH - Apo", apo_out[a_cols]))

            if n_holo:
                if n_apo:
                    fh.write("\n")
                holo_out = holo_rows
                h_cols = [c for c in preferred_cols if c in holo_out.columns] + [
                    c for c in holo_out.columns if c not in preferred_cols
                ]
                fh.write("Holo\n")
                holo_out[h_cols].to_csv(fh, sep="\t", index=False, float_format="%.3f")
                sections.append(("POST_DOCKED_SCORCH - Holo", holo_out[h_cols]))

            leftover_mask = (
                ~(apo_mask | holo_mask)
                if (apo_mask is not None and holo_mask is not None)
                else pd.Series(False, index=df.index)
            )
            leftover_rows = (
                df.loc[leftover_mask].sort_values("pdb_id")
                if not df.empty
                else pd.DataFrame()
            )
            if not leftover_rows.empty:
                if n_apo or n_holo:
                    fh.write("\n")
                lo_base = leftover_rows
                l_cols = [c for c in preferred_cols if c in lo_base.columns] + [
                    c for c in lo_base.columns if c not in preferred_cols
                ]
                lo_base[l_cols].to_csv(fh, sep="\t", index=False)
                sections.append(("POST_DOCKED_SCORCH - Unlabeled", lo_base[l_cols]))

            dbg(
                "INFO",
                "reranked.summary",
                f"sections=Apo:{n_apo} Holo:{n_holo} other={len(leftover_rows)} out={summary_path}",
            )
            if pretty_enabled and sections:
                pretty_path = summary_path.with_name(summary_path.stem + "_pretty.txt")
                _write_pretty_summary(sections, pretty_path)
                print(f"[post_docked_scorch] pretty_summary_out={pretty_path}")
        else:
            _df = df.copy()
            cols = [c for c in preferred_cols if c in _df.columns] + [
                c for c in _df.columns if c not in preferred_cols
            ]
            _df[cols].to_csv(fh, sep="\t", index=False)
            dbg("INFO", "reranked.summary", f"out={summary_path} rows={len(_df)}")
            if pretty_enabled:
                pretty_path = summary_path.with_name(summary_path.stem + "_pretty.txt")
                _write_pretty_summary([("POST_DOCKED_SCORCH", _df[cols])], pretty_path)
                print(f"[post_docked_scorch] pretty_summary_out={pretty_path}")

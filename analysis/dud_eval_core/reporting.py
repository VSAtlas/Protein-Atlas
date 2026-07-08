from __future__ import annotations

from pathlib import Path
from typing import Optional, Set

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve

import analysis.dud_eval_core.log as dud_eval_log
from analysis.dud_eval_core.log import dbg
from analysis.dud_eval_core.metrics import bedroc, ef_at_fractions, log_auc_from_roc, pr_auc
from analysis.dud_eval_core.types import TargetEvaluation

def _fmt_counts_and_round(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with big counts comma-formatted and floats rounded to 3 dp (for display only)."""
    out = df.copy()
    # Comma-format counts if present
    for c in ("N", "n_actives"):
        if c in out.columns:
            out[c] = out[c].map(lambda v: f"{int(v):,}" if pd.notna(v) else "")
    # Round floats to 3 decimals (display only)
    for c in out.select_dtypes(include=["float", "float64"]).columns:
        if c not in (
            "actives_fraction",
        ):  # optional: keep full precision here if you prefer
            out[c] = out[c].map(lambda v: f"{v:.3f}" if pd.notna(v) else "")
    return out

def _df_to_pretty_text(df: pd.DataFrame, title: str | None = None) -> str:
    """
    Produce an aligned, readable table using pandas' to_string.
    Assumes df is already ordered and formatted for display.
    """
    # Choose widths implicitly; pandas to_string aligns numbers right by default
    body = df.to_string(index=False)
    return (title + "\n" + body) if title else body

def _write_pretty_summary(
    sections: list[tuple[str | None, pd.DataFrame]], out_path: Path
) -> None:
    """
    sections: list of (title, df) pairs. title can be None for single-table case.
    Writes a single text file with optional section headers separated by blank lines.
    """
    parts = []
    for title, df in sections:
        if df is None or df.empty:
            continue
        df_disp = _fmt_counts_and_round(df)
        parts.append(_df_to_pretty_text(df_disp, title=title))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n\n".join(parts) + ("\n" if parts else ""))

def _write_pretty_table_noformat(
    df: pd.DataFrame, out_path: Path, title: str | None = None
) -> None:
    """
    Write a readable, aligned text table without changing any values.
    This preserves rows, columns, order, and raw numeric formatting.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        if title:
            f.write(f"{title}\n")
        f.write(df.to_string(index=False))
        f.write("\n")

def _compute_metrics_and_plots(
    *,
    pdb_id: str,
    best: pd.DataFrame,
    out_dir: Path,
    bedroc_alpha: float,
    logauc_lambda: float,
    score_high_is_better: bool,
    missing_name_detected: int,
    drop_nonfinite: int,
    dropped_token: int,
    drop_nan_best: int,
    token_regex: str,
    variant_label: Optional[str],
    run_id: Optional[str],
    ligand_basenames: Set[str],
    has_run_id_col: bool,
    run_ids_present: Set[str],
    hist_xlabel: str,
    title_suffix: str = "",
    mode_label: str = "docking",
    status_reason: str = "ok",
) -> TargetEvaluation:
    selected_cols = [
        col for col in ("lig_id", "is_active", "best_score") if col in best.columns
    ]
    z_selected = best[selected_cols].copy() if selected_cols else best.copy()

    y_true = best["is_active"].astype(int).to_numpy()
    if score_high_is_better:
        y_high = best["best_score"].to_numpy()
        y_low = -y_high
    else:
        y_low = best["best_score"].to_numpy()
        y_high = -y_low

    N = len(best)
    n_act = int(y_true.sum())
    dbg(
        "DEBUG",
        "screen",
        f"pdb={pdb_id} mode={mode_label} missing_name_detected={missing_name_detected} kept_rows={len(best)} ligands={N}",
    )
    if N == 0 or n_act == 0 or n_act == N:
        dbg(
            "WARN",
            "metrics",
            f"pdb={pdb_id} mode={mode_label} degenerate_set N={N} actives={n_act}",
        )
        reason = status_reason if status_reason != "ok" else "degenerate_set"
        return TargetEvaluation(
            metrics=None,
            ligand_basenames=ligand_basenames,
            has_run_id_column=has_run_id_col,
            run_ids=run_ids_present,
            status_reason=reason,
            z_selected=z_selected,
        )

    ef = ef_at_fractions(y_true, y_low, fractions=(0.01, 0.02, 0.05, 0.10))
    dbg(
        "DEBUG",
        "metrics",
        f"pdb={pdb_id} mode={mode_label} start N={N} n_actives={n_act}",
    )
    try:
        rocAUC = float(roc_auc_score(y_true, y_high))
    except Exception as exc:
        dbg(
            "WARN",
            "metrics",
            f"pdb={pdb_id} mode={mode_label} metric=ROC_AUC err={exc}",
        )
        rocAUC = float("nan")
    try:
        fpr, tpr, _ = roc_curve(y_true, y_high)
    except Exception as exc:
        dbg(
            "WARN",
            "metrics",
            f"pdb={pdb_id} mode={mode_label} metric=ROC_curve err={exc}",
        )
        fpr = np.array([0.0, 1.0])
        tpr = np.array([0.0, 1.0])
    try:
        prAUC, precision, recall = pr_auc(y_true, y_high)
    except Exception as exc:
        dbg(
            "WARN", "metrics", f"pdb={pdb_id} mode={mode_label} metric=PR_AUC err={exc}"
        )
        prAUC = float("nan")
        precision = np.array([0.0, 1.0])
        recall = np.array([0.0, 1.0])
    try:
        lAUC = log_auc_from_roc(fpr, tpr, lam=logauc_lambda)
    except Exception as exc:
        dbg(
            "WARN", "metrics", f"pdb={pdb_id} mode={mode_label} metric=logAUC err={exc}"
        )
        lAUC = float("nan")
    lAUC_adj = lAUC - 0.14462
    try:
        bed = bedroc(y_true, y_high, alpha=bedroc_alpha)
    except Exception as exc:
        dbg(
            "WARN", "metrics", f"pdb={pdb_id} mode={mode_label} metric=BEDROC err={exc}"
        )
        bed = float("nan")

    out_dir.mkdir(parents=True, exist_ok=True)
    if not dud_eval_log._BACKEND_LOGGED:
        try:
            backend = plt.get_backend()
        except Exception as exc:
            dbg("WARN", "mpl", f"pdb={pdb_id} backend_detect_failed err={exc}")
        else:
            dbg("DEBUG", "mpl", f"backend={backend}")
        dud_eval_log._BACKEND_LOGGED = True

    # EF curve
    order = np.argsort(y_low)
    y_sorted = y_true[order]
    cum_pos = np.cumsum(y_sorted)
    k = np.arange(1, N + 1)
    ef_curve = (N / max(n_act, 1)) * (cum_pos / k)
    frac = k / N
    plt.figure()
    plt.plot(frac, ef_curve)
    plt.xlabel("Fraction screened")
    plt.ylabel("Enrichment factor (EF)")
    plt.title(f"{pdb_id} - EF curve{title_suffix}")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / "ef_curve.png", dpi=200)
    plt.close()

    # ROC
    plt.figure()
    plt.plot(fpr, tpr, label=f"AUC={rocAUC:.3f}")
    plt.plot([0, 1], [0, 1], "k--", lw=1)
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    plt.title(f"{pdb_id} - ROC{title_suffix}")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / "roc.png", dpi=200)
    plt.close()

    # PR
    base = n_act / N
    dbg(
        "INFO",
        "screen",
        f"pdb={pdb_id} mode={mode_label} drop_nonfinite={drop_nonfinite} missing_name_detected={missing_name_detected} drop_missing_token={dropped_token} drop_nan_best_score={drop_nan_best} kept_ligands={N} actives={n_act} active_fraction={base:.3f} token_regex='{token_regex}'",
    )
    plt.figure()
    plt.plot(recall, precision, label=f"PR-AUC={prAUC:.3f}")
    plt.hlines(
        base,
        0,
        1,
        colors="k",
        linestyles="--",
        linewidth=1,
        label=f"Baseline={base:.3f}",
    )
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title(f"{pdb_id} - Precision–Recall{title_suffix}")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / "pr.png", dpi=200)
    plt.close()

    # Score distributions (always plotted as low-is-better values)
    mask_active = best["is_active"] == 1
    mask_decoy = best["is_active"] == 0
    plt.figure()
    plt.hist(y_low[mask_active.to_numpy()], bins=40, alpha=1, label="Actives", zorder=2)
    plt.hist(y_low[mask_decoy.to_numpy()], bins=40, alpha=0.5, label="Decoys", zorder=1)
    plt.xlabel(hist_xlabel)
    plt.ylabel("Count")
    plt.title(f"{pdb_id} - Score distributions{title_suffix}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / "score_hist.png", dpi=200)
    plt.close()

    row = {
        "run_id": str(run_id) if run_id else "(none)",
        "pdb_id": pdb_id,
        "N": N,
        "n_actives": n_act,
        "actives_fraction": base,
        "ROC_AUC": rocAUC,
        "PR_AUC": prAUC,
        "logAUC": lAUC,
        "logAUC_adj": lAUC_adj,
        f"BEDROC_alpha_{bedroc_alpha:g}": bed,
        **ef,
    }
    row["status_reason"] = status_reason if status_reason else "ok"
    if variant_label:
        row["variant"] = variant_label
    pd.DataFrame([row]).to_csv(out_dir / "metrics.tsv", sep="\t", index=False)
    dbg(
        "DEBUG",
        "metrics",
        f"pdb={pdb_id} mode={mode_label} ROC_AUC={rocAUC:.3f} PR_AUC={prAUC:.3f} BEDROC={bed:.3f}",
    )
    return TargetEvaluation(
        metrics=pd.Series(row),
        ligand_basenames=ligand_basenames,
        has_run_id_column=has_run_id_col,
        run_ids=run_ids_present,
        status_reason=str(row["status_reason"]),
        z_selected=z_selected,
    )

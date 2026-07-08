from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.labels import truthy_series
from analysis.ml.decoy_bias_audit import audit_decoy_bias
from analysis.statistics import ranked_binary_metrics


DEFAULT_SCORE_COLUMNS = (
    "z_selected",
    "final_score",
    "SCORCH_score_used",
    "z_vs_decoys_blend",
    "consensus_score",
)


def _score_direction(direction: str | None) -> int:
    return -1 if str(direction or "").strip().lower() in {"lower", "lower_is_better", "min"} else 1


def _prepare_decoy_labels(
    df: pd.DataFrame,
    *,
    positive_col: str,
    decoy_col: str,
    include_fda_background: bool,
) -> pd.DataFrame:
    if positive_col not in df.columns:
        raise ValueError(f"positive column {positive_col!r} not found")
    if decoy_col not in df.columns:
        raise ValueError(f"decoy column {decoy_col!r} not found")
    positive = truthy_series(df[positive_col])
    decoy = truthy_series(df[decoy_col])
    if include_fda_background:
        keep = positive | ~decoy
        negative = ~positive
    else:
        keep = positive | decoy
        negative = decoy
    out = df.loc[keep].copy()
    out["decoy_benchmark_label"] = 0
    out.loc[positive.loc[keep], "decoy_benchmark_label"] = 1
    out["decoy_benchmark_role"] = "negative_decoy"
    out.loc[positive.loc[keep], "decoy_benchmark_role"] = "positive_control"
    if include_fda_background:
        out.loc[negative.loc[keep] & ~decoy.loc[keep], "decoy_benchmark_role"] = "negative_fda_background"
    return out


def _metric_row(
    frame: pd.DataFrame,
    *,
    score_col: str,
    score_multiplier: int,
    group_name: str,
    group_value: str,
) -> dict[str, Any]:
    work = frame[[score_col, "decoy_benchmark_label"]].copy()
    work[score_col] = pd.to_numeric(work[score_col], errors="coerce")
    work = work.dropna(subset=[score_col])
    labels = work["decoy_benchmark_label"].astype(int).tolist()
    scores = (work[score_col] * score_multiplier).astype(float).tolist()
    positives = int(sum(labels))
    negatives = int(len(labels) - positives)
    row: dict[str, Any] = {
        "group": group_name,
        "group_value": group_value,
        "score_col": score_col,
        "n_rows": int(len(labels)),
        "n_positive_controls": positives,
        "n_negative_decoys_or_background": negatives,
        "positive_rate": float(positives / len(labels)) if labels else 0.0,
        "metric_status": "ok" if positives and negatives else "insufficient_classes",
    }
    if positives and negatives:
        row.update(ranked_binary_metrics(scores, labels, [0.01, 0.05, 0.10]))
    else:
        row.update({"EF@1%": pd.NA, "EF@5%": pd.NA, "EF@10%": pd.NA, "AUPRC": pd.NA, "AUROC": pd.NA})
    return row


def _write_score_distribution(rows: pd.DataFrame, score_cols: list[str], out_dir: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return
    available = [col for col in score_cols if col in rows.columns and pd.to_numeric(rows[col], errors="coerce").notna().any()]
    if not available:
        return
    n = len(available)
    fig, axes = plt.subplots(n, 1, figsize=(7, max(3, 2.5 * n)), squeeze=False)
    for ax, score_col in zip(axes[:, 0], available):
        for role, group in rows.groupby("decoy_benchmark_role"):
            values = pd.to_numeric(group[score_col], errors="coerce").dropna()
            if values.empty:
                continue
            ax.hist(values, bins=40, alpha=0.5, label=str(role), density=True)
        ax.set_title(score_col)
        ax.set_xlabel("score")
        ax.set_ylabel("density")
        ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_dir / "decoy_benchmark_score_distributions.png", dpi=160)
    plt.close(fig)


def run_decoy_benchmark(
    input_path: str | Path,
    out_dir: str | Path,
    *,
    score_cols: list[str] | None = None,
    score_directions: dict[str, str] | None = None,
    positive_col: str = "is_control",
    decoy_col: str = "is_decoy",
    group_col: str = "pdb_id",
    include_fda_background: bool = False,
) -> dict[str, pd.DataFrame]:
    df = pd.read_csv(input_path, low_memory=False)
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    score_cols = [col for col in (score_cols or list(DEFAULT_SCORE_COLUMNS)) if col in df.columns]
    if not score_cols:
        raise ValueError("No requested score columns are present in input table")
    score_directions = score_directions or {}
    labeled = _prepare_decoy_labels(
        df,
        positive_col=positive_col,
        decoy_col=decoy_col,
        include_fda_background=include_fda_background,
    )
    labeled.to_csv(out_path / "decoy_benchmark_rows.csv", index=False)
    summary_rows = []
    target_rows = []
    for score_col in score_cols:
        multiplier = _score_direction(score_directions.get(score_col))
        summary_rows.append(
            _metric_row(
                labeled,
                score_col=score_col,
                score_multiplier=multiplier,
                group_name="overall",
                group_value="all",
            )
        )
        if group_col in labeled.columns:
            for group_value, group in labeled.groupby(group_col, dropna=False):
                target_rows.append(
                    _metric_row(
                        group,
                        score_col=score_col,
                        score_multiplier=multiplier,
                        group_name=group_col,
                        group_value=str(group_value),
                    )
                )
    summary = pd.DataFrame(summary_rows)
    by_target = pd.DataFrame(target_rows)
    summary.to_csv(out_path / "decoy_benchmark_summary.csv", index=False)
    by_target.to_csv(out_path / "decoy_benchmark_by_target.csv", index=False)
    manifest = {
        "input_path": str(input_path),
        "score_cols": score_cols,
        "positive_col": positive_col,
        "decoy_col": decoy_col,
        "group_col": group_col,
        "include_fda_background": include_fda_background,
        "interpretation": "Decoy benchmark negatives are technical docking controls and are not FDA/off-target ML labels.",
    }
    (out_path / "decoy_benchmark_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    _write_score_distribution(labeled, score_cols, out_path)
    bias_manifest = audit_decoy_bias(out_path / "decoy_benchmark_rows.csv", out_path / "decoy_bias_audit")
    manifest["decoy_bias_audit"] = bias_manifest.get("outputs", str(out_path / "decoy_bias_audit"))
    (out_path / "decoy_benchmark_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return {"summary": summary, "by_target": by_target, "rows": labeled}

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import pandas as pd

from analysis.run_enrichment import PERMUTATION_COLUMNS, SUMMARY_COLUMNS, run_enrichment
from analysis._common import write_csv_rows


DEFAULT_SCORE_COLS = [
    "atlas_score",
    "mmgbsa_score",
    "exposure_plausibility",
    "target_adr_evidence",
    "pathway_evidence",
    "priority_score",
    "mechanism_graph_score",
    "calibrated_activity_probability",
    "ml_prediction_score",
]
DEFAULT_LABEL_COLS = [
    "literature_supported_label",
    "spd_exposure_relevant",
    "toxcast_active",
    "papyrus_active",
    "chembl_active",
]


def run_enrichment_suite(
    pair_table_path: str | Path,
    score_cols: Sequence[str],
    label_cols: Sequence[str],
    top_fractions: Sequence[float],
    n_permutations: int,
    out_dir: str | Path,
    seed: int = 42,
    score_directions: dict[str, str] | None = None,
) -> pd.DataFrame:
    rows = pd.read_csv(pair_table_path).to_dict("records")
    columns = set(rows[0]) if rows else set()
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    summary_rows: list[dict[str, object]] = []
    permutation_rows: list[dict[str, object]] = []
    top_ks = [max(1, int(round(frac * len(rows)))) for frac in top_fractions]
    for score in score_cols:
        if score not in columns:
            summary_rows.append({"score": score, "label": "", "metric": "skipped_missing_score", "observed": ""})
            continue
        for label in label_cols:
            if label not in columns:
                summary_rows.append({"score": score, "label": label, "metric": "skipped_missing_label", "observed": ""})
                continue
            summary, permutations = run_enrichment(
                rows,
                score=score,
                label=label,
                n_permutations=n_permutations,
                n_bootstraps=min(1000, n_permutations),
                top_ks=top_ks,
                seed=seed,
                direction=(score_directions or {}).get(score, "lower" if score == "mmgbsa_score" else "higher"),
            )
            summary_rows.extend(summary)
            permutation_rows.extend(permutations)
    write_csv_rows(out_path / "enrichment_summary.csv", summary_rows, SUMMARY_COLUMNS)
    write_csv_rows(out_path / "permutation_results.csv", permutation_rows, PERMUTATION_COLUMNS)
    return pd.DataFrame(summary_rows)


def run_ablation_analysis(pair_table_path: str | Path, out_path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(pair_table_path)
    from analysis.priority_score import DEFAULT_WEIGHTS, compute_priority_score
    from analysis.statistics import ranked_binary_metrics

    rows: list[dict[str, object]] = []
    score_tables: dict[str, pd.Series] = {}
    if "priority_score" in df:
        score_tables["full_score"] = pd.to_numeric(df["priority_score"], errors="coerce")
    else:
        score_tables["full_score"] = pd.to_numeric(compute_priority_score(df)["priority_score"], errors="coerce")
    for component in ["exposure_plausibility", "tissue_expression", "target_adr_evidence", "pathway_evidence", "mmgbsa_score", "structure_quality"]:
        weights = {key: value for key, value in DEFAULT_WEIGHTS.items() if key != component}
        score_tables[f"minus_{component}"] = pd.to_numeric(compute_priority_score(df, weights=weights)["priority_score"], errors="coerce")
    if "atlas_score" in df:
        score_tables["atlas_score_only"] = pd.to_numeric(df["atlas_score"], errors="coerce")
    if "mechanism_graph_score" in df:
        score_tables["mechanism_graph_only"] = pd.to_numeric(df["mechanism_graph_score"], errors="coerce")

    label_col = next((col for col in ("spd_exposure_relevant", "literature_supported_label", "toxcast_active", "papyrus_active", "chembl_active") if col in df.columns), "")
    for name, scores in score_tables.items():
        row: dict[str, object] = {"ablation": name, "score_col": name, "n_non_missing": int(scores.notna().sum()), "label": label_col}
        if label_col:
            sub = pd.DataFrame({"score": scores, "label": df[label_col]}).dropna()
            if not sub.empty:
                metrics = ranked_binary_metrics(sub["score"].astype(float).tolist(), sub["label"].astype(bool).astype(int).tolist(), [0.01, 0.05, 0.10])
                row.update(metrics)
                row["n_labeled"] = len(sub)
                row["n_positive"] = int(sub["label"].astype(bool).sum())
        rows.append(row)
    out = pd.DataFrame(rows)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    return out

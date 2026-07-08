from __future__ import annotations

from pathlib import Path
from typing import Sequence

import pandas as pd


STATISTICAL_FORMULAS = (
    {
        "title": "Matched-control p-value",
        "formula": "p_matched = (1 + #{control_score >= observed_score}) / (1 + n_controls)",
        "notes": "Controls are matched on configured fields such as protein class and ligand chemotype. For lower-is-better scores, the inequality is reversed.",
        "source": "analysis.controls.matched_controls.run_matched_controls",
    },
)


def _score_direction_higher(direction: str | None, score_col: str) -> bool:
    if direction:
        return direction == "higher"
    return score_col not in {"mmgbsa_score", "vina_score", "docking_score"}


def _candidate_pool(df: pd.DataFrame, row: pd.Series, match_on: Sequence[str]) -> pd.DataFrame:
    pool = df[~((df["drug_id"] == row["drug_id"]) & (df["target_id"] == row["target_id"]))]
    for field in match_on:
        if field in pool.columns and field in row.index and pd.notna(row[field]):
            matched = pool[pool[field] == row[field]]
            if not matched.empty:
                pool = matched
    return pool


def run_matched_controls(
    pair_table_path: str | Path,
    score_col: str,
    top_n: int,
    n_controls_per_pair: int,
    match_on: list[str],
    out_path: str | Path,
    seed: int = 42,
    direction: str | None = None,
) -> pd.DataFrame:
    df = pd.read_csv(pair_table_path)
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    scored = df.dropna(subset=[score_col]).copy()
    higher = _score_direction_higher(direction, score_col)
    top = scored.sort_values(score_col, ascending=not higher).head(top_n)
    controls: list[dict[str, object]] = []
    significance: list[dict[str, object]] = []
    for idx, (_row_idx, row) in enumerate(top.iterrows()):
        pool = _candidate_pool(scored, row, match_on)
        if pool.empty:
            pool = scored[~((scored["drug_id"] == row["drug_id"]) & (scored["target_id"] == row["target_id"]))]
        sampled = pool.sample(n=min(n_controls_per_pair, len(pool)), replace=len(pool) < n_controls_per_pair, random_state=seed + idx) if not pool.empty else pool
        control_scores = sampled[score_col].tolist()
        observed = float(row[score_col])
        extreme = sum(1 for score in control_scores if score >= observed) if higher else sum(1 for score in control_scores if score <= observed)
        p_value = (extreme + 1.0) / (len(control_scores) + 1.0)
        for control_idx, (_ctrl_i, control) in enumerate(sampled.iterrows(), start=1):
            controls.append(
                {
                    "drug_id": row["drug_id"],
                    "target_id": row["target_id"],
                    "control_index": control_idx,
                    "control_drug_id": control.get("drug_id"),
                    "control_target_id": control.get("target_id"),
                    "control_score": control.get(score_col),
                    "match_on": ";".join(match_on),
                }
            )
        significance.append(
            {
                "drug_id": row["drug_id"],
                "target_id": row["target_id"],
                "pdb_id": row.get("pdb_id", ""),
                "score_col": score_col,
                "observed_score": observed,
                "n_controls": len(control_scores),
                "matched_empirical_p": p_value,
                "match_on": ";".join(match_on),
            }
        )
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    sig_df = pd.DataFrame(significance)
    sig_df.to_csv(out, index=False)
    pd.DataFrame(controls).to_csv(out.with_name("matched_controls.csv"), index=False)
    return sig_df


def write_matched_control_plot(significance_path: str | Path, controls_path: str | Path, out_path: str | Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sig = pd.read_csv(significance_path)
    controls = pd.read_csv(controls_path)
    fig, ax = plt.subplots(figsize=(6, 4), dpi=160)
    if "control_score" in controls:
        pd.to_numeric(controls["control_score"], errors="coerce").dropna().plot.hist(ax=ax, bins=40, alpha=0.7, label="controls")
    if "observed_score" in sig:
        ax.axvline(pd.to_numeric(sig["observed_score"], errors="coerce").median(), color="black", label="median observed")
    ax.set_xlabel("score")
    ax.legend()
    fig.tight_layout()
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)

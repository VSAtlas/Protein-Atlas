from __future__ import annotations

from pathlib import Path

import pandas as pd


def compare_pbas_atlas_profiles(
    pbas_pairs_path: str | Path,
    atlas_pairs_path: str | Path,
    out_dir: str | Path,
    *,
    atlas_score_col: str = "atlas_score",
    top_k: int = 20,
) -> dict[str, int]:
    pbas = pd.read_csv(pbas_pairs_path)
    atlas = pd.read_csv(atlas_pairs_path)
    atlas = atlas.rename(columns={atlas_score_col: "atlas_score"})
    common_cols = ["drug_id", "target_id"]
    merged = pbas.merge(atlas, on=common_cols, how="inner", suffixes=("_pbas", "_atlas"))
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    summary = {
        "pbas_drugs": int(pbas["drug_id"].nunique()),
        "atlas_drugs": int(atlas["drug_id"].nunique()),
        "drug_overlap": int(len(set(pbas["drug_id"].astype(str)) & set(atlas["drug_id"].astype(str)))),
        "pbas_targets": int(pbas["target_id"].nunique()),
        "atlas_targets": int(atlas["target_id"].nunique()),
        "target_overlap": int(len(set(pbas["target_id"].astype(str)) & set(atlas["target_id"].astype(str)))),
        "pair_overlap": len(merged),
    }
    pd.DataFrame([summary]).to_csv(out / "overlap_summary.csv", index=False)
    corr = pd.DataFrame(
        [
            {
                "score_x": "pbas_score",
                "score_y": "atlas_score",
                "n_overlap": len(merged),
                "pearson": pd.to_numeric(merged.get("pbas_score"), errors="coerce").corr(pd.to_numeric(merged.get("atlas_score"), errors="coerce")),
                "spearman": pd.to_numeric(merged.get("pbas_score"), errors="coerce").corr(pd.to_numeric(merged.get("atlas_score"), errors="coerce"), method="spearman"),
            }
        ]
    )
    corr.to_csv(out / "score_correlation.csv", index=False)
    _topk_overlap(pbas, atlas, "drug_id", "target_id", top_k).to_csv(out / "topk_overlap_by_drug.csv", index=False)
    _topk_overlap(pbas, atlas, "target_id", "drug_id", top_k).to_csv(out / "topk_overlap_by_target.csv", index=False)
    _hit_sets(pbas, atlas, merged, top_k, out)
    _scatter(merged, out / "pbas_vs_atlas_score_scatter.png")
    return summary


def _topk_overlap(pbas: pd.DataFrame, atlas: pd.DataFrame, group_col: str, item_col: str, top_k: int) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for group in sorted(set(pbas[group_col].astype(str)) & set(atlas[group_col].astype(str))):
        p_items = set(pbas[pbas[group_col].astype(str) == group].nlargest(top_k, "pbas_score")[item_col].astype(str))
        a_items = set(atlas[atlas[group_col].astype(str) == group].nlargest(top_k, "atlas_score")[item_col].astype(str))
        union = p_items | a_items
        rows.append({group_col: group, "pbas_topk": len(p_items), "atlas_topk": len(a_items), "overlap": len(p_items & a_items), "jaccard": len(p_items & a_items) / len(union) if union else 0.0})
    return pd.DataFrame(rows)


def _hit_sets(pbas: pd.DataFrame, atlas: pd.DataFrame, merged: pd.DataFrame, top_k: int, out: Path) -> None:
    p_hi = pbas.nlargest(top_k, "pbas_score")[["drug_id", "target_id", "pbas_score"]]
    a_hi = atlas.nlargest(top_k, "atlas_score")[["drug_id", "target_id", "atlas_score"]]
    concordant = merged.nlargest(top_k, "atlas_score")
    pbas_only = p_hi.merge(a_hi, on=["drug_id", "target_id"], how="left", indicator=True)
    atlas_only = a_hi.merge(p_hi, on=["drug_id", "target_id"], how="left", indicator=True)
    pbas_only[pbas_only["_merge"] == "left_only"].drop(columns=["_merge"]).to_csv(out / "pbas_only_high_scoring_pairs.csv", index=False)
    atlas_only[atlas_only["_merge"] == "left_only"].drop(columns=["_merge"]).to_csv(out / "atlas_only_high_scoring_pairs.csv", index=False)
    concordant.to_csv(out / "concordant_hits.csv", index=False)
    pd.concat([pbas_only[pbas_only["_merge"] == "left_only"], atlas_only[atlas_only["_merge"] == "left_only"]], ignore_index=True).to_csv(out / "discordant_hits.csv", index=False)


def _scatter(df: pd.DataFrame, out_path: Path) -> None:
    if df.empty:
        return
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5.5, 4.2), dpi=160)
    ax.scatter(pd.to_numeric(df["pbas_score"], errors="coerce"), pd.to_numeric(df["atlas_score"], errors="coerce"), s=10, alpha=0.45)
    ax.set_xlabel("PBAS score")
    ax.set_ylabel("Atlas score")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)

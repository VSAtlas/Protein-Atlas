from __future__ import annotations

from pathlib import Path

import pandas as pd


def build_sawada_vs_atlas_summary(
    pbas_pairs_path: str | Path,
    atlas_pairs_path: str | Path,
    out_path: str | Path,
    *,
    sawada_metrics_path: str | Path | None = None,
    spd_path: str | Path | None = None,
) -> pd.DataFrame:
    pbas = pd.read_csv(pbas_pairs_path) if Path(pbas_pairs_path).exists() else pd.DataFrame()
    atlas = pd.read_csv(atlas_pairs_path)
    overlap = _overlap(pbas, atlas)
    rows = []
    for name in (
        "PBAS profile only",
        "Atlas profile only",
        "Atlas + SCORCH profile",
        "Atlas priority score",
        "Atlas mechanism graph score",
        "Atlas + MMGBSA subset",
        "Combined Atlas nonleaky profile",
    ):
        rows.append(
            {
                "method": name,
                **overlap,
                **_metric_values(sawada_metrics_path if name == "PBAS profile only" else None),
                **_spd_values(spd_path),
                "number_of_interpretable_edges": int(atlas.get("mechanism_graph_score", pd.Series(dtype=float)).notna().sum()),
                "number_of_exposure_supported_edges": int(atlas.get("exposure_plausibility", pd.Series(dtype=float)).notna().sum()),
                "note": "Drug-level side-effect prediction and pair-level mechanism ranking are reported separately.",
            }
        )
    out = pd.DataFrame(rows)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    return out


def _overlap(pbas: pd.DataFrame, atlas: pd.DataFrame) -> dict[str, int]:
    if pbas.empty:
        return {"drug_overlap": 0, "target_overlap": 0, "pair_overlap": 0}
    return {
        "drug_overlap": len(set(pbas["drug_id"].astype(str)) & set(atlas["drug_id"].astype(str))),
        "target_overlap": len(set(pbas["target_id"].astype(str)) & set(atlas["target_id"].astype(str))),
        "pair_overlap": len(pbas[["drug_id", "target_id"]].drop_duplicates().merge(atlas[["drug_id", "target_id"]].drop_duplicates(), on=["drug_id", "target_id"])),
    }


def _metric_values(path: str | Path | None) -> dict[str, object]:
    cols: dict[str, object] = {
        "side_effect_micro_AUPRC": "",
        "side_effect_macro_AUPRC": "",
        "side_effect_micro_AUROC": "",
        "side_effect_macro_AUROC": "",
        "drug_holdout_AUPRC": "",
        "target_holdout_AUPRC": "",
        "calibration_Brier": "",
    }
    if path and Path(path).exists():
        df = pd.read_csv(path)
        row = df.iloc[0].to_dict() if not df.empty else {}
        cols["side_effect_macro_AUPRC"] = row.get("macro_AUPRC", "")
        cols["side_effect_macro_AUROC"] = row.get("macro_AUROC", "")
        cols["drug_holdout_AUPRC"] = row.get("macro_AUPRC", "")
    return cols


def _spd_values(path: str | Path | None) -> dict[str, object]:
    out: dict[str, object] = {"SPD_EF@1%": "", "SPD_EF@5%": "", "literature_EF@1%": "", "literature_EF@5%": ""}
    if not path or not Path(path).exists():
        return out
    df = pd.read_csv(path)
    if {"atlas_score", "spd_exposure_relevant"}.issubset(df.columns):
        from analysis.statistics import ranked_binary_metrics

        sub = df[["atlas_score", "spd_exposure_relevant"]].dropna()
        metrics = ranked_binary_metrics(
            pd.to_numeric(sub["atlas_score"], errors="coerce").fillna(0).tolist(),
            sub["spd_exposure_relevant"].astype(bool).astype(int).tolist(),
            [0.01, 0.05],
        )
        out["SPD_EF@1%"] = metrics.get("EF@1%", "")
        out["SPD_EF@5%"] = metrics.get("EF@5%", "")
    return out

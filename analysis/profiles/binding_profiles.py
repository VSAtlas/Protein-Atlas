from __future__ import annotations

from pathlib import Path

import pandas as pd


def build_binding_profile(
    pair_table_path: str | Path,
    out_path: str | Path,
    *,
    value_col: str,
    target_space_path: str | Path | None = None,
    dense: bool = False,
    protein_class: str | None = None,
    min_structure_quality: float | None = None,
    pdb_panel: list[str] | None = None,
) -> pd.DataFrame:
    df = pd.read_csv(pair_table_path)
    if protein_class and "protein_class" in df:
        df = df[df["protein_class"].astype(str) == protein_class]
    if min_structure_quality is not None and "structure_quality" in df:
        df = df[pd.to_numeric(df["structure_quality"], errors="coerce") >= min_structure_quality]
    if pdb_panel and "pdb_id" in df:
        df = df[df["pdb_id"].astype(str).isin({str(v) for v in pdb_panel})]
    required = {"drug_id", "target_id", value_col}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"missing required profile columns: {sorted(missing)}")
    scored = df[["drug_id", "target_id", value_col]].copy()
    scored[value_col] = pd.to_numeric(scored[value_col], errors="coerce")
    scored = scored.sort_values(value_col, ascending=False).drop_duplicates(["drug_id", "target_id"])
    profile = scored.pivot(index="drug_id", columns="target_id", values=value_col)
    if target_space_path is not None and Path(target_space_path).exists():
        target_space = pd.read_csv(target_space_path)
        target_col = "target_id" if "target_id" in target_space else target_space.columns[0]
        profile = profile.reindex(columns=[str(v) for v in target_space[target_col].dropna().astype(str).unique()])
    if dense:
        profile = profile.fillna(0.0)
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.suffix.lower() == ".parquet":
        profile.to_parquet(out)
    else:
        profile.to_csv(out)
    return profile

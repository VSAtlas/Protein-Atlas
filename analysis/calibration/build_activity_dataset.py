from __future__ import annotations

from pathlib import Path

import pandas as pd


def build_activity_dataset(benchmark_table_path: str | Path, label_col: str, score_col: str, out_path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(benchmark_table_path)
    out = df[[col for col in [score_col, label_col, "drug_id", "target_id", "protein_class", "ligand_chemotype"] if col in df.columns]].copy()
    out = out.dropna(subset=[score_col, label_col])
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    return out


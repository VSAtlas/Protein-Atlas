from __future__ import annotations

from pathlib import Path

import pandas as pd

from analysis.external.mapping import normalize_columns


def load_sider(path: str | Path) -> pd.DataFrame:
    sep = "\t" if Path(path).suffix.lower() in {".tsv", ".tab"} else ","
    return normalize_columns(
        pd.read_csv(path, sep=sep),
        {"drug_id": ["drug_id", "stitch_id", "drug"], "adr": ["adr", "side_effect", "meddra_term"], "source": ["source"]},
    )


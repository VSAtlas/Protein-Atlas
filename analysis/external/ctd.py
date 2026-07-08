from __future__ import annotations

from pathlib import Path

import pandas as pd

from analysis.external.mapping import normalize_columns


def load_ctd(path: str | Path) -> pd.DataFrame:
    sep = "\t" if Path(path).suffix.lower() in {".tsv", ".tab"} else ","
    return normalize_columns(
        pd.read_csv(path, sep=sep, low_memory=False),
        {
            "chemical_id": ["chemical_id", "chemical", "drug_id"],
            "gene_id": ["gene_id", "gene", "target_id"],
            "disease_id": ["disease_id", "disease"],
            "edge_type": ["edge_type", "relationship"],
            "pubmed_ids": ["pubmed_ids", "pmids"],
            "source": ["source"],
        },
    )

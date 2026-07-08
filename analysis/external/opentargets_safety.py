from __future__ import annotations

from pathlib import Path

import pandas as pd

from analysis.external.mapping import normalize_columns


def load_opentargets_safety(path: str | Path) -> pd.DataFrame:
    sep = "\t" if Path(path).suffix.lower() in {".tsv", ".tab"} else ","
    return normalize_columns(
        pd.read_csv(path, sep=sep),
        {
            "target_id": ["target_id", "uniprot", "gene", "gene_symbol"],
            "adr": ["adr", "event", "toxicity", "safety_event"],
            "confidence": ["confidence", "score"],
            "pubmed_ids": ["pubmed_ids", "pmids"],
            "source": ["source"],
        },
    )


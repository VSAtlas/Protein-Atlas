from __future__ import annotations

from pathlib import Path

import pandas as pd

from analysis.external.mapping import normalize_columns


def load_reactome(path: str | Path) -> pd.DataFrame:
    sep = "\t" if Path(path).suffix.lower() in {".tsv", ".tab"} else ","
    return normalize_columns(
        pd.read_csv(path, sep=sep),
        {
            "target_id": ["target_id", "uniprot", "gene", "gene_symbol"],
            "pathway": ["pathway", "pathway_name", "reactome_pathway"],
            "pathway_id": ["pathway_id", "reactome_id"],
            "adr": ["adr", "linked_adr", "toxicity"],
            "source": ["source"],
        },
    )


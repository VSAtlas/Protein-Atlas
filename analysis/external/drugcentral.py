from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from analysis.external.compound_index import load_atlas_compound_index


def _clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def _name_key(value: Any) -> str:
    return _clean(value).lower()


def _atlas_name_map(mapping_path: str | Path) -> pd.DataFrame:
    atlas = load_atlas_compound_index(mapping_path)
    rows: list[dict[str, Any]] = []
    for row in atlas.to_dict("records"):
        for field in ("drug_name",):
            name = _clean(row.get(field))
            if name:
                rows.append({"drug_name_key": _name_key(name), "drug_id": row["drug_id"], "atlas_drug_name": name})
    if not rows:
        return pd.DataFrame(columns=["drug_name_key", "drug_id", "atlas_drug_name"])
    return pd.DataFrame(rows).drop_duplicates("drug_name_key")


def _activity_nm_from_pactivity(value: pd.Series) -> pd.Series:
    pvalue = pd.to_numeric(value, errors="coerce")
    return 10 ** (9.0 - pvalue)


def stage_drugcentral_activity(
    interactions_path: str | Path,
    mapping_path: str | Path,
    out_path: str | Path,
) -> pd.DataFrame:
    """Stage DrugCentral target interactions as secondary activity evidence.

    DrugCentral ACT_VALUE is a pActivity-like value in the local file, not a raw nM
    concentration. We preserve ACT_SOURCE lineage and convert to nM only for the
    generic four-state activity threshold logic.
    """

    raw = pd.read_csv(interactions_path, sep="\t", low_memory=False)
    human = raw[raw.get("ORGANISM", "").astype(str).str.contains("Homo sapiens", case=False, na=False)].copy()
    measured = human[human["ACT_VALUE"].notna() & human["ACCESSION"].notna()].copy()
    measured["drug_name_key"] = measured["DRUG_NAME"].map(_name_key)
    measured = measured.merge(_atlas_name_map(mapping_path), on="drug_name_key", how="left")
    measured["drug_id"] = measured["drug_id"].fillna("drugcentral:" + measured["STRUCT_ID"].astype(str))
    measured["activity_nM"] = _activity_nm_from_pactivity(measured["ACT_VALUE"])
    relation = measured["RELATION"].fillna("=")
    relation = relation.mask(relation.astype(str).eq("~"), "=")
    relation = relation.mask(~relation.astype(str).isin(["=", "<", "<=", ">", ">="]), pd.NA)
    measured["activity_relation"] = relation
    measured["assay_id"] = (
        "DrugCentral:"
        + measured["STRUCT_ID"].astype(str)
        + ":"
        + measured["ACCESSION"].astype(str)
        + ":"
        + measured["ACT_TYPE"].astype(str)
    )
    out = pd.DataFrame(
        {
            "drug_id": measured["drug_id"],
            "drug_name": measured["atlas_drug_name"].fillna(measured["DRUG_NAME"]),
            "target_id": measured["ACCESSION"],
            "gene_symbol": measured["GENE"],
            "uniprot": measured["ACCESSION"],
            "assay_id": measured["assay_id"],
            "activity_nM": measured["activity_nM"].replace([np.inf, -np.inf], pd.NA),
            "activity_relation": measured["activity_relation"],
            "activity_type": measured["ACT_TYPE"],
            "activity_units": "nM",
            "pchembl_value": pd.to_numeric(measured["ACT_VALUE"], errors="coerce"),
            "source": "DrugCentral:" + measured["ACT_SOURCE"].fillna("unknown").astype(str),
            "activity_document_ids": measured["ACT_SOURCE_URL"],
            "data_validity_comment": measured["ACT_COMMENT"],
            "drugcentral_struct_id": measured["STRUCT_ID"],
            "drugcentral_source": measured["ACT_SOURCE"],
        }
    )
    out_file = Path(out_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_file, sep="\t", index=False)
    manifest = {
        "source": "DrugCentral target interactions",
        "rows": int(len(out)),
        "atlas_drug_mapped_rows": int((~out["drug_id"].astype(str).str.startswith("drugcentral:")).sum()),
        "human_measured_input_rows": int(len(measured)),
        "output": str(out_file),
        "warning": "DrugCentral is a derived wrapper source; source-family holdout is required for claims.",
    }
    out_file.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return out

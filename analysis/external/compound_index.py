from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd


def atlas_drug_id(row: dict[str, Any] | pd.Series) -> str:
    existing = row.get("drug_id")
    if pd.notna(existing) and str(existing).strip():
        return str(existing).strip()
    scheme = str(row.get("scheme") or "drug").strip()
    file_num = row.get("file_num")
    try:
        if file_num is None or pd.isna(file_num):
            raise ValueError("missing file_num")
        return f"{scheme}_{int(str(file_num)):07d}"
    except Exception:
        fallback = row.get("display_name") or row.get("generic_name") or row.get("pubchem_cid_resolved")
        return str(fallback).strip()


def load_atlas_compound_index(mapping_path: str | Path) -> pd.DataFrame:
    mapping = pd.read_csv(mapping_path, low_memory=False)
    rows: list[dict[str, Any]] = []
    for row in mapping.to_dict("records"):
        pubchem_cid = pd.to_numeric(pd.Series([row.get("pubchem_cid_resolved") or row.get("pubchem_cid")]), errors="coerce").iloc[0]
        rows.append(
            {
                "drug_id": atlas_drug_id(row),
                "drug_name": row.get("display_name") or row.get("generic_name") or row.get("pubchem_name"),
                "inchikey": _clean_identifier(row.get("inchikey") or row.get("remark_inchikey")),
                "pubchem_cid": int(pubchem_cid) if pd.notna(pubchem_cid) else pd.NA,
                "smiles": row.get("smiles") or row.get("remark_smiles"),
                "scheme": row.get("scheme"),
                "file_num": row.get("file_num"),
            }
        )
    out = pd.DataFrame(rows).drop_duplicates()
    if not out.empty:
        out["inchikey_connectivity"] = out["inchikey"].map(inchikey_connectivity)
    return out


def inchikey_connectivity(value: Any) -> str:
    text = _clean_identifier(value)
    if not text:
        return ""
    return text.split("-", 1)[0]


def _clean_identifier(value: Any) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def map_by_inchikey(source: pd.DataFrame, atlas_index: pd.DataFrame, *, source_inchikey_col: str = "inchikey") -> pd.DataFrame:
    out = source.copy()
    if source_inchikey_col not in out.columns or atlas_index.empty:
        out["drug_id"] = pd.NA
        out["drug_name"] = pd.NA
        return out
    out = out.rename(
        columns={
            col: f"source_{col}"
            for col in ("drug_id", "drug_name", "atlas_smiles", "atlas_pubchem_cid")
            if col in out.columns
        }
    )
    out["_inchikey_full"] = out[source_inchikey_col].map(_clean_identifier)
    out["_inchikey_connectivity"] = out["_inchikey_full"].map(inchikey_connectivity)

    full_map = (
        atlas_index[atlas_index["inchikey"].astype(str).ne("")]
        .drop_duplicates("inchikey")[["inchikey", "drug_id", "drug_name", "smiles", "pubchem_cid"]]
        .rename(columns={"inchikey": "_inchikey_full", "smiles": "atlas_smiles", "pubchem_cid": "atlas_pubchem_cid"})
    )
    out = out.merge(full_map, on="_inchikey_full", how="left")

    missing = out["drug_id"].isna()
    if missing.any():
        conn_map = (
            atlas_index[atlas_index["inchikey_connectivity"].astype(str).ne("")]
            .drop_duplicates("inchikey_connectivity")[["inchikey_connectivity", "drug_id", "drug_name", "smiles", "pubchem_cid"]]
            .rename(
                columns={
                    "inchikey_connectivity": "_inchikey_connectivity",
                    "drug_id": "drug_id_connectivity",
                    "drug_name": "drug_name_connectivity",
                    "smiles": "atlas_smiles_connectivity",
                    "pubchem_cid": "atlas_pubchem_cid_connectivity",
                }
            )
        )
        out = out.merge(conn_map, on="_inchikey_connectivity", how="left")
        out["drug_id"] = out["drug_id"].astype("object").combine_first(out["drug_id_connectivity"])
        out["drug_name"] = out["drug_name"].astype("object").combine_first(out["drug_name_connectivity"])
        out["atlas_smiles"] = out["atlas_smiles"].astype("object").combine_first(out["atlas_smiles_connectivity"])
        out["atlas_pubchem_cid"] = out["atlas_pubchem_cid"].astype("object").combine_first(
            out["atlas_pubchem_cid_connectivity"]
        )
        out = out.drop(
            columns=[
                "drug_id_connectivity",
                "drug_name_connectivity",
                "atlas_smiles_connectivity",
                "atlas_pubchem_cid_connectivity",
            ],
            errors="ignore",
        )

    return out.drop(columns=["_inchikey_full", "_inchikey_connectivity"], errors="ignore")


def map_by_pubchem_cid(source: pd.DataFrame, atlas_index: pd.DataFrame, *, source_cid_col: str = "pubchem_cid") -> pd.DataFrame:
    out = source.copy()
    if source_cid_col not in out.columns or atlas_index.empty:
        out["drug_id"] = pd.NA
        out["drug_name"] = pd.NA
        return out
    out["_pubchem_cid"] = pd.to_numeric(out[source_cid_col], errors="coerce").astype("Int64")
    cid_map = (
        atlas_index[atlas_index["pubchem_cid"].notna()]
        .drop_duplicates("pubchem_cid")[["pubchem_cid", "drug_id", "drug_name", "inchikey", "smiles"]]
        .rename(columns={"pubchem_cid": "_pubchem_cid", "inchikey": "atlas_inchikey", "smiles": "atlas_smiles"})
    )
    return out.merge(cid_map, on="_pubchem_cid", how="left").drop(columns=["_pubchem_cid"], errors="ignore")

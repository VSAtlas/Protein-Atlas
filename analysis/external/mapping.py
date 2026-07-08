from __future__ import annotations

import re
from pathlib import Path
from typing import Mapping

import pandas as pd


def read_optional_mapping(path: str | Path | None) -> pd.DataFrame:
    if path is None or str(path).strip() == "":
        return pd.DataFrame()
    map_path = Path(path)
    if not map_path.exists():
        return pd.DataFrame()
    sep = "\t" if map_path.suffix.lower() in {".tsv", ".tab"} else ","
    return pd.read_csv(map_path, sep=sep)


def _normalize_alias_key(value: object) -> str:
    token = str(value or "").strip().lower()
    token = re.sub(r"[^a-z0-9]+", "_", token)
    token = token.strip("_")
    return token


def normalize_columns(df: pd.DataFrame, aliases: Mapping[str, list[str]]) -> pd.DataFrame:
    out = df.copy()
    normalized_to_original = {_normalize_alias_key(col): col for col in out.columns}
    for canonical, names in aliases.items():
        if canonical in out.columns:
            continue
        for name in names:
            original = normalized_to_original.get(_normalize_alias_key(name))
            if original is not None:
                out[canonical] = out[original]
                break
        if canonical not in out.columns:
            out[canonical] = pd.NA
    return out


def apply_pair_mapping(df: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    if mapping.empty:
        return df
    out = df.copy()
    aliases = normalize_columns(
        mapping,
        {
            "source_drug_id": ["source_drug_id", "external_drug_id", "drug_name", "compound_id"],
            "source_target_id": ["source_target_id", "external_target_id", "assay_target", "gene_symbol"],
            "drug_id": ["drug_id", "atlas_drug_id"],
            "target_id": ["target_id", "atlas_target_id", "uniprot"],
        },
    )
    if {"source_drug_id", "drug_id"}.issubset(aliases.columns) and "drug_id" in out.columns:
        drug_map = aliases.dropna(subset=["source_drug_id", "drug_id"]).drop_duplicates("source_drug_id")
        out = out.merge(drug_map[["source_drug_id", "drug_id"]], left_on="drug_id", right_on="source_drug_id", how="left", suffixes=("", "_mapped"))
        out["drug_id"] = out["drug_id_mapped"].fillna(out["drug_id"])
        out = out.drop(columns=[col for col in ("source_drug_id", "drug_id_mapped") if col in out.columns])
    if {"source_target_id", "target_id"}.issubset(aliases.columns) and "target_id" in out.columns:
        target_map = aliases.dropna(subset=["source_target_id", "target_id"]).drop_duplicates("source_target_id")
        out = out.merge(target_map[["source_target_id", "target_id"]], left_on="target_id", right_on="source_target_id", how="left", suffixes=("", "_mapped"))
        out["target_id"] = out["target_id_mapped"].fillna(out["target_id"])
        out = out.drop(columns=[col for col in ("source_target_id", "target_id_mapped") if col in out.columns])
    return out

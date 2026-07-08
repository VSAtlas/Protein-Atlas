from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from analysis.external.compound_index import load_atlas_compound_index


def load_adrecs_target(path: str | Path) -> pd.DataFrame:
    """Load ADReCS-Target ADR-protein associations as target safety evidence."""

    source = Path(path)
    if source.suffix.lower() in {".tsv", ".txt"}:
        raw = pd.read_csv(source, sep="\t", low_memory=False)
    elif source.suffix.lower() == ".csv":
        raw = pd.read_csv(source, low_memory=False)
    else:
        raw = pd.read_excel(source)
    lower = {str(col).strip().lower(): col for col in raw.columns}

    def pick(*names: str) -> pd.Series:
        for name in names:
            col = lower.get(name.lower())
            if col is not None:
                return raw[col]
        return pd.Series(pd.NA, index=raw.index)

    out = pd.DataFrame(
        {
            "target_id": pick("Uniprot AC", "uniprot", "uniprot_ac"),
            "adr": pick("ADR Term", "adr", "adr_term"),
            "confidence": 1.0,
            "pubmed_ids": "",
            "source": "ADReCS-Target",
            "_source_badd_tid": pick("BADD_TID"),
            "_source_adr_id": pick("ADR_ID"),
            "_source_adrecs_id": pick("ADReCS ID"),
            "_source_drug_name": pick("Drug_Name"),
        }
    )
    out["target_id"] = out["target_id"].astype(str).str.strip()
    out["adr"] = out["adr"].astype(str).str.strip()
    out = out[out["target_id"].ne("") & out["adr"].ne("") & out["target_id"].ne("nan") & out["adr"].ne("nan")]
    return out.drop_duplicates()


def stage_adrecs_target(
    association_path: str | Path,
    out_path: str | Path,
    *,
    protein_info_path: str | Path | None = None,
    adr_info_path: str | Path | None = None,
) -> pd.DataFrame:
    """Normalize ADReCS-Target associations as target-ADR mechanism positives."""

    assoc = load_adrecs_target(association_path)
    assoc = assoc.rename(columns={"adr": "adr_id"})
    assoc["adr_term"] = assoc["adr_id"]
    assoc["source_name"] = "ADReCS-Target"
    assoc["evidence_namespace"] = "target_adr_mechanism"
    assoc["label_state"] = 1
    assoc["label_confidence"] = assoc["confidence"].fillna(1.0)
    assoc["training_allowed"] = True
    assoc["benchmark_only"] = False
    assoc["source_role"] = "curated_target_adr_positive"
    assoc["provenance_priority"] = 2

    if protein_info_path is not None and Path(protein_info_path).exists():
        protein = pd.read_excel(protein_info_path)
        protein = protein.rename(
            columns={
                "UNIPROT_AC": "target_id",
                "Gene.names": "gene_symbol",
                "GeneID": "entrez_gene_id",
                "Protein.names": "target_name",
            }
        )
        keep = [col for col in ["target_id", "gene_symbol", "entrez_gene_id", "target_name"] if col in protein.columns]
        assoc = assoc.merge(protein[keep].drop_duplicates("target_id"), on="target_id", how="left")

    if adr_info_path is not None and Path(adr_info_path).exists():
        adr = pd.read_excel(adr_info_path)
        adr = adr.rename(
            columns={
                "ADR_ID": "_source_adr_id",
                "ADRECS_ID": "adrecs_id",
                "ADR_TERM": "adrecs_adr_term",
                "DATA_SOURCE": "pubmed_ids",
            }
        )
        keep = [col for col in ["_source_adr_id", "adrecs_id", "adrecs_adr_term", "pubmed_ids"] if col in adr.columns]
        assoc = assoc.merge(adr[keep].drop_duplicates("_source_adr_id"), on="_source_adr_id", how="left")
        assoc["adr_term"] = assoc["adrecs_adr_term"].fillna(assoc["adr_term"])
    elif "pubmed_ids" not in assoc.columns:
        assoc["pubmed_ids"] = ""

    columns = [
        "target_id",
        "gene_symbol",
        "entrez_gene_id",
        "target_name",
        "adr_id",
        "adr_term",
        "source_name",
        "evidence_namespace",
        "label_state",
        "label_confidence",
        "training_allowed",
        "benchmark_only",
        "source_role",
        "provenance_priority",
        "pubmed_ids",
        "_source_badd_tid",
        "_source_adr_id",
        "_source_adrecs_id",
        "_source_drug_name",
    ]
    for col in columns:
        if col not in assoc.columns:
            assoc[col] = pd.NA
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    normalized = assoc[columns].drop_duplicates()
    normalized.to_csv(out, sep="\t", index=False)
    manifest = {
        "source": "ADReCS-Target",
        "input": str(association_path),
        "output": str(out),
        "rows": int(len(normalized)),
        "unique_targets": int(normalized["target_id"].nunique()),
        "unique_adrs": int(normalized["adr_id"].nunique()),
        "label_policy": "Curated ADReCS target-ADR associations are strict positive mechanism evidence only; absence is unknown.",
    }
    out.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return normalized


def _drug_name_map(mapping_path: str | Path) -> pd.DataFrame:
    mapping = pd.read_csv(mapping_path, low_memory=False)
    atlas = load_atlas_compound_index(mapping_path)
    id_to_drug = atlas.drop_duplicates("drug_id").set_index("drug_id")
    name_cols = [
        "display_name",
        "generic_name",
        "pubchem_name",
        "rxnorm_generic_name",
        "drugcentral_generic_name",
        "sdf_title",
    ]
    rows: list[dict[str, str]] = []
    for idx, row in mapping.iterrows():
        drug_id = str(atlas.iloc[idx]["drug_id"]) if idx < len(atlas) else ""
        for col in name_cols:
            value = row.get(col)
            if pd.isna(value):
                continue
            text = str(value).strip()
            if not text:
                continue
            rows.append(
                {
                    "_drug_name_key": text.lower(),
                    "drug_id": drug_id,
                    "mapped_drug_name": str(id_to_drug.at[drug_id, "drug_name"]) if drug_id in id_to_drug.index else text,
                }
            )
    return pd.DataFrame(rows).drop_duplicates()


def stage_adrecs_drug_target_adr(
    target_adr_path: str | Path,
    mapping_path: str | Path,
    out_path: str | Path,
) -> pd.DataFrame:
    """Map ADReCS drug-target-ADR assertions to Atlas FDA ligand IDs."""

    target_adr = pd.read_csv(target_adr_path, sep="\t", low_memory=False)
    drug_map = _drug_name_map(mapping_path)
    if target_adr.empty or drug_map.empty:
        out = pd.DataFrame()
    else:
        work = target_adr.copy()
        work["_drug_name_key"] = work["_source_drug_name"].fillna("").astype(str).str.strip().str.lower()
        out = work.merge(drug_map, on="_drug_name_key", how="inner")
    if out.empty:
        columns = [
            "drug_id",
            "drug_name",
            "target_id",
            "adr_id",
            "adr_term",
            "evidence_scope",
            "mechanism_label",
            "positive_evidence_type",
            "positive_source",
            "positive_confidence",
            "positive_label_status",
            "document_ids",
        ]
        normalized = pd.DataFrame(columns=columns)
    else:
        normalized = pd.DataFrame(
            {
                "drug_id": out["drug_id"],
                "drug_name": out["mapped_drug_name"],
                "target_id": out["target_id"],
                "adr_id": out["adr_id"],
                "adr_term": out["adr_term"],
                "evidence_scope": "drug_target",
                "mechanism_label": 1,
                "positive_evidence_type": "curated_adrecs_drug_target_adr",
                "positive_source": "ADReCS-Target",
                "positive_confidence": out.get("label_confidence", pd.Series(0.9, index=out.index)),
                "positive_label_status": "strict_positive_drug_target_adr",
                "document_ids": out.get("pubmed_ids", pd.Series(pd.NA, index=out.index)),
                "_source_drug_name": out["_source_drug_name"],
                "_source_badd_tid": out["_source_badd_tid"],
                "_source_adr_id": out["_source_adr_id"],
                "_source_adrecs_id": out["_source_adrecs_id"],
            }
        ).drop_duplicates()
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    normalized.to_csv(path, sep="\t", index=False)
    manifest = {
        "source": "ADReCS-Target",
        "input": str(target_adr_path),
        "mapping": str(mapping_path),
        "output": str(path),
        "rows": int(len(normalized)),
        "unique_drugs": int(normalized["drug_id"].nunique()) if "drug_id" in normalized else 0,
        "unique_targets": int(normalized["target_id"].nunique()) if "target_id" in normalized else 0,
        "unique_adrs": int(normalized["adr_id"].nunique()) if "adr_id" in normalized else 0,
        "label_policy": "ADReCS rows with mapped FDA drug name plus target plus ADR are strict positive mechanism evidence; absence remains unknown.",
    }
    path.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return normalized

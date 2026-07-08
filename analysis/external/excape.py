from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.external.compound_index import inchikey_connectivity, load_atlas_compound_index


EXCAPE_COLUMNS = [
    "Ambit_InchiKey",
    "Original_Entry_ID",
    "Entrez_ID",
    "Activity_Flag",
    "pXC50",
    "DB",
    "Original_Assay_ID",
    "Tax_ID",
    "Gene_Symbol",
    "InChI",
    "SMILES",
]


def _clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def _target_map(pair_table_path: str | Path) -> pd.DataFrame:
    pairs = pd.read_csv(pair_table_path, low_memory=False)
    gene_col = next((col for col in ("target_gene", "gene_symbol", "target_symbol") if col in pairs.columns), None)
    target_col = next((col for col in ("target_id", "target_uniprot", "uniprot") if col in pairs.columns), None)
    if gene_col is None or target_col is None:
        return pd.DataFrame(columns=["gene_symbol", "target_id"])
    out = pairs[[gene_col, target_col]].dropna().drop_duplicates()
    out.columns = ["gene_symbol", "target_id"]
    out["gene_symbol"] = out["gene_symbol"].astype(str).str.upper().str.strip()
    out["target_id"] = out["target_id"].astype(str).str.strip()
    return out.drop_duplicates("gene_symbol")


def _drug_map(mapping_path: str | Path) -> pd.DataFrame:
    atlas = load_atlas_compound_index(mapping_path)
    rows: list[dict[str, Any]] = []
    for row in atlas.to_dict("records"):
        inchikey = _clean(row.get("inchikey"))
        if not inchikey:
            continue
        rows.append(
            {
                "ambit_inchikey": inchikey,
                "inchikey_connectivity": inchikey_connectivity(inchikey),
                "drug_id": row["drug_id"],
                "drug_name": row.get("drug_name"),
            }
        )
    return pd.DataFrame(rows).drop_duplicates("ambit_inchikey")


def _activity_label(flag: Any) -> str:
    text = _clean(flag).upper()
    if text == "A":
        return "active"
    if text == "N":
        return "inactive"
    return "inconclusive"


def stage_excape_for_atlas(
    excape_path: str | Path,
    pair_table_path: str | Path,
    mapping_path: str | Path,
    out_path: str | Path,
    *,
    chunksize: int = 500_000,
) -> pd.DataFrame:
    """Stream-filter ExCAPE/PubChem+ChEMBL to Atlas drugs and targets.

    The source file can be tens of GiB. This function never loads it all into
    memory and only emits rows mapping to the supplied Atlas pair table target
    genes and FDA/Atlas ligand InChIKeys.
    """

    target_map = _target_map(pair_table_path)
    drug_map = _drug_map(mapping_path)
    if target_map.empty or drug_map.empty:
        raise ValueError("ExCAPE staging requires nonempty target and drug maps")
    target_genes = set(target_map["gene_symbol"].dropna().astype(str))
    full_key_map = drug_map.drop_duplicates("ambit_inchikey").set_index("ambit_inchikey")
    conn_key_map = drug_map[drug_map["inchikey_connectivity"].ne("")].drop_duplicates("inchikey_connectivity").set_index(
        "inchikey_connectivity"
    )
    out_file = Path(out_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    first = True
    rows_in = 0
    rows_human_target = 0
    rows_drug_mapped = 0
    reader = pd.read_csv(
        excape_path,
        sep="\t",
        usecols=lambda col: col in EXCAPE_COLUMNS,
        chunksize=chunksize,
        low_memory=False,
    )
    for chunk in reader:
        rows_in += len(chunk)
        chunk["Gene_Symbol"] = chunk["Gene_Symbol"].astype("object").where(chunk["Gene_Symbol"].notna(), "").astype(str).str.upper()
        human_target = chunk[pd.to_numeric(chunk["Tax_ID"], errors="coerce").eq(9606) & chunk["Gene_Symbol"].isin(target_genes)].copy()
        if human_target.empty:
            continue
        rows_human_target += len(human_target)
        human_target["ambit_inchikey"] = human_target["Ambit_InchiKey"].map(_clean)
        human_target["inchikey_connectivity"] = human_target["ambit_inchikey"].map(inchikey_connectivity)
        mapped = human_target.merge(
            full_key_map[["drug_id", "drug_name"]],
            left_on="ambit_inchikey",
            right_index=True,
            how="left",
        )
        missing = mapped["drug_id"].isna()
        if missing.any():
            conn = mapped.loc[missing].merge(
                conn_key_map[["drug_id", "drug_name"]],
                left_on="inchikey_connectivity",
                right_index=True,
                how="left",
                suffixes=("", "_conn"),
            )
            mapped.loc[missing, "drug_id"] = conn["drug_id_conn"].to_numpy()
            mapped.loc[missing, "drug_name"] = conn["drug_name_conn"].to_numpy()
        mapped = mapped[mapped["drug_id"].notna()].copy()
        if mapped.empty:
            continue
        rows_drug_mapped += len(mapped)
        mapped = mapped.merge(target_map, left_on="Gene_Symbol", right_on="gene_symbol", how="left")
        staged = pd.DataFrame(
            {
                "drug_id": mapped["drug_id"],
                "drug_name": mapped["drug_name"],
                "target_id": mapped["target_id"],
                "gene_symbol": mapped["Gene_Symbol"],
                "uniprot": mapped["target_id"],
                "assay_id": "ExCAPE:" + mapped["DB"].astype(str) + ":" + mapped["Original_Assay_ID"].astype(str),
                "activity_nM": pd.NA,
                "activity_relation": "=",
                "activity_type": "pXC50",
                "activity_units": pd.NA,
                "pchembl_value": pd.to_numeric(mapped["pXC50"], errors="coerce"),
                "source": "ExCAPE-DB:" + mapped["DB"].fillna("unknown").astype(str),
                "source_specific_activity_label": mapped["Activity_Flag"].map(_activity_label),
                "activity_document_ids": mapped["Original_Entry_ID"],
                "pubchem_or_chembl_id": mapped["Original_Entry_ID"],
                "entrez_id": mapped["Entrez_ID"],
                "tax_id": mapped["Tax_ID"],
                "inchikey": mapped["ambit_inchikey"],
                "smiles": mapped["SMILES"],
            }
        )
        staged.to_csv(out_file, sep="\t", index=False, mode="w" if first else "a", header=first)
        first = False
    if first:
        pd.DataFrame(
            columns=[
                "drug_id",
                "drug_name",
                "target_id",
                "gene_symbol",
                "uniprot",
                "assay_id",
                "activity_nM",
                "activity_relation",
                "activity_type",
                "activity_units",
                "pchembl_value",
                "source",
                "source_specific_activity_label",
                "activity_document_ids",
                "pubchem_or_chembl_id",
                "entrez_id",
                "tax_id",
                "inchikey",
                "smiles",
            ]
        ).to_csv(out_file, sep="\t", index=False)
    manifest = {
        "source": "ExCAPE-DB PubChem/ChEMBL activity table",
        "input": str(excape_path),
        "output": str(out_file),
        "rows_scanned": int(rows_in),
        "rows_matching_human_atlas_targets": int(rows_human_target),
        "rows_mapped_to_atlas_drugs": int(rows_drug_mapped),
        "target_genes": int(len(target_genes)),
        "atlas_drug_keys": int(len(drug_map)),
        "warning": "ExCAPE is aggregated from PubChem/ChEMBL lineage; use as training-only or sensitivity evidence with source-holdout controls.",
    }
    out_file.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return pd.read_csv(out_file, sep="\t", low_memory=False)

from __future__ import annotations

import json
import re
import zipfile
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.external.compound_index import load_atlas_compound_index, map_by_inchikey


ENDPOINT_PRIORITY = ("Kd (nM)", "Ki (nM)", "IC50 (nM)", "EC50 (nM)")
BASE_COLUMNS = [
    "BindingDB Reactant_set_id",
    "Ligand SMILES",
    "Ligand InChI Key",
    "BindingDB MonomerID",
    "BindingDB Ligand Name",
    "Target Name",
    "Target Source Organism According to Curator or DataSource",
    "Curation/DataSource",
    "PMID",
    "PubChem AID",
    "Date of publication",
    *ENDPOINT_PRIORITY,
]


def _clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def _target_ids(pair_table_path: str | Path) -> set[str]:
    pairs = pd.read_csv(pair_table_path, low_memory=False)
    ids: set[str] = set()
    for col in ("target_id", "target_uniprot", "uniprot"):
        if col in pairs.columns:
            ids.update(pairs[col].dropna().astype(str).str.strip())
    return {value for value in ids if value and value.lower() != "nan"}


def _bindingdb_usecols(header: list[str]) -> list[str]:
    target_cols = [
        col
        for col in header
        if re.match(r"UniProt \((SwissProt|TrEMBL)\) Primary ID of Target Chain \d+", col)
    ]
    return [col for col in BASE_COLUMNS if col in header] + target_cols


def _zip_member(path: str | Path) -> str:
    with zipfile.ZipFile(path) as zf:
        members = [name for name in zf.namelist() if name.lower().endswith(".tsv")]
    if not members:
        raise ValueError(f"BindingDB ZIP has no TSV member: {path}")
    return members[0]


def _header(path: str | Path) -> list[str]:
    member = _zip_member(path)
    with zipfile.ZipFile(path) as zf, zf.open(member) as handle:
        line = handle.readline().decode("utf-8", errors="replace").rstrip("\n")
    return line.split("\t")


def _first_matching_target(chunk: pd.DataFrame, target_set: set[str]) -> pd.Series:
    target = pd.Series(pd.NA, index=chunk.index, dtype="object")
    target_cols = [col for col in chunk.columns if col.startswith("UniProt (") and "Primary ID" in col]
    for col in target_cols:
        values = chunk[col].map(_clean)
        match = target.isna() & values.isin(target_set)
        target = target.mask(match, values)
    return target


def _parse_endpoint(value: Any) -> tuple[str, float | None]:
    text = _clean(value).replace(",", "")
    if not text:
        return "", None
    relation = "="
    for prefix in ("<=", ">=", "<", ">"):
        if text.startswith(prefix):
            relation = prefix
            text = text[len(prefix) :].strip()
            break
    match = re.search(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", text)
    if not match:
        return relation, None
    return relation, float(match.group(0))


def _best_endpoint(chunk: pd.DataFrame) -> pd.DataFrame:
    endpoint = pd.Series(pd.NA, index=chunk.index, dtype="object")
    relation = pd.Series("", index=chunk.index, dtype="object")
    value = pd.Series(pd.NA, index=chunk.index, dtype="Float64")
    for col in ENDPOINT_PRIORITY:
        if col not in chunk.columns:
            continue
        missing = value.isna()
        if not missing.any():
            break
        parsed = chunk.loc[missing, col].map(_parse_endpoint)
        parsed_relation = parsed.map(lambda item: item[0])
        parsed_value = parsed.map(lambda item: item[1])
        has_value = parsed_value.notna()
        idx = parsed_value[has_value].index
        endpoint.loc[idx] = col.replace(" (nM)", "")
        relation.loc[idx] = parsed_relation.loc[idx]
        value.loc[idx] = parsed_value.loc[idx].astype(float)
    return pd.DataFrame({"activity_type": endpoint, "activity_relation": relation, "activity_nM": value})


def _activity_label(relation: pd.Series, value_nm: pd.Series) -> pd.Series:
    relation_text = relation.fillna("").astype(str).str.strip()
    values = pd.to_numeric(value_nm, errors="coerce")
    active = relation_text.isin({"", "=", "<", "<="}) & values.le(1000.0)
    exact = relation_text.isin({"", "="})
    lower_bound = relation_text.isin({">", ">="})
    inactive = (exact | lower_bound) & values.ge(10000.0)
    label = pd.Series("inconclusive", index=value_nm.index, dtype="object")
    label = label.mask(active, "active")
    label = label.mask(inactive, "inactive")
    return label


def stage_bindingdb_for_atlas(
    bindingdb_zip: str | Path,
    pair_table_path: str | Path,
    mapping_path: str | Path,
    out_path: str | Path,
    *,
    chunksize: int = 200_000,
) -> pd.DataFrame:
    """Stream-filter BindingDB to Atlas FDA drugs and selected target UniProt IDs."""

    target_set = _target_ids(pair_table_path)
    atlas_index = load_atlas_compound_index(mapping_path)
    if not target_set:
        raise ValueError("BindingDB staging requires target IDs from the pair table")
    if atlas_index.empty:
        raise ValueError("BindingDB staging requires a nonempty Atlas compound mapping")

    header = _header(bindingdb_zip)
    usecols = _bindingdb_usecols(header)
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    first = True
    rows_scanned = 0
    rows_matching_targets = 0
    rows_mapped_drugs = 0
    rows_written = 0

    reader = pd.read_csv(
        bindingdb_zip,
        sep="\t",
        compression="zip",
        usecols=usecols,
        chunksize=chunksize,
        low_memory=False,
    )
    for chunk in reader:
        rows_scanned += len(chunk)
        chunk["_target_id"] = _first_matching_target(chunk, target_set)
        matched = chunk[chunk["_target_id"].notna()].copy()
        if matched.empty:
            continue
        rows_matching_targets += len(matched)
        mapped = map_by_inchikey(matched.rename(columns={"Ligand InChI Key": "inchikey"}), atlas_index)
        mapped = mapped[mapped["drug_id"].notna()].copy()
        if mapped.empty:
            continue
        rows_mapped_drugs += len(mapped)
        endpoints = _best_endpoint(mapped)
        staged = pd.DataFrame(
            {
                "drug_id": mapped["drug_id"],
                "drug_name": mapped["drug_name"],
                "target_id": mapped["_target_id"],
                "uniprot": mapped["_target_id"],
                "assay_id": "BindingDB:" + mapped["BindingDB Reactant_set_id"].astype(str),
                "activity_nM": endpoints["activity_nM"],
                "activity_relation": endpoints["activity_relation"],
                "activity_type": endpoints["activity_type"],
                "activity_units": "nM",
                "source": "BindingDB:" + mapped["Curation/DataSource"].fillna("unknown").astype(str),
                "source_specific_activity_label": _activity_label(
                    endpoints["activity_relation"], endpoints["activity_nM"]
                ),
                "bindingdb_reactant_set_id": mapped["BindingDB Reactant_set_id"],
                "bindingdb_monomerid": mapped["BindingDB MonomerID"],
                "bindingdb_ligand_name": mapped["BindingDB Ligand Name"],
                "target_name": mapped["Target Name"],
                "target_organism": mapped["Target Source Organism According to Curator or DataSource"],
                "document_ids": mapped["PMID"],
                "pubchem_aid": mapped["PubChem AID"],
                "publication_year": pd.to_datetime(
                    mapped["Date of publication"], errors="coerce"
                ).dt.year,
                "inchikey": mapped["inchikey"],
                "smiles": mapped["Ligand SMILES"],
            }
        )
        staged = staged[staged["activity_nM"].notna()].drop_duplicates()
        if staged.empty:
            continue
        rows_written += len(staged)
        staged.to_csv(out, sep="\t", index=False, mode="w" if first else "a", header=first)
        first = False

    if first:
        pd.DataFrame(
            columns=[
                "drug_id",
                "drug_name",
                "target_id",
                "uniprot",
                "assay_id",
                "activity_nM",
                "activity_relation",
                "activity_type",
                "activity_units",
                "source",
                "source_specific_activity_label",
                "bindingdb_reactant_set_id",
                "bindingdb_monomerid",
                "bindingdb_ligand_name",
                "target_name",
                "target_organism",
                "document_ids",
                "pubchem_aid",
                "publication_year",
                "inchikey",
                "smiles",
            ]
        ).to_csv(out, sep="\t", index=False)

    manifest = {
        "source": "BindingDB",
        "input": str(bindingdb_zip),
        "output": str(out),
        "rows_scanned": int(rows_scanned),
        "rows_matching_atlas_targets": int(rows_matching_targets),
        "rows_mapped_to_atlas_drugs": int(rows_mapped_drugs),
        "rows_written": int(rows_written),
        "target_ids": int(len(target_set)),
        "atlas_drug_rows": int(len(atlas_index)),
        "label_policy": "active if best Kd/Ki/IC50/EC50 <= 1000 nM; inactive if exact or censored >/>= measured values are >= 10000 nM; otherwise gray/unknown.",
    }
    out.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return pd.read_csv(out, sep="\t", low_memory=False)

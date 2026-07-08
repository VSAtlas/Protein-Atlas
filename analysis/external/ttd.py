from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from analysis.external.compound_index import (
    load_atlas_compound_index,
    map_by_inchikey,
    map_by_pubchem_cid,
)

Chem: Any
try:  # pragma: no cover - exercised only when RDKit is installed.
    from rdkit import Chem
except Exception:  # pragma: no cover
    Chem = None


ACTIVITY_RE = re.compile(
    r"^\s*(?P<endpoint>[A-Za-z0-9_ /.-]+?)\s*"
    r"(?P<relation><=|>=|=|<|>|~)\s*"
    r"(?P<value>[0-9.eE+-]+)\s*"
    r"(?P<unit>[A-Za-zµuUMnp]+)\s*$"
)


def _clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def _inchi_to_inchikey(value: Any) -> str:
    text = _clean(value)
    if not text or Chem is None:
        return ""
    mol = Chem.MolFromInchi(text if text.startswith("InChI=") else f"InChI={text}")
    if mol is None:
        return ""
    return str(Chem.MolToInchiKey(mol) or "")


def _smiles_to_inchikey(value: Any) -> str:
    text = _clean(value)
    if not text or Chem is None:
        return ""
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return ""
    return str(Chem.MolToInchiKey(mol) or "")


def _parse_ttd_key_value_file(path: str | Path, id_fields: set[str]) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    current: dict[str, Any] = {}
    current_id = ""
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith("TTD -") or text.startswith("Title -") or text.startswith("Version "):
                continue
            fields = [part.strip() for part in line.rstrip("\n").split("\t")]
            fields = [part for part in fields if part != ""]
            if len(fields) < 2:
                continue
            if fields[0] in id_fields:
                if fields[1].lower().startswith(("ttd ", "drug ", "target ")):
                    continue
                if current:
                    records.append(current)
                current_id = fields[1]
                current = {fields[0]: current_id}
                continue
            if fields[0].startswith("T") and len(fields) >= 3:
                row_id, key, value = fields[0], fields[1], " ".join(fields[2:])
                if row_id != current_id:
                    if current:
                        records.append(current)
                    current_id = row_id
                    current = {"TTD_ID": row_id}
                current[key] = value
                continue
            if current and len(fields) >= 2:
                current[fields[0]] = " ".join(fields[1:])
    if current:
        records.append(current)
    return pd.DataFrame(records)


def _first_token(value: Any, separators: tuple[str, ...] = (";", "|", ",")) -> str:
    text = _clean(value)
    for sep in separators:
        if sep in text:
            text = text.split(sep, 1)[0]
    return text.strip()


def load_ttd_targets(target_path: str | Path) -> pd.DataFrame:
    targets = _parse_ttd_key_value_file(target_path, {"TARGETID"})
    if targets.empty:
        return pd.DataFrame(columns=["ttd_target_id", "gene_symbol", "target_name", "swissprot"])
    id_col = "TARGETID" if "TARGETID" in targets.columns else "TTD_ID"
    out = pd.DataFrame(
        {
            "ttd_target_id": targets[id_col].map(_clean),
            "gene_symbol": targets.get("GENENAME", pd.Series(pd.NA, index=targets.index)).map(_first_token),
            "target_name": targets.get("TARGNAME", pd.Series(pd.NA, index=targets.index)).map(_clean),
            "swissprot": targets.get("UNIPROID", pd.Series(pd.NA, index=targets.index)).map(_first_token),
        }
    )
    return out.drop_duplicates("ttd_target_id")


def load_ttd_approved_structures(path: str | Path) -> pd.DataFrame:
    structs = _parse_ttd_key_value_file(path, {"DRUG__ID"})
    if structs.empty:
        return pd.DataFrame(columns=["ttd_drug_id", "ttd_drug_name", "smiles", "inchikey"])
    out = pd.DataFrame(
        {
            "ttd_drug_id": structs["DRUG__ID"].map(_clean),
            "ttd_drug_name": structs.get("DRUGNAME", pd.Series(pd.NA, index=structs.index)).map(_clean),
            "smiles": structs.get("DRUGSMIL", pd.Series(pd.NA, index=structs.index)).map(_clean),
            "inchi": structs.get("DRUGINCH", pd.Series(pd.NA, index=structs.index)).map(_clean),
        }
    )
    out["inchikey"] = out["inchi"].map(_inchi_to_inchikey)
    missing = out["inchikey"].eq("")
    out.loc[missing, "inchikey"] = out.loc[missing, "smiles"].map(_smiles_to_inchikey)
    return out.drop_duplicates("ttd_drug_id")


def load_ttd_crossmatch(path: str | Path) -> pd.DataFrame:
    rows = _parse_ttd_key_value_file(path, set())
    if rows.empty:
        return pd.DataFrame(columns=["ttd_drug_id", "pubchem_cid", "ttd_drug_name"])
    id_col = "TTD_ID" if "TTD_ID" in rows.columns else rows.columns[0]
    out = pd.DataFrame(
        {
            "ttd_drug_id": rows[id_col].map(_clean),
            "pubchem_cid": pd.to_numeric(
                rows.get("PUBCHCID", pd.Series(pd.NA, index=rows.index)).map(_first_token),
                errors="coerce",
            ).astype("Int64"),
            "ttd_drug_name": rows.get("DRUGNAME", pd.Series(pd.NA, index=rows.index)).map(_clean),
        }
    )
    return out.drop_duplicates("ttd_drug_id")


def _parse_activity(text: Any) -> dict[str, Any]:
    value = _clean(text)
    match = ACTIVITY_RE.match(value)
    if not match:
        return {
            "activity_type": pd.NA,
            "activity_relation": pd.NA,
            "activity_value": pd.NA,
            "activity_units": pd.NA,
        }
    payload = match.groupdict()
    endpoint = payload["endpoint"].strip().replace(" ", "_")
    return {
        "activity_type": endpoint,
        "activity_relation": payload["relation"],
        "activity_value": float(payload["value"]),
        "activity_units": payload["unit"],
    }


def _to_nm(value: pd.Series, unit: pd.Series) -> pd.Series:
    out = pd.to_numeric(value, errors="coerce")
    units = unit.astype("object").where(unit.notna(), "").astype(str).str.lower()
    out = out.where(~units.isin({"um", "µm", "uM".lower(), "micromolar"}), out * 1000.0)
    out = out.where(~units.isin({"mm", "millimolar"}), out * 1_000_000.0)
    out = out.where(~units.isin({"pm", "picomolar"}), out / 1000.0)
    return out


def _pchembl_from_nm(value_nm: pd.Series) -> pd.Series:
    value = pd.to_numeric(value_nm, errors="coerce")
    return 9.0 - value.apply(lambda item: pd.NA if pd.isna(item) or item <= 0 else np.log10(item))


def _target_map_from_feature_table(feature_table_path: str | Path | None) -> pd.DataFrame:
    if feature_table_path is None or not Path(feature_table_path).exists():
        return pd.DataFrame(columns=["gene_symbol", "target_id"])
    df = pd.read_csv(feature_table_path, low_memory=False)
    gene_col = next((col for col in ("target_gene", "gene_symbol", "target_symbol") if col in df.columns), None)
    target_col = next((col for col in ("target_id", "target_uniprot", "uniprot") if col in df.columns), None)
    if gene_col is None or target_col is None:
        return pd.DataFrame(columns=["gene_symbol", "target_id"])
    out = df[[gene_col, target_col]].dropna().drop_duplicates()
    out.columns = ["gene_symbol", "target_id"]
    out["gene_symbol"] = out["gene_symbol"].astype(str).str.upper().str.strip()
    out["target_id"] = out["target_id"].astype(str).str.strip()
    return out.drop_duplicates("gene_symbol")


def stage_ttd_activity(
    activity_path: str | Path,
    target_path: str | Path,
    crossmatch_path: str | Path,
    structures_path: str | Path,
    mapping_path: str | Path,
    out_path: str | Path,
    *,
    feature_table_path: str | Path | None = None,
) -> pd.DataFrame:
    """Stage TTD measured target-compound activity as Atlas four-state input."""

    activity = pd.read_csv(activity_path, sep="\t", low_memory=False)
    parsed = pd.DataFrame([_parse_activity(item) for item in activity["Activity"]])
    out = pd.concat([activity, parsed], axis=1)
    out["activity_nM"] = _to_nm(out["activity_value"], out["activity_units"])
    out["pchembl_value"] = _pchembl_from_nm(out["activity_nM"])
    out = out.rename(
        columns={
            "TTD Target ID": "ttd_target_id",
            "TTD Drug/Compound ID": "ttd_drug_id",
            "Pubchem CID": "pubchem_cid",
        }
    )
    targets = load_ttd_targets(target_path)
    out = out.merge(targets, on="ttd_target_id", how="left")
    target_map = _target_map_from_feature_table(feature_table_path)
    out["gene_key"] = out["gene_symbol"].astype("object").where(out["gene_symbol"].notna(), "").astype(str).str.upper()
    if not target_map.empty:
        out = out.merge(target_map, left_on="gene_key", right_on="gene_symbol", how="left", suffixes=("", "_atlas"))
    if "target_id" not in out.columns:
        out["target_id"] = pd.NA
    out["target_id"] = out["target_id"].fillna(out["gene_symbol"]).fillna(out["swissprot"])

    cross = load_ttd_crossmatch(crossmatch_path)
    structs = load_ttd_approved_structures(structures_path)
    out = out.merge(cross, on="ttd_drug_id", how="left", suffixes=("", "_cross"))
    out["pubchem_cid"] = pd.to_numeric(out["pubchem_cid"], errors="coerce").astype("Int64")
    out["pubchem_cid"] = out["pubchem_cid"].fillna(out["pubchem_cid_cross"])
    out = out.merge(structs, on="ttd_drug_id", how="left", suffixes=("", "_structure"))

    atlas = load_atlas_compound_index(mapping_path)
    out = map_by_pubchem_cid(out, atlas, source_cid_col="pubchem_cid")
    missing = out["drug_id"].isna()
    if missing.any():
        mapped_missing = map_by_inchikey(out.loc[missing].copy(), atlas, source_inchikey_col="inchikey")
        for column in ("drug_id", "drug_name", "atlas_smiles", "atlas_pubchem_cid"):
            if column in mapped_missing.columns:
                out.loc[missing, column] = mapped_missing[column].to_numpy()
    out["drug_id"] = out["drug_id"].fillna("ttd:" + out["ttd_drug_id"].astype(str))
    cross_name = out["ttd_drug_name_cross"] if "ttd_drug_name_cross" in out.columns else pd.Series(pd.NA, index=out.index)
    structure_name = out["ttd_drug_name"] if "ttd_drug_name" in out.columns else pd.Series(pd.NA, index=out.index)
    out["drug_name"] = out["drug_name"].fillna(structure_name).fillna(cross_name)
    out["assay_id"] = "TTD:" + out["ttd_target_id"].astype(str) + ":" + out["ttd_drug_id"].astype(str)
    out["source"] = "TTD"
    staged = pd.DataFrame(
        {
            "drug_id": out["drug_id"],
            "drug_name": out["drug_name"],
            "target_id": out["target_id"],
            "gene_symbol": out["gene_symbol"],
            "uniprot": out["target_id"],
            "assay_id": out["assay_id"],
            "activity_nM": out["activity_nM"],
            "activity_relation": out["activity_relation"],
            "activity_type": out["activity_type"],
            "activity_units": "nM",
            "pchembl_value": out["pchembl_value"],
            "source": out["source"],
            "activity_document_ids": out["ttd_target_id"].astype(str) + ";" + out["ttd_drug_id"].astype(str),
            "source_specific_activity_label": pd.NA,
            "ttd_target_id": out["ttd_target_id"],
            "ttd_drug_id": out["ttd_drug_id"],
            "pubchem_cid": out["pubchem_cid"],
            "smiles": out.get("smiles"),
            "inchikey": out.get("inchikey"),
        }
    )
    out_file = Path(out_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    staged.to_csv(out_file, sep="\t", index=False)
    manifest = {
        "source": "TTD target-compound activity",
        "rows": int(len(staged)),
        "atlas_drug_mapped_rows": int((~staged["drug_id"].astype(str).str.startswith("ttd:")).sum()) if not staged.empty else 0,
        "target_mapped_rows": int(staged["target_id"].notna().sum()),
        "output": str(out_file),
    }
    out_file.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return staged

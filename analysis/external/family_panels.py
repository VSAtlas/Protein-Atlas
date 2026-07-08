from __future__ import annotations

import io
import json
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.external.compound_index import atlas_drug_id, load_atlas_compound_index


CHEMBL_MOLECULE_API = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"


def _json(url: str) -> dict[str, Any]:
    request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "AtlasAnalysis/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode("utf-8"))


def resolve_atlas_chembl_molecules(
    mapping_path: str | Path,
    out_path: str | Path,
    *,
    sleep_sec: float = 0.05,
    overwrite: bool = False,
) -> pd.DataFrame:
    """Resolve Atlas/FDA InChIKeys to ChEMBL molecule IDs for family-panel joins."""

    out = Path(out_path)
    if out.exists() and not overwrite:
        return pd.read_csv(out, low_memory=False)
    mapping = pd.read_csv(mapping_path, low_memory=False)
    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    if out.exists() and overwrite:
        try:
            existing = pd.read_csv(out, low_memory=False)
            rows = existing.to_dict("records")
            seen = {
                (str(row.get("drug_id") or ""), str(row.get("inchikey") or ""))
                for row in rows
                if str(row.get("drug_id") or "") and str(row.get("inchikey") or "")
            }
        except Exception:
            rows = []
            seen = set()
    pending: list[dict[str, Any]] = []
    for row in mapping.to_dict("records"):
        inchikey = str(row.get("inchikey") or row.get("remark_inchikey") or "").strip()
        if not inchikey:
            continue
        drug_id = atlas_drug_id(row)
        key = (drug_id, inchikey)
        if key in seen:
            continue
        pending.append(dict(row, _atlas_drug_id=drug_id, _atlas_inchikey=inchikey))

    for start in range(0, len(pending), 100):
        chunk = pending[start : start + 100]
        key_rows = {str(row["_atlas_inchikey"]): row for row in chunk}
        inchikeys = sorted(key_rows)
        query = urllib.parse.urlencode(
            {"molecule_structures__standard_inchi_key__in": ",".join(inchikeys), "limit": 1000}
        )
        returned: set[str] = set()
        try:
            payload = _json(f"{CHEMBL_MOLECULE_API}?{query}")
        except Exception as exc:
            for row in chunk:
                inchikey = str(row["_atlas_inchikey"])
                drug_id = str(row["_atlas_drug_id"])
                rows.append(
                    {
                        "drug_id": drug_id,
                        "inchikey": inchikey,
                        "molecule_chembl_id": pd.NA,
                        "chembl_pref_name": pd.NA,
                        "mapping_status": "chembl_query_failed",
                        "mapping_error": str(exc),
                    }
                )
                seen.add((drug_id, inchikey))
            continue
        for molecule in payload.get("molecules", []):
            structures = molecule.get("molecule_structures") or {}
            inchikey = str(structures.get("standard_inchi_key") or "").strip()
            source_row = key_rows.get(inchikey)
            if source_row is None:
                continue
            returned.add(inchikey)
            rows.append(
                {
                    "drug_id": source_row["_atlas_drug_id"],
                    "inchikey": inchikey,
                    "drug_name": source_row.get("display_name")
                    or source_row.get("generic_name")
                    or source_row.get("pubchem_name"),
                    "molecule_chembl_id": molecule.get("molecule_chembl_id"),
                    "chembl_pref_name": molecule.get("pref_name"),
                    "mapping_status": "exact_inchikey_match",
                    "mapping_error": "",
                }
            )
            seen.add((str(source_row["_atlas_drug_id"]), inchikey))
        for inchikey in sorted(set(inchikeys) - returned):
            source_row = key_rows[inchikey]
            rows.append(
                {
                    "drug_id": source_row["_atlas_drug_id"],
                    "inchikey": inchikey,
                    "drug_name": source_row.get("display_name")
                    or source_row.get("generic_name")
                    or source_row.get("pubchem_name"),
                    "molecule_chembl_id": pd.NA,
                    "chembl_pref_name": pd.NA,
                    "mapping_status": "no_chembl_molecule_for_inchikey",
                    "mapping_error": "",
                }
            )
            seen.add((str(source_row["_atlas_drug_id"]), inchikey))
        out.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).drop_duplicates().to_csv(out, index=False)
        time.sleep(max(0.0, sleep_sec))
    result = pd.DataFrame(rows).drop_duplicates()
    out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out, index=False)
    return result

def _read_kiba_workbook(kiba_zip_path: str | Path) -> bytes:
    with zipfile.ZipFile(kiba_zip_path) as archive:
        return archive.read("KiBA.xlsx")


def _kiba_matrix_to_rows(workbook: bytes, chembl_map: pd.DataFrame) -> pd.DataFrame:
    matrix = pd.read_excel(io.BytesIO(workbook), sheet_name="KIBA")
    matrix = matrix.rename(columns={matrix.columns[0]: "molecule_chembl_id"})
    long = matrix.melt(id_vars=["molecule_chembl_id"], var_name="target_id", value_name="kiba_score")
    long = long.dropna(subset=["kiba_score"])
    long["molecule_chembl_id"] = long["molecule_chembl_id"].astype(str).str.strip()
    mapped = long.merge(
        chembl_map.dropna(subset=["molecule_chembl_id"]).drop_duplicates(["molecule_chembl_id", "drug_id"]),
        on="molecule_chembl_id",
        how="inner",
    )
    if mapped.empty:
        return pd.DataFrame(columns=_BIOACTIVITY_COLUMNS)
    score = pd.to_numeric(mapped["kiba_score"], errors="coerce")
    label = pd.Series(pd.NA, index=mapped.index, dtype="Int64")
    label = label.mask(score.ge(12.1), 1)
    label = label.mask(score.le(4.0), 0)
    return pd.DataFrame(
        {
            "drug_id": mapped["drug_id"],
            "target_id": mapped["target_id"],
            "activity_nM": pd.NA,
            "pchembl_value": pd.NA,
            "activity_relation": "",
            "activity_type": "KIBA integrated kinase bioactivity score",
            "activity_units": "KIBA_score",
            "assay_id": "KIBA",
            "source": "KIBA",
            "source_specific_activity_label": label,
            "activity_document_ids": "Tang2014_JCIM;Zenodo5105698",
            "activity_publication_year": 2014,
            "database_release_year": 2021,
            "kiba_score": score,
            "label_status": pd.Series("kiba_gray_zone", index=mapped.index)
            .mask(label.eq(1).fillna(False), "kiba_active_score_ge_12_1")
            .mask(label.eq(0).fillna(False), "kiba_inactive_score_le_4_0"),
            "molecule_chembl_id": mapped["molecule_chembl_id"],
        },
        columns=_BIOACTIVITY_COLUMNS + ["kiba_score", "label_status", "molecule_chembl_id"],
    )


def _direct_sheet_to_rows(workbook: bytes, sheet_name: str, chembl_map: pd.DataFrame) -> pd.DataFrame:
    raw = pd.read_excel(io.BytesIO(workbook), sheet_name=sheet_name)
    if "compound_id" not in raw.columns or "protein_id" not in raw.columns:
        return pd.DataFrame(columns=_BIOACTIVITY_COLUMNS)
    value_col = next((col for col in ("ki (nM)", "kd (nM)", "Ki (nM)", "Kd (nM)") if col in raw.columns), None)
    if value_col is None:
        return pd.DataFrame(columns=_BIOACTIVITY_COLUMNS)
    raw["molecule_chembl_id"] = raw["compound_id"].astype(str).str.strip()
    mapped = raw.merge(
        chembl_map.dropna(subset=["molecule_chembl_id"]).drop_duplicates(["molecule_chembl_id", "drug_id"]),
        on="molecule_chembl_id",
        how="inner",
    )
    if mapped.empty:
        return pd.DataFrame(columns=_BIOACTIVITY_COLUMNS)
    values = pd.to_numeric(mapped[value_col], errors="coerce")
    values = values.where(values.gt(0))
    return pd.DataFrame(
        {
            "drug_id": mapped["drug_id"],
            "target_id": mapped["protein_id"].astype(str).str.strip(),
            "activity_nM": values,
            "pchembl_value": pd.NA,
            "activity_relation": "=",
            "activity_type": value_col.split()[0].upper(),
            "activity_units": "nM",
            "assay_id": sheet_name,
            "source": f"KIBA/{sheet_name}",
            "source_specific_activity_label": pd.NA,
            "activity_document_ids": "Tang2014_JCIM;Zenodo5105698",
            "activity_publication_year": 2014,
            "database_release_year": 2021,
        },
        columns=_BIOACTIVITY_COLUMNS,
    ).dropna(subset=["activity_nM"])


_BIOACTIVITY_COLUMNS = [
    "drug_id",
    "target_id",
    "activity_nM",
    "pchembl_value",
    "activity_relation",
    "activity_type",
    "activity_units",
    "assay_id",
    "source",
    "source_specific_activity_label",
    "activity_document_ids",
    "activity_publication_year",
    "database_release_year",
]


def stage_kiba_family_panel(
    kiba_zip_path: str | Path,
    mapping_path: str | Path,
    out_path: str | Path,
    *,
    chembl_mapping_path: str | Path | None = None,
    resolve_chembl: bool = False,
    sleep_sec: float = 0.05,
) -> pd.DataFrame:
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    chembl_map_path = Path(chembl_mapping_path) if chembl_mapping_path is not None else out.with_suffix(".chembl_atlas_mapping.csv")
    if resolve_chembl or not chembl_map_path.exists():
        chembl_map = resolve_atlas_chembl_molecules(mapping_path, chembl_map_path, sleep_sec=sleep_sec)
    else:
        chembl_map = pd.read_csv(chembl_map_path, low_memory=False)
    workbook = _read_kiba_workbook(kiba_zip_path)
    parts = [_kiba_matrix_to_rows(workbook, chembl_map)]
    for sheet in ("Davis_Metz", "Davis_Anastassiadis", "Metz_Anastassiadis"):
        parts.append(_direct_sheet_to_rows(workbook, sheet, chembl_map))
    result = pd.concat(parts, ignore_index=True).drop_duplicates()
    result.to_csv(out, sep="\t", index=False)
    manifest = {
        "source": "KIBA kinase family panel",
        "input": str(kiba_zip_path),
        "output": str(out),
        "chembl_mapping": str(chembl_map_path),
        "rows": int(len(result)),
        "mapped_atlas_drugs": int(result["drug_id"].nunique()) if not result.empty else 0,
        "mapped_targets": int(result["target_id"].nunique()) if not result.empty else 0,
        "label_policy": "KIBA score >= 12.1 -> active; KIBA score <= 4.0 -> inactive; direct Ki/Kd sheets use nM thresholds in four-state builder.",
    }
    out.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return result


CARDIAC_TARGETS = {
    "herg": {"target_id": "Q12809", "gene_symbol": "KCNH2", "assay_id": "hERG"},
    "cav1.2": {"target_id": "Q13936", "gene_symbol": "CACNA1C", "assay_id": "Cav1.2"},
    "nav1.5": {"target_id": "Q14524", "gene_symbol": "SCN5A", "assay_id": "Nav1.5"},
}


def stage_cardiac_ion_channel_panel(
    extracted_root: str | Path,
    mapping_path: str | Path,
    out_path: str | Path,
) -> pd.DataFrame:
    atlas_index = load_atlas_compound_index(mapping_path)
    root = Path(extracted_root)
    rows: list[pd.DataFrame] = []
    for csv_path in sorted(root.rglob("*.csv")):
        target_key = next((key for key in CARDIAC_TARGETS if key in str(csv_path).lower()), "")
        if not target_key:
            continue
        target = CARDIAC_TARGETS[target_key]
        raw = pd.read_csv(csv_path, low_memory=False)
        inchikey_col = next((col for col in raw.columns if col.lower().replace("l", "i").startswith("inchi key")), "InChl Key")
        raw = raw.rename(columns={inchikey_col: "inchikey", "pIC50": "pIC50"})
        from analysis.external.compound_index import map_by_inchikey

        mapped = map_by_inchikey(raw, atlas_index, source_inchikey_col="inchikey")
        mapped = mapped[mapped["drug_id"].notna()].copy()
        if mapped.empty:
            continue
        p_ic50 = pd.to_numeric(mapped["pIC50"], errors="coerce")
        activity_nm = 10 ** (9 - p_ic50)
        rows.append(
            pd.DataFrame(
                {
                    "drug_id": mapped["drug_id"],
                    "target_id": target["target_id"],
                    "activity_nM": activity_nm,
                    "pchembl_value": p_ic50,
                    "activity_relation": "=",
                    "activity_type": "IC50",
                    "activity_units": "nM",
                    "assay_id": target["assay_id"],
                    "source": "Cardiac ion-channel panel",
                    "source_specific_activity_label": pd.NA,
                    "activity_document_ids": "Arab2023_Zenodo8245086",
                    "activity_publication_year": 2023,
                    "database_release_year": 2023,
                    "cardiac_gene_symbol": target["gene_symbol"],
                    "cardiac_source_file": str(csv_path),
                    "cardiac_original_source": mapped.get("Source", pd.Series(pd.NA, index=mapped.index)),
                },
                columns=_BIOACTIVITY_COLUMNS + ["cardiac_gene_symbol", "cardiac_source_file", "cardiac_original_source"],
            )
        )
    result = pd.concat(rows, ignore_index=True).drop_duplicates() if rows else pd.DataFrame(columns=_BIOACTIVITY_COLUMNS)
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out, sep="\t", index=False)
    manifest = {
        "source": "cardiac ion-channel panel",
        "input_root": str(root),
        "output": str(out),
        "rows": int(len(result)),
        "mapped_atlas_drugs": int(result["drug_id"].nunique()) if not result.empty else 0,
        "mapped_targets": int(result["target_id"].nunique()) if not result.empty else 0,
        "target_mapping": CARDIAC_TARGETS,
        "label_policy": "pIC50 converted to IC50 nM; four-state builder applies <=1 uM active and >=10 uM inactive thresholds.",
    }
    out.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return result

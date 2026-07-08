from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.external.papyrus_chembl import build_bioactivity_benchmark
from analysis.ml.build_ml_dataset import build_ml_dataset
from analysis.ml.splits import make_split, split_overlap_summary


PAPYRUS_PLUS_URL = "https://zenodo.org/api/records/13987985/files/05.7++_combined_set_without_stereochemistry.tsv.xz/content"
CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
TOXCAST_SUMMARY_URL = "https://clowder.edap-cluster.com/files/68af6b70e4b02565fc7c3a98/blob"
TOXCAST_SUMMARY_SHA256 = "c61d46cf685d5d84e72d0d4f6f46746ab1714d6f924e43480d390925ae827c59"


def _text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _http_json(url: str, *, sleep_sec: float = 0.08) -> dict[str, Any]:
    time.sleep(max(0.0, sleep_sec))
    request = urllib.request.Request(url, headers={"User-Agent": "AtlasAnalysis/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.load(response)


def _download(url: str, out_path: Path, *, overwrite: bool = False) -> Path:
    if out_path.exists() and not overwrite:
        return out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "AtlasAnalysis/1.0"})
    with urllib.request.urlopen(request, timeout=300) as response:
        with out_path.open("wb") as handle:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
    return out_path


def resolve_fda_mapping_path(repo_root: str | Path) -> Path:
    root = Path(repo_root)
    candidates = [
        root / "chemdb" / "data" / "fda_mapping_from_pdbqt.csv",
        root / "fda_mapping_from_pdbqt.csv",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"No FDA PDBQT mapping found under {root}")


def parse_pdb_header_mapping(input_pdbs_dir: str | Path) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for path in sorted(Path(input_pdbs_dir).glob("*.pdb")):
        pdb_id = path.stem.upper()
        gene = ""
        uniprot = ""
        target_name = ""
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if line.startswith("ATOM"):
                    break
                if line.startswith("COMPND") and "MOLECULE:" in line and not target_name:
                    target_name = line.split("MOLECULE:", 1)[1].split(";", 1)[0].strip()
                elif line.startswith("SOURCE") and "GENE:" in line and not gene:
                    gene = line.split("GENE:", 1)[1].split(";", 1)[0].strip().split(",")[0]
                elif line.startswith("DBREF") and " UNP " in line and not uniprot:
                    parts = line.split()
                    if "UNP" in parts:
                        idx = parts.index("UNP")
                        if idx + 1 < len(parts):
                            uniprot = parts[idx + 1].strip()
        rows.append(
            {
                "pdb_id": pdb_id,
                "target_gene": gene.upper(),
                "target_uniprot": uniprot.upper(),
                "target_name": target_name,
                "target_id": uniprot.upper() or gene.upper() or pdb_id,
            }
        )
    return pd.DataFrame(rows)


def build_pilot_pair_table(
    master_rows_path: str | Path,
    fda_mapping_path: str | Path,
    input_pdbs_dir: str | Path,
    out_path: str | Path,
    pk_table_path: str | Path | None = None,
) -> pd.DataFrame:
    master = pd.read_csv(master_rows_path)
    master = master[master.get("is_decoy", "").fillna("").astype(str).str.lower().isin({"", "0", "0.0", "false"})].copy()
    mapping = pd.read_csv(fda_mapping_path)
    if "ligand_base" not in mapping.columns:
        mapping["ligand_base"] = mapping.apply(
            lambda row: f"{_text(row.get('scheme'))}_{int(float(row.get('file_num'))):07d}"
            if _text(row.get("scheme")) and _text(row.get("file_num"))
            else "",
            axis=1,
        )
    pdb_map = parse_pdb_header_mapping(input_pdbs_dir)
    out = master.merge(
        mapping[
            [
                col
                for col in (
                    "ligand_base",
                    "display_name",
                    "generic_name",
                    "pubchem_name",
                    "inchikey",
                    "smiles",
                    "drugcentral_id",
                )
                if col in mapping.columns
            ]
        ].drop_duplicates("ligand_base"),
        on="ligand_base",
        how="left",
    )
    out = out.merge(pdb_map, on="pdb_id", how="left", suffixes=("", "_pdb"))
    out["drug_id"] = out["ligand_base"]
    out["target_id"] = out["target_uniprot"].fillna("").where(out["target_uniprot"].fillna("").astype(str).str.len() > 0, out["target_gene"])
    out["ligand_chemotype"] = out["inchikey"].fillna("").astype(str).str.split("-").str[0]
    out["ligand_chemotype"] = out["ligand_chemotype"].where(
        out["ligand_chemotype"].astype(str).str.len() > 0,
        out["ligand_base"],
    )
    out["atlas_score"] = pd.to_numeric(out.get("z_selected"), errors="coerce")
    out["free_cmax_um"] = pd.to_numeric(out.get("free_cmax_um", out.get("free_cmax_uM")), errors="coerce")
    if pk_table_path is not None and Path(pk_table_path).exists():
        pk = pd.read_csv(pk_table_path)
        pk = add_standard_pk_columns(pk)
        pk_cols = [
            col
            for col in (
                "drug_id",
                "free_cmax_um",
                "cmax_um",
                "fraction_unbound_plasma",
                "pk_source",
                "pk_missing_reason",
            )
            if col in pk.columns
        ]
        pk = pk[pk_cols].drop_duplicates("drug_id")
        out = out.merge(pk, on="drug_id", how="left", suffixes=("", "_pk"))
        for col in ("free_cmax_um", "cmax_um", "fraction_unbound_plasma"):
            pk_col = f"{col}_pk"
            if pk_col in out.columns:
                out[col] = pd.to_numeric(out.get(col), errors="coerce").fillna(
                    pd.to_numeric(out[pk_col], errors="coerce")
                )
                out = out.drop(columns=[pk_col])
    cols = [
        col
        for col in (
            "drug_id",
            "target_id",
            "pdb_id",
            "variant",
            "ph_label",
            "ligand_base",
            "display_name",
            "generic_name",
            "inchikey",
            "ligand_chemotype",
            "smiles",
            "target_gene",
            "target_uniprot",
            "target_name",
            "atlas_score",
            "consensus_score",
            "free_cmax_um",
            "cmax_um",
            "fraction_unbound_plasma",
            "pk_source",
            "pk_missing_reason",
        )
        if col in out.columns
    ]
    out = out[cols].copy()
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    return out


def add_standard_pk_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "drug_id" not in out.columns:
        for source in ("ligand_base", "dedup_drug_key", "drugcentral_id", "mapped_drug_name"):
            if source in out.columns:
                out["drug_id"] = out[source]
                break
    if "free_cmax_um" not in out.columns:
        for source in ("free_cmax_uM", "free_cmax", "free_cmax_umol_l"):
            if source in out.columns:
                out["free_cmax_um"] = pd.to_numeric(out[source], errors="coerce")
                break
    if "cmax_um" not in out.columns:
        for source in ("cmax_uM", "total_cmax_um", "total_cmax_uM"):
            if source in out.columns:
                out["cmax_um"] = pd.to_numeric(out[source], errors="coerce")
                break
    if "fraction_unbound_plasma" not in out.columns:
        for source in ("fu", "fu_plasma", "fraction_unbound", "fraction_unbound_plasma"):
            if source in out.columns:
                out["fraction_unbound_plasma"] = pd.to_numeric(out[source], errors="coerce")
                break
    if "pk_source" not in out.columns:
        out["pk_source"] = "dedicated_pk_table"
    out["pk_missing_reason"] = out.get("pk_missing_reason", pd.Series("", index=out.index))
    return out


def _activity_type_from_papyrus(row: pd.Series) -> str:
    for col, label in (
        ("type_IC50", "IC50"),
        ("type_Ki", "Ki"),
        ("type_KD", "Kd"),
        ("type_EC50", "EC50"),
    ):
        if str(row.get(col, "")).split(";")[0] == "1":
            return label
    return "activity"


def build_papyrus_pilot_extract(
    papyrus_xz_path: str | Path,
    pilot_pair_table_path: str | Path,
    out_path: str | Path,
    *,
    chunksize: int = 250_000,
) -> pd.DataFrame:
    pair = pd.read_csv(pilot_pair_table_path)
    pair["inchikey_connectivity"] = pair["inchikey"].fillna("").astype(str).str.split("-").str[0]
    ligand_map = pair[["inchikey_connectivity", "drug_id"]].dropna().drop_duplicates()
    ligand_map = ligand_map[ligand_map["inchikey_connectivity"].astype(str).str.len() > 0]
    targets = {str(v).upper() for v in pair["target_id"].dropna().unique()}
    ligand_keys = set(ligand_map["inchikey_connectivity"].astype(str))
    pieces: list[pd.DataFrame] = []
    for chunk in pd.read_csv(papyrus_xz_path, sep="\t", compression="xz", chunksize=chunksize):
        if "accession" not in chunk.columns or "InChIKey" not in chunk.columns:
            continue
        conn = chunk["InChIKey"].fillna("").astype(str).str.split("-").str[0]
        keep = chunk["accession"].fillna("").astype(str).str.upper().isin(targets) & conn.isin(ligand_keys)
        if not keep.any():
            continue
        kept = chunk.loc[keep].copy()
        kept["inchikey_connectivity"] = conn.loc[keep].values
        kept = kept.merge(ligand_map, on="inchikey_connectivity", how="left")
        kept["target_id"] = kept["accession"].astype(str).str.upper()
        kept["activity_type"] = kept.apply(_activity_type_from_papyrus, axis=1)
        kept["pchembl_value"] = pd.to_numeric(kept.get("pchembl_value_Mean"), errors="coerce")
        kept["activity_relation"] = kept.get("relation", "")
        kept["activity_nM"] = pd.NA
        kept["activity_units"] = "nM"
        kept["assay_id"] = kept.get("AID", "")
        kept["source"] = kept.get("source", "Papyrus++")
        kept["activity_document_ids"] = kept.get("all_doc_ids", kept.get("doc_id", ""))
        kept["activity_publication_year"] = pd.to_numeric(kept.get("Year"), errors="coerce")
        kept["database_release_year"] = 2024
        pieces.append(
            kept[
                [
                    "drug_id",
                    "target_id",
                    "activity_nM",
                    "pchembl_value",
                    "activity_relation",
                    "activity_type",
                    "activity_units",
                    "assay_id",
                    "source",
                    "activity_document_ids",
                    "activity_publication_year",
                    "database_release_year",
                ]
            ]
        )
    out = pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame(
        columns=[
            "drug_id",
            "target_id",
            "activity_nM",
            "pchembl_value",
            "activity_relation",
            "activity_type",
            "activity_units",
            "assay_id",
            "source",
            "activity_document_ids",
            "activity_publication_year",
            "database_release_year",
        ]
    )
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, sep="\t", index=False)
    return out


def _chembl_molecule_names(row: pd.Series) -> list[str]:
    names = []
    for col in ("display_name", "generic_name", "ligand_display"):
        value = _text(row.get(col))
        if value and not value.lower().startswith("rdk_"):
            names.append(value)
    return list(dict.fromkeys(names))


def _name_key(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", " ", _text(value).lower()).strip()


def _cas_key(value: Any) -> str:
    return re.sub(r"[^0-9]+", "", _text(value))


def _first_existing(columns: list[str], aliases: tuple[str, ...]) -> str | None:
    def key(value: str) -> str:
        return re.sub(r"[^a-z0-9]+", "", value.lower())

    lowered = {key(col): col for col in columns}
    for alias in aliases:
        if key(alias) in lowered:
            return lowered[key(alias)]
    return None


def _read_zip_table(zip_path: str | Path, member_name: str) -> pd.DataFrame:
    with zipfile.ZipFile(zip_path) as archive:
        with archive.open(member_name) as handle:
            if member_name.lower().endswith(".xlsx"):
                return pd.read_excel(handle)
            return pd.read_csv(handle)


def _toxcast_member(zip_path: str | Path, contains: str, suffixes: tuple[str, ...]) -> str:
    contains_l = contains.lower()
    with zipfile.ZipFile(zip_path) as archive:
        matches = [
            name
            for name in archive.namelist()
            if contains_l in Path(name).name.lower() and Path(name).name.lower().endswith(suffixes)
        ]
    if not matches:
        raise FileNotFoundError(f"No ToxCast archive member matching {contains!r} and {suffixes!r}")
    csv_matches = [name for name in matches if name.lower().endswith(".csv")]
    return sorted(csv_matches or matches)[0]


def _build_toxcast_drug_map(pair: pd.DataFrame) -> tuple[dict[str, str], dict[str, str]]:
    name_map: dict[str, str] = {}
    cas_map: dict[str, str] = {}
    for _, row in pair.drop_duplicates("drug_id").iterrows():
        drug_id = _text(row.get("drug_id"))
        if not drug_id:
            continue
        for col in (
            "display_name",
            "generic_name",
            "pubchem_name",
            "pubchem_record_title",
            "sdf_title",
            "remark_name",
        ):
            key = _name_key(row.get(col))
            if key and not key.startswith("rdk "):
                name_map.setdefault(key, drug_id)
        for col in ("brand_names", "pubchem_synonyms", "rxnorm_brand_names", "drugcentral_brand_names"):
            for value in _text(row.get(col)).replace("|", ";").split(";"):
                key = _name_key(value)
                if key:
                    name_map.setdefault(key, drug_id)
        cas = _cas_key(row.get("cas"))
        if cas:
            cas_map.setdefault(cas, drug_id)
    return name_map, cas_map


def _build_toxcast_dsstox_map(pair: pd.DataFrame, zip_path: str | Path) -> dict[str, str]:
    member = _toxcast_member(zip_path, "analytical_qc", (".xlsx", ".csv"))
    qc = _read_zip_table(zip_path, member)
    name_to_dsstox: dict[str, str] = {}
    cas_to_dsstox: dict[str, str] = {}
    dsstox_col = _first_existing(list(qc.columns), ("dsstox_substance_id", "dtxsid", "dsstox_id"))
    name_col = _first_existing(list(qc.columns), ("chnm", "chemical_name", "name"))
    cas_col = _first_existing(list(qc.columns), ("casn", "casrn", "cas"))
    if not dsstox_col:
        return {}
    for _, row in qc.iterrows():
        dsstox_id = _text(row.get(dsstox_col))
        if not dsstox_id:
            continue
        if name_col:
            key = _name_key(row.get(name_col))
            if key:
                name_to_dsstox.setdefault(key, dsstox_id)
        if cas_col:
            key = _cas_key(row.get(cas_col))
            if key:
                cas_to_dsstox.setdefault(key, dsstox_id)
    mapping: dict[str, str] = {}
    for _, row in pair.drop_duplicates("drug_id").iterrows():
        drug_id = _text(row.get("drug_id"))
        if not drug_id:
            continue
        cas = _cas_key(row.get("cas"))
        if cas and cas in cas_to_dsstox:
            mapping[drug_id] = cas_to_dsstox[cas]
            continue
        for col in (
            "display_name",
            "generic_name",
            "pubchem_name",
            "pubchem_record_title",
            "sdf_title",
            "remark_name",
        ):
            key = _name_key(row.get(col))
            if key and key in name_to_dsstox:
                mapping[drug_id] = name_to_dsstox[key]
                break
    return mapping


def _load_toxcast_annotation_filters(zip_path: str | Path, cache_dir: Path) -> pd.DataFrame:
    member = _toxcast_member(zip_path, "assay_annotations", (".xlsx", ".csv"))
    raw = _read_zip_table(zip_path, member)
    aeid_col = _first_existing(list(raw.columns), ("aeid", "assay_endpoint_id"))
    if not aeid_col:
        raise ValueError("ToxCast assay annotations lack aeid")
    export_col = _first_existing(list(raw.columns), ("export_ready",))
    usability_col = _first_existing(list(raw.columns), ("data_usability",))
    viability_col = _first_existing(list(raw.columns), ("cell_viability_assay",))
    burst_col = _first_existing(list(raw.columns), ("burst_assay",))
    usable = pd.Series(True, index=raw.index)
    if export_col:
        usable &= pd.to_numeric(raw[export_col], errors="coerce").fillna(0).eq(1)
    if usability_col:
        usable &= pd.to_numeric(raw[usability_col], errors="coerce").fillna(0).eq(1)
    if viability_col:
        usable &= ~pd.to_numeric(raw[viability_col], errors="coerce").fillna(0).eq(1)
    out = pd.DataFrame(
        {
            "aeid": raw[aeid_col].astype(str),
            "toxcast_assay_usable": usable,
            "toxcast_export_ready": pd.to_numeric(raw[export_col], errors="coerce") if export_col else pd.NA,
            "toxcast_data_usability": pd.to_numeric(raw[usability_col], errors="coerce") if usability_col else pd.NA,
            "toxcast_cell_viability_assay": pd.to_numeric(raw[viability_col], errors="coerce") if viability_col else pd.NA,
            "toxcast_burst_assay": pd.to_numeric(raw[burst_col], errors="coerce") if burst_col else pd.NA,
        }
    ).drop_duplicates("aeid")
    out.to_csv(cache_dir / "toxcast_assay_annotation_filters.csv", index=False)
    return out


def _load_toxcast_chemical_qc(zip_path: str | Path, cache_dir: Path) -> pd.DataFrame:
    member = _toxcast_member(zip_path, "analytical_qc", (".xlsx", ".csv"))
    raw = _read_zip_table(zip_path, member)
    dsstox_col = _first_existing(list(raw.columns), ("dsstox_substance_id", "dtxsid", "dsstox_id"))
    pass_col = _first_existing(list(raw.columns), ("pass_or_caution", "qc_level"))
    if not dsstox_col:
        return pd.DataFrame(columns=["dsstox_substance_id", "toxcast_chemical_qc_status", "toxcast_chemical_qc_usable"])
    status = raw[pass_col].astype(str).str.lower() if pass_col else pd.Series("unknown", index=raw.index)
    usable = status.isin({"pass", "caution", "pass_or_caution"})
    out = pd.DataFrame(
        {
            "dsstox_substance_id": raw[dsstox_col].astype(str),
            "toxcast_chemical_qc_status": status,
            "toxcast_chemical_qc_usable": usable,
        }
    )
    out = (
        out.sort_values(["dsstox_substance_id", "toxcast_chemical_qc_usable"], ascending=[True, False])
        .drop_duplicates("dsstox_substance_id")
        .copy()
    )
    out.to_csv(cache_dir / "toxcast_chemical_qc.csv", index=False)
    return out


def _load_toxcast_cytotox(zip_path: str | Path, cache_dir: Path) -> pd.DataFrame:
    member = _toxcast_member(zip_path, "cytotox", (".xlsx", ".csv"))
    raw = _read_zip_table(zip_path, member)
    dsstox_col = _first_existing(list(raw.columns), ("dsstox_substance_id", "dtxsid", "dsstox_id"))
    lower_col = _first_existing(list(raw.columns), ("cytotox_lower_bound_um", "cytotox_lower_bound"))
    median_col = _first_existing(list(raw.columns), ("cytotox_median_um", "cytotox_median"))
    if not dsstox_col or not lower_col:
        return pd.DataFrame(columns=["dsstox_substance_id", "cytotox_lower_bound_um", "cytotox_median_um"])
    out = pd.DataFrame(
        {
            "dsstox_substance_id": raw[dsstox_col].astype(str),
            "cytotox_lower_bound_um": pd.to_numeric(raw[lower_col], errors="coerce"),
            "cytotox_median_um": pd.to_numeric(raw[median_col], errors="coerce") if median_col else pd.NA,
        }
    ).drop_duplicates("dsstox_substance_id")
    out.to_csv(cache_dir / "toxcast_cytotox_thresholds.csv", index=False)
    return out


def _build_toxcast_endpoint_map(zip_path: str | Path, target_lookup: dict[str, str], out_path: Path) -> pd.DataFrame:
    member = _toxcast_member(zip_path, "assay_target_mappings", (".xlsx", ".csv"))
    raw = _read_zip_table(zip_path, member)
    endpoint_col = _first_existing(
        list(raw.columns),
        (
            "aenm",
            "assay_component_endpoint_name",
            "assay_endpoint_name",
            "endpoint_name",
            "assay",
        ),
    )
    aeid_col = _first_existing(list(raw.columns), ("aeid", "assay_endpoint_id"))
    uniprot_col = _first_existing(
        list(raw.columns),
        (
            "intended_target_uniprot_accession_number",
            "target_uniprot_accession_number",
            "uniprot_accession_number",
            "uniprot",
            "accession",
            "target_id",
        ),
    )
    gene_col = _first_existing(
        list(raw.columns),
        (
            "intended_target_gene_symbol",
            "target_gene_symbol",
            "gene_symbol",
            "official_symbol",
            "symbol",
            "gene",
        ),
    )
    if not endpoint_col and not aeid_col:
        raise ValueError("ToxCast assay target mapping lacks an endpoint name or aeid column")
    if not uniprot_col and not gene_col:
        raise ValueError("ToxCast assay target mapping lacks UniProt or gene target columns")
    rows = []
    for _, row in raw.iterrows():
        targets: list[str] = []
        if uniprot_col:
            targets.extend(str(row.get(uniprot_col, "")).replace("|", ";").replace(",", ";").split(";"))
        if gene_col:
            targets.extend(str(row.get(gene_col, "")).replace("|", ";").replace(",", ";").split(";"))
        for target in targets:
            target_key = _text(target).upper()
            target_id = target_lookup.get(target_key, "")
            if not target_key or target_key in {"NAN", "NONE"} or not target_id:
                continue
            rows.append(
                {
                    "aenm": _text(row.get(endpoint_col)) if endpoint_col else "",
                    "aeid": _text(row.get(aeid_col)) if aeid_col else "",
                    "target_id": target_id,
                }
            )
    endpoint_map = pd.DataFrame(rows).drop_duplicates() if rows else pd.DataFrame(columns=["aenm", "aeid", "target_id"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    endpoint_map.to_csv(out_path, index=False)
    return endpoint_map


def _chembl_lookup_molecules(pair: pd.DataFrame, cache_path: Path) -> dict[str, str]:
    if cache_path.exists():
        return json.loads(cache_path.read_text(encoding="utf-8"))
    approved_name_to_chembl: dict[str, str] = {}
    next_url = f"{CHEMBL_API}/molecule.json?max_phase=4&limit=1000"
    while next_url:
        data = _http_json(next_url if next_url.startswith("http") else f"https://www.ebi.ac.uk{next_url}")
        for molecule in data.get("molecules", []):
            chembl_id = molecule.get("molecule_chembl_id")
            if not chembl_id:
                continue
            for name in [molecule.get("pref_name"), *[syn.get("molecule_synonym") for syn in molecule.get("molecule_synonyms", [])]]:
                key = _name_key(name)
                if key:
                    approved_name_to_chembl.setdefault(key, chembl_id)
        next_url = data.get("page_meta", {}).get("next")
    mapping: dict[str, str] = {}
    ligand_rows = pair.drop_duplicates("drug_id")
    for _, row in ligand_rows.iterrows():
        drug_id = _text(row.get("drug_id"))
        for name in _chembl_molecule_names(row):
            chembl_id = approved_name_to_chembl.get(_name_key(name))
            if chembl_id:
                mapping[drug_id] = chembl_id
                break
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(mapping, indent=2, sort_keys=True), encoding="utf-8")
    return mapping


def _chembl_lookup_targets(pair: pd.DataFrame, cache_path: Path) -> dict[str, str]:
    if cache_path.exists():
        return json.loads(cache_path.read_text(encoding="utf-8"))
    mapping: dict[str, str] = {}
    for target_id in sorted({str(v).upper() for v in pair["target_id"].dropna().unique()}):
        if not target_id:
            continue
        query = urllib.parse.urlencode({"target_components__accession": target_id, "limit": 20})
        data = _http_json(f"{CHEMBL_API}/target.json?{query}")
        targets = [
            item
            for item in data.get("targets", [])
            if item.get("organism") == "Homo sapiens" and item.get("target_type") in {"SINGLE PROTEIN", "PROTEIN FAMILY"}
        ]
        if targets:
            mapping[target_id] = targets[0]["target_chembl_id"]
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(mapping, indent=2, sort_keys=True), encoding="utf-8")
    return mapping


def build_chembl_pilot_extract(
    pilot_pair_table_path: str | Path,
    out_path: str | Path,
    cache_dir: str | Path,
) -> pd.DataFrame:
    pair = pd.read_csv(pilot_pair_table_path)
    cache = Path(cache_dir)
    molecule_map = _chembl_lookup_molecules(pair, cache / "chembl_molecule_map.json")
    target_map = _chembl_lookup_targets(pair, cache / "chembl_target_map.json")
    reverse_mol = {value: key for key, value in molecule_map.items()}
    reverse_target = {value: key for key, value in target_map.items()}
    rows: list[dict[str, Any]] = []
    molecule_ids = sorted(set(molecule_map.values()))
    for target_chembl, target_id in reverse_target.items():
        for start in range(0, len(molecule_ids), 200):
            chunk = molecule_ids[start : start + 200]
            chunk_cache = cache / "chembl_pair_activity" / f"{target_chembl}_{start}.jsonl"
            if chunk_cache.exists():
                activities = [
                    json.loads(line)
                    for line in chunk_cache.read_text(encoding="utf-8").splitlines()
                    if line.strip()
                ]
            else:
                activities = []
                query = urllib.parse.urlencode(
                    {
                        "target_chembl_id": target_chembl,
                        "molecule_chembl_id__in": ",".join(chunk),
                        "standard_type__in": "IC50,Ki,Kd,EC50",
                        "limit": 1000,
                    }
                )
                next_url = f"{CHEMBL_API}/activity.json?{query}"
                while next_url:
                    data = _http_json(
                        next_url if next_url.startswith("http") else f"https://www.ebi.ac.uk{next_url}",
                        sleep_sec=0.01,
                    )
                    activities.extend(data.get("activities", []))
                    next_url = data.get("page_meta", {}).get("next")
                chunk_cache.parent.mkdir(parents=True, exist_ok=True)
                chunk_cache.write_text(
                    "\n".join(json.dumps(item, sort_keys=True) for item in activities),
                    encoding="utf-8",
                )
            for item in activities:
                mol = item.get("molecule_chembl_id")
                rows.append(
                    {
                        "drug_id": reverse_mol.get(mol, mol),
                        "target_id": target_id,
                        "activity_nM": item.get("standard_value"),
                        "pchembl_value": item.get("pchembl_value"),
                        "activity_relation": item.get("standard_relation"),
                        "activity_type": item.get("standard_type"),
                        "activity_units": item.get("standard_units"),
                        "assay_id": item.get("assay_chembl_id"),
                        "source": "ChEMBL",
                        "activity_document_ids": item.get("document_chembl_id"),
                        "activity_publication_year": item.get("document_year"),
                        "database_release_year": 2026,
                        "chembl_src_id": item.get("src_id"),
                        "chembl_data_validity_comment": item.get("data_validity_comment"),
                        "chembl_potential_duplicate": item.get("potential_duplicate"),
                    }
                )
    out = pd.DataFrame(rows)
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, sep="\t", index=False)
    (cache / "chembl_extract_summary.json").write_text(
        json.dumps(
            {
                "mapped_molecules": len(molecule_map),
                "mapped_targets": len(target_map),
                "activity_rows": len(out),
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return out


def build_toxcast_pilot_extract(
    toxcast_zip_path: str | Path,
    pilot_pair_table_path: str | Path,
    fda_mapping_path: str | Path,
    out_path: str | Path,
    cache_dir: str | Path,
    *,
    active_threshold_nM: float = 10000.0,
    chunksize: int = 250_000,
) -> pd.DataFrame:
    pair = pd.read_csv(pilot_pair_table_path)
    fda_mapping = pd.read_csv(fda_mapping_path)
    if "ligand_base" not in fda_mapping.columns:
        fda_mapping["ligand_base"] = fda_mapping.apply(
            lambda row: f"{_text(row.get('scheme'))}_{int(float(row.get('file_num'))):07d}"
            if _text(row.get("scheme")) and _text(row.get("file_num"))
            else "",
            axis=1,
        )
    drug_meta = pair[["drug_id"]].drop_duplicates().merge(
        fda_mapping.drop_duplicates("ligand_base"),
        left_on="drug_id",
        right_on="ligand_base",
        how="left",
    )
    name_map, cas_map = _build_toxcast_drug_map(drug_meta)
    dsstox_map = _build_toxcast_dsstox_map(drug_meta, toxcast_zip_path)
    target_lookup: dict[str, str] = {}
    for _, row in pair.drop_duplicates("target_id").iterrows():
        target_id = _text(row.get("target_id")).upper()
        if not target_id:
            continue
        for col in ("target_id", "target_uniprot", "target_gene"):
            key = _text(row.get(col)).upper()
            if key:
                target_lookup.setdefault(key, target_id)
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    annotation_filters = _load_toxcast_annotation_filters(toxcast_zip_path, cache)
    assay_usable = annotation_filters.set_index("aeid")["toxcast_assay_usable"].to_dict()
    chemical_qc = _load_toxcast_chemical_qc(toxcast_zip_path, cache)
    chemical_qc_usable = chemical_qc.set_index("dsstox_substance_id")["toxcast_chemical_qc_usable"].to_dict()
    cytotox = _load_toxcast_cytotox(toxcast_zip_path, cache)
    cytotox_lower = cytotox.set_index("dsstox_substance_id")["cytotox_lower_bound_um"].to_dict()
    endpoint_map = _build_toxcast_endpoint_map(
        toxcast_zip_path,
        target_lookup,
        cache / "toxcast_endpoint_target_map.csv",
    )
    if endpoint_map.empty:
        out = pd.DataFrame(
            columns=[
                "drug_id",
                "target_id",
                "activity_nM",
                "pchembl_value",
                "activity_relation",
                "activity_type",
                "activity_units",
                "assay_id",
                "source",
                "toxcast_hit_call",
            ]
        )
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(out_path, sep="\t", index=False)
        return out

    endpoint_to_targets = endpoint_map.groupby("aenm")["target_id"].apply(lambda s: sorted(set(s))).to_dict()
    aeid_to_targets = endpoint_map.groupby("aeid")["target_id"].apply(lambda s: sorted(set(s))).to_dict()
    member = _toxcast_member(toxcast_zip_path, "mc5-6_winning_model_fits", (".csv",))
    rows: list[pd.DataFrame] = []
    with zipfile.ZipFile(toxcast_zip_path) as archive:
        with archive.open(member) as handle:
            for chunk in pd.read_csv(handle, chunksize=chunksize, low_memory=False):
                cols = list(chunk.columns)
                aenm_col = _first_existing(cols, ("aenm", "assay_component_endpoint_name", "assay_endpoint_name"))
                aeid_col = _first_existing(cols, ("aeid", "assay_endpoint_id"))
                hit_col = _first_existing(cols, ("hitc", "hit_call", "hitcall"))
                modl_ga_col = _first_existing(cols, ("modl_ga", "gain_ac50", "ac50"))
                chemical_col = _first_existing(cols, ("chnm", "chemical_name", "name", "preferred_name"))
                cas_col = _first_existing(cols, ("casn", "casrn", "cas", "cas_number"))
                dsstox_col = _first_existing(cols, ("dsstox_substance_id", "dtxsid", "dsstox_id"))
                if not hit_col or (not aenm_col and not aeid_col) or (not chemical_col and not cas_col):
                    raise ValueError(
                        "ToxCast mc5-6 file lacks required hit, endpoint, or chemical identifier columns"
                    )
                drug_by_name = (
                    chunk[chemical_col].map(lambda value: name_map.get(_name_key(value)))
                    if chemical_col
                    else pd.Series(pd.NA, index=chunk.index)
                )
                drug_by_cas = (
                    chunk[cas_col].map(lambda value: cas_map.get(_cas_key(value)))
                    if cas_col
                    else pd.Series(pd.NA, index=chunk.index)
                )
                drug_id = drug_by_name.fillna(drug_by_cas)
                drug_dsstox = drug_id.map(lambda value: dsstox_map.get(_text(value)))
                if dsstox_col:
                    row_dsstox = chunk[dsstox_col].map(_text)
                    drug_dsstox = drug_dsstox.fillna(row_dsstox)
                if aenm_col:
                    target_lists = chunk[aenm_col].map(lambda value: endpoint_to_targets.get(_text(value), []))
                else:
                    target_lists = pd.Series([[] for _ in range(len(chunk))], index=chunk.index)
                if aeid_col:
                    target_lists = target_lists.where(
                        target_lists.map(bool),
                        chunk[aeid_col].map(lambda value: aeid_to_targets.get(_text(value), [])),
                    )
                keep = drug_id.notna() & target_lists.map(bool)
                if not keep.any():
                    continue
                kept = chunk.loc[keep].copy()
                kept["_drug_id"] = drug_id.loc[keep].values
                kept["_dsstox_substance_id"] = drug_dsstox.loc[keep].values
                kept["_targets"] = target_lists.loc[keep].values
                kept = kept.explode("_targets")
                hit = pd.to_numeric(kept[hit_col], errors="coerce")
                modl_ga = pd.to_numeric(kept[modl_ga_col], errors="coerce") if modl_ga_col else pd.Series(pd.NA, index=kept.index)
                if modl_ga_col and modl_ga_col.lower() == "ac50":
                    ac50_nM = modl_ga * 1000.0
                else:
                    # Older tcpl matrix fields store model-fit concentrations as log10 micromolar.
                    ac50_nM = (10.0**modl_ga) * 1000.0
                active = hit.eq(1)
                inactive = hit.eq(0)
                assay_id = kept[aeid_col].astype(str) if aeid_col else kept[aenm_col].astype(str)
                assay_ok = assay_id.map(lambda value: bool(assay_usable.get(str(value), False)))
                chemical_ok = kept["_dsstox_substance_id"].map(
                    lambda value: bool(chemical_qc_usable.get(_text(value), True))
                )
                cytotox_lower_um = pd.to_numeric(
                    kept["_dsstox_substance_id"].map(lambda value: cytotox_lower.get(_text(value))),
                    errors="coerce",
                )
                cytotox_confounded = active & cytotox_lower_um.notna() & (ac50_nM >= (cytotox_lower_um * 1000.0))
                usable_label = assay_ok & chemical_ok & ~cytotox_confounded
                activity_nM = ac50_nM.where(active & ac50_nM.notna(), active_threshold_nM + 1.0)
                activity_relation = pd.Series(">", index=kept.index)
                activity_relation.loc[active] = "="
                out_chunk = pd.DataFrame(
                    {
                        "drug_id": kept["_drug_id"].astype(str),
                        "target_id": kept["_targets"].astype(str).str.upper(),
                        "toxcast_dsstox_substance_id": kept["_dsstox_substance_id"].astype(str),
                        "activity_nM": activity_nM,
                        "pchembl_value": pd.NA,
                        "activity_relation": activity_relation,
                        "activity_type": "ToxCast hit_call",
                        "activity_units": "nM",
                        "assay_id": assay_id,
                        "source": "ToxCast invitrodb v4.3",
                        "toxcast_source_role": "in_vitro_activity_benchmark_not_clinical_adr_validation",
                        "toxcast_hit_call": hit,
                        "toxcast_assay_quality_pass": assay_ok,
                        "toxcast_chemical_qc_pass": chemical_ok,
                        "toxcast_cytotox_lower_bound_um": cytotox_lower_um,
                        "toxcast_cytotoxicity_confounded": cytotox_confounded,
                        "toxcast_label_usable": usable_label,
                        "database_release_year": 2025,
                    }
                )
                rows.append(out_chunk.loc[(active | inactive) & usable_label])
    out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=[
            "drug_id",
            "target_id",
            "activity_nM",
            "pchembl_value",
            "activity_relation",
            "activity_type",
            "activity_units",
            "assay_id",
            "source",
            "toxcast_source_role",
            "toxcast_hit_call",
            "toxcast_assay_quality_pass",
            "toxcast_chemical_qc_pass",
            "toxcast_cytotox_lower_bound_um",
            "toxcast_cytotoxicity_confounded",
            "toxcast_label_usable",
            "database_release_year",
        ]
    )
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, sep="\t", index=False)
    (cache / "toxcast_extract_summary.json").write_text(
        json.dumps(
            {
                "archive": str(toxcast_zip_path),
                "archive_sha256_expected": TOXCAST_SUMMARY_SHA256,
                "mc56_member": member,
                "mapped_name_keys": len(name_map),
                "mapped_cas_keys": len(cas_map),
                "mapped_dsstox_ids": len(dsstox_map),
                "mapped_target_keys": len(target_lookup),
                "mapped_endpoint_targets": len(endpoint_map),
                "usable_assay_endpoints": int(pd.Series(assay_usable).eq(True).sum()),
                "chemical_qc_records": len(chemical_qc),
                "cytotox_records": len(cytotox),
                "activity_rows": len(out),
                "active_rows": int(out["toxcast_hit_call"].eq(1).sum()) if "toxcast_hit_call" in out else 0,
                "inactive_rows": int(out["toxcast_hit_call"].eq(0).sum()) if "toxcast_hit_call" in out else 0,
                "filters": {
                    "assay": "export_ready==1 and data_usability==1 and cell_viability_assay!=1",
                    "chemical_qc": "analytical_qc pass_or_caution in {pass,caution}",
                    "cytotoxicity": "exclude active hits with ac50_uM >= cytotox_lower_bound_um",
                },
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return out


def write_holdout_summaries(dataset_path: str | Path, out_dir: str | Path) -> pd.DataFrame:
    df = pd.read_csv(dataset_path)
    rows: list[dict[str, Any]] = []
    for split in ("drug_holdout", "target_holdout", "scaffold_holdout"):
        try:
            train_idx, test_idx = make_split(df, split_mode=split, seed=42)
            summary = split_overlap_summary(df.loc[train_idx], df.loc[test_idx], split)
            rows.append({**summary, "overlaps": json.dumps(summary["overlaps"], sort_keys=True)})
        except Exception as exc:
            rows.append({"split_mode": split, "n_train": 0, "n_test": 0, "overlaps": "{}", "passes_holdout": False, "error": str(exc)})
    out = pd.DataFrame(rows)
    path = Path(out_dir) / "ml_holdout_readiness_summary.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    return out


def _source_flag(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(False, index=df.index)
    return df[col].notna()


def _true_flag(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(False, index=df.index)
    return df[col].eq(True).fillna(False)


def _year_value(value: Any) -> float | None:
    year = pd.to_numeric(value, errors="coerce")
    if pd.isna(year):
        return None
    return float(year)


def _assign_label_provenance(labeled: pd.DataFrame) -> pd.DataFrame:
    out = labeled.copy()
    source_flags = {
        "chembl": _source_flag(out, "chembl_active") | _source_flag(out, "chembl_inactive"),
        "papyrus": _source_flag(out, "papyrus_active") | _source_flag(out, "papyrus_inactive"),
        "toxcast": _source_flag(out, "toxcast_active") | _source_flag(out, "toxcast_inactive"),
    }
    for source, flags in source_flags.items():
        out[f"{source}_label_available"] = flags
    provenance: list[str] = []
    release_years: list[float | None] = []
    publication_years: list[float | None] = []
    evidence_year_sources: list[str] = []
    for idx in out.index:
        sources = [source for source, flags in source_flags.items() if bool(flags.loc[idx])]
        provenance.append(";".join(sources) if sources else "unlabeled")
        row_release_years: list[float] = []
        row_publication_years: list[tuple[str, float]] = []
        if "chembl" in sources:
            row_release_years.append(
                _year_value(out.get("chembl_database_release_year", pd.Series(pd.NA, index=out.index)).loc[idx])
                or 2026.0
            )
            year = _year_value(out.get("chembl_activity_publication_year_min", pd.Series(pd.NA, index=out.index)).loc[idx])
            if year is not None:
                row_publication_years.append(("chembl", year))
        if "papyrus" in sources:
            row_release_years.append(
                _year_value(out.get("papyrus_database_release_year", pd.Series(pd.NA, index=out.index)).loc[idx])
                or 2024.0
            )
            year = _year_value(out.get("papyrus_activity_publication_year_min", pd.Series(pd.NA, index=out.index)).loc[idx])
            if year is not None:
                row_publication_years.append(("papyrus", year))
        if "toxcast" in sources:
            row_release_years.append(
                _year_value(out.get("toxcast_database_release_year", pd.Series(pd.NA, index=out.index)).loc[idx])
                or 2025.0
            )
        release_years.append(max(row_release_years) if row_release_years else None)
        if row_publication_years:
            source, year = min(row_publication_years, key=lambda item: item[1])
            publication_years.append(year)
            evidence_year_sources.append(source)
        else:
            publication_years.append(None)
            evidence_year_sources.append("missing_publication_year")
    out["label_source"] = provenance
    out["label_source_count"] = [0 if value == "unlabeled" else len(value.split(";")) for value in provenance]
    out["database_release_year"] = release_years
    out["activity_publication_year"] = publication_years
    out["activity_publication_year_source"] = evidence_year_sources
    return out


def stage_pilot_bioactivity_sources(
    pilot_dir: str | Path,
    repo_root: str | Path,
    *,
    external_dir: str | Path | None = None,
    fda_mapping_path: str | Path | None = None,
    pk_table_path: str | Path | None = None,
    overwrite_downloads: bool = False,
) -> dict[str, str]:
    pilot = Path(pilot_dir)
    root = Path(repo_root)
    external = Path(external_dir) if external_dir is not None else root / "data" / "external"
    resolved_fda_mapping_path = Path(fda_mapping_path) if fda_mapping_path is not None else resolve_fda_mapping_path(root)
    pair_table = build_pilot_pair_table(
        pilot / "master_rows.csv",
        resolved_fda_mapping_path,
        root / "input_pdbs",
        pilot / "bioactivity_pair_table.csv",
        pk_table_path=pk_table_path,
    )
    papyrus_xz = _download(
        PAPYRUS_PLUS_URL,
        external / "papyrus" / "05.7++_combined_set_without_stereochemistry.tsv.xz",
        overwrite=overwrite_downloads,
    )
    papyrus_source = build_papyrus_pilot_extract(
        papyrus_xz,
        pilot / "bioactivity_pair_table.csv",
        external / "papyrus" / "bioactivity.tsv",
    )
    chembl_source = build_chembl_pilot_extract(
        pilot / "bioactivity_pair_table.csv",
        external / "chembl" / "bioactivity.tsv",
        pilot / "external_query_cache",
    )
    toxcast_zip = _download(
        TOXCAST_SUMMARY_URL,
        external / "toxcast" / "INVITRODB_SUMMARY.zip",
        overwrite=overwrite_downloads,
    )
    toxcast_source = build_toxcast_pilot_extract(
        toxcast_zip,
        pilot / "bioactivity_pair_table.csv",
        resolved_fda_mapping_path,
        external / "toxcast" / "bioactivity.tsv",
        pilot / "external_query_cache",
    )
    papyrus_benchmark = build_bioactivity_benchmark(
        pilot / "bioactivity_pair_table.csv",
        external / "papyrus" / "bioactivity.tsv",
        None,
        pilot / "papyrus_atlas_benchmark.csv",
        label_prefix="papyrus",
    )
    chembl_benchmark = build_bioactivity_benchmark(
        pilot / "bioactivity_pair_table.csv",
        external / "chembl" / "bioactivity.tsv",
        None,
        pilot / "chembl_atlas_benchmark.csv",
        label_prefix="chembl",
    )
    toxcast_benchmark = build_bioactivity_benchmark(
        pilot / "bioactivity_pair_table.csv",
        external / "toxcast" / "bioactivity.tsv",
        None,
        pilot / "toxcast_atlas_benchmark.csv",
        label_prefix="toxcast",
    )
    labeled = chembl_benchmark.merge(
        papyrus_benchmark[
            [col for col in papyrus_benchmark.columns if col in {"drug_id", "target_id"} or col.startswith("papyrus_")]
        ].drop_duplicates(["drug_id", "target_id"]),
        on=["drug_id", "target_id"],
        how="left",
    )
    labeled = labeled.merge(
        toxcast_benchmark[
            [col for col in toxcast_benchmark.columns if col in {"drug_id", "target_id"} or col.startswith("toxcast_")]
        ].drop_duplicates(["drug_id", "target_id"]),
        on=["drug_id", "target_id"],
        how="left",
    )
    labeled = _assign_label_provenance(labeled)
    labeled["bioactivity_ml_label"] = pd.NA
    active = _true_flag(labeled, "chembl_active") | _true_flag(labeled, "papyrus_active") | _true_flag(labeled, "toxcast_active")
    inactive = (
        _true_flag(labeled, "chembl_inactive")
        | _true_flag(labeled, "papyrus_inactive")
        | _true_flag(labeled, "toxcast_inactive")
    )
    conflicting = (
        _true_flag(labeled, "chembl_conflict")
        | _true_flag(labeled, "papyrus_conflict")
        | _true_flag(labeled, "toxcast_conflict")
    )
    cross_source_disagreement = active & inactive & ~conflicting
    ambiguous = conflicting | cross_source_disagreement
    active = active & ~ambiguous
    inactive = inactive & ~ambiguous
    labeled.loc[inactive, "bioactivity_ml_label"] = 0
    labeled.loc[active, "bioactivity_ml_label"] = 1
    labeled.loc[conflicting, "bioactivity_ml_label_missing_reason"] = "conflicting_active_inactive_evidence"
    labeled.loc[cross_source_disagreement, "bioactivity_ml_label_missing_reason"] = "cross_source_active_inactive_disagreement"
    labeled.loc[cross_source_disagreement, "bioactivity_ml_cross_source_disagreement"] = True
    ml_source = pilot / "bioactivity_ml_source.csv"
    labeled.to_csv(ml_source, index=False)
    ml_table = pilot / "ml_pair_table_bioactivity_dedup_leakage_controlled.csv"
    build_ml_dataset(
        ml_source,
        "bioactivity_ml_label",
        "pilot_nonleaky",
        ml_table,
        deduplicate_drug_target=True,
    )
    holdout = write_holdout_summaries(ml_table, pilot)
    summary = {
        "bioactivity_pair_table_rows": len(pair_table),
        "fda_mapping_path": str(resolved_fda_mapping_path),
        "free_cmax_nonmissing_rows": int(pd.to_numeric(pair_table.get("free_cmax_um"), errors="coerce").notna().sum()),
        "free_cmax_missing_rows": int(pd.to_numeric(pair_table.get("free_cmax_um"), errors="coerce").isna().sum()),
        "free_cmax_missing_policy": "unknown_not_negative; no exposure imputation applied",
        "pk_table_path": str(pk_table_path) if pk_table_path else "",
        "pk_join_policy": "dedicated_pk_table_only; no exposure inference",
        "papyrus_source_rows": len(papyrus_source),
        "chembl_source_rows": len(chembl_source),
        "toxcast_source_rows": len(toxcast_source),
        "papyrus_labeled_pairs": int(papyrus_benchmark.get("papyrus_active", pd.Series(dtype=object)).notna().sum()),
        "chembl_labeled_pairs": int(chembl_benchmark.get("chembl_active", pd.Series(dtype=object)).notna().sum()),
        "toxcast_labeled_pairs": int(toxcast_benchmark.get("toxcast_active", pd.Series(dtype=object)).notna().sum()),
        "bioactivity_conflicting_pairs": int(conflicting.sum()),
        "bioactivity_cross_source_disagreement_pairs": int(cross_source_disagreement.sum()),
        "bioactivity_ml_rows": int(pd.read_csv(ml_table).shape[0]) if ml_table.exists() else 0,
        "holdout_summaries": holdout.to_dict(orient="records"),
        "toxcast_status": "ingested_from_full_invitrodb_v4_3_summary_archive",
        "toxcast_role": "in_vitro_activity_benchmark_not_clinical_adr_validation",
        "toxcast_filters": "export_ready/data_usability/cell_viability/chemical_qc/cytotoxicity_lower_bound",
        "toxcast_archive": str(toxcast_zip),
    }
    summary_path = pilot / "external_bioactivity_stage_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "pair_table": str(pilot / "bioactivity_pair_table.csv"),
        "papyrus_source": str(external / "papyrus" / "bioactivity.tsv"),
        "chembl_source": str(external / "chembl" / "bioactivity.tsv"),
        "toxcast_source": str(external / "toxcast" / "bioactivity.tsv"),
        "papyrus_benchmark": str(pilot / "papyrus_atlas_benchmark.csv"),
        "chembl_benchmark": str(pilot / "chembl_atlas_benchmark.csv"),
        "toxcast_benchmark": str(pilot / "toxcast_atlas_benchmark.csv"),
        "bioactivity_ml_table": str(ml_table),
        "holdout_summary": str(pilot / "ml_holdout_readiness_summary.csv"),
        "summary": str(summary_path),
    }

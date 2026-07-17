from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.external.compound_index import load_atlas_compound_index, map_by_inchikey


PUG_REST = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

PUBCHEM_BIOASSAY_COLUMNS = [
    "drug_id",
    "drug_name",
    "inchikey",
    "smiles",
    "target_id",
    "accession",
    "assay_id",
    "activity_outcome",
    "activity_type",
    "activity_relation",
    "activity_units",
    "source",
    "pubchem_cid",
]


def _cid_properties_url(cids: list[int]) -> str:
    joined = ",".join(str(int(cid)) for cid in cids)
    return f"{PUG_REST}/compound/cid/{joined}/property/Title,InChIKey,CanonicalSMILES/JSON"


def _json(url: str) -> dict[str, Any]:
    request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "AtlasAnalysis/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode("utf-8"))


def _list_from_identifier_list(payload: dict[str, Any], key: str) -> list[Any]:
    return list(payload.get("IdentifierList", {}).get(key, []))


def aids_for_uniprot(accession: str) -> list[int]:
    url = f"{PUG_REST}/assay/target/accession/{urllib.parse.quote(str(accession))}/aids/JSON"
    try:
        return [int(v) for v in _list_from_identifier_list(_json(url), "AID")]
    except Exception:
        return []


def cids_for_aid(aid: int, cids_type: str) -> set[int]:
    query = urllib.parse.urlencode({"cids_type": cids_type})
    url = f"{PUG_REST}/assay/aid/{int(aid)}/cids/JSON?{query}"
    try:
        return {int(v) for v in _list_from_identifier_list(_json(url), "CID")}
    except Exception:
        return set()


def properties_for_cids(cids: list[int]) -> pd.DataFrame:
    if not cids:
        return pd.DataFrame(columns=["pubchem_cid", "pubchem_title", "inchikey", "smiles"])
    rows: list[dict[str, Any]] = []
    for start in range(0, len(cids), 100):
        chunk = sorted({int(cid) for cid in cids[start : start + 100]})
        try:
            payload = _json(_cid_properties_url(chunk))
        except Exception:
            continue
        for item in payload.get("PropertyTable", {}).get("Properties", []):
            rows.append(
                {
                    "pubchem_cid": item.get("CID"),
                    "inchikey": item.get("InChIKey"),
                    "pubchem_title": item.get("Title"),
                    "smiles": item.get("CanonicalSMILES")
                    or item.get("ConnectivitySMILES")
                    or item.get("SMILES"),
                }
            )
        time.sleep(0.1)
    return pd.DataFrame(rows).drop_duplicates()


def _read_pubchem_aid_csv(path: str | Path) -> pd.DataFrame:
    """Read PubChem AID full or concise CSV into a row-level table."""

    raw = pd.read_csv(path, low_memory=False)
    if "PUBCHEM_RESULT_TAG" in raw.columns:
        tag = pd.to_numeric(raw["PUBCHEM_RESULT_TAG"], errors="coerce")
        raw = raw[tag.notna()].copy()
        return pd.DataFrame(
            {
                "pubchem_cid": pd.to_numeric(raw.get("PUBCHEM_CID"), errors="coerce"),
                "smiles": raw.get("PUBCHEM_EXT_DATASOURCE_SMILES"),
                "activity_outcome": raw.get("PUBCHEM_ACTIVITY_OUTCOME"),
                "activity_value_nM": pd.to_numeric(raw.get("Standard Value"), errors="coerce"),
                "activity_relation": raw.get("Standard Relation"),
                "activity_type": raw.get("Standard Type"),
                "activity_units": raw.get("Standard Units"),
                "sid": raw.get("PUBCHEM_SID"),
            }
        )
    return pd.DataFrame(
        {
            "pubchem_cid": pd.to_numeric(raw.get("CID"), errors="coerce"),
            "smiles": pd.NA,
            "activity_outcome": raw.get("Activity Outcome"),
            "activity_value_nM": pd.to_numeric(raw.get("Activity Value [uM]"), errors="coerce") * 1000.0,
            "activity_relation": pd.NA,
            "activity_type": raw.get("Activity Name"),
            "activity_units": "nM",
            "sid": raw.get("SID"),
            "target_id": raw.get("Target Accession"),
            "assay_name": raw.get("Assay Name"),
            "assay_type": raw.get("Assay Type"),
            "pubmed_id": raw.get("PubMed ID"),
        }
    )


def _target_accessions_from_summary(path: str | Path | None) -> list[str]:
    if path is None or not Path(path).exists():
        return []
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    accessions: list[str] = []
    for summary in payload.get("AssaySummaries", {}).get("AssaySummary", []):
        for target in summary.get("Target", []) or []:
            accession = str(target.get("Accession") or "").strip()
            if accession:
                accessions.append(accession)
    return sorted(set(accessions))


def stage_pubchem_aid_table(
    aid_csv_path: str | Path,
    mapping_path: str | Path,
    out_path: str | Path,
    *,
    aid: int,
    summary_json_path: str | Path | None = None,
    cid_properties_cache: str | Path | None = None,
    keep_unmapped: bool = True,
) -> pd.DataFrame:
    """Normalize one downloaded PubChem AID table into Atlas source-table schema."""

    aid_rows = _read_pubchem_aid_csv(aid_csv_path)
    aid_rows = aid_rows[aid_rows["pubchem_cid"].notna()].copy()
    aid_rows["pubchem_cid"] = aid_rows["pubchem_cid"].astype(int)
    target_accessions = _target_accessions_from_summary(summary_json_path)
    if "target_id" not in aid_rows.columns or aid_rows["target_id"].isna().all():
        aid_rows["target_id"] = target_accessions[0] if target_accessions else pd.NA
    properties_cache_path = (
        Path(cid_properties_cache)
        if cid_properties_cache is not None
        else Path(out_path).with_suffix(".cid_properties.csv")
    )
    if properties_cache_path.exists():
        properties = pd.read_csv(properties_cache_path, low_memory=False)
    else:
        properties = properties_for_cids(sorted(aid_rows["pubchem_cid"].dropna().astype(int).unique()))
        properties_cache_path.parent.mkdir(parents=True, exist_ok=True)
        properties.to_csv(properties_cache_path, index=False)
    aid_rows = aid_rows.merge(properties, on="pubchem_cid", how="left", suffixes=("", "_resolved"))
    aid_rows["smiles"] = aid_rows["smiles"].fillna(aid_rows.get("smiles_resolved"))
    atlas_index = load_atlas_compound_index(mapping_path)
    mapped = map_by_inchikey(aid_rows, atlas_index, source_inchikey_col="inchikey")
    mapped["atlas_mapped"] = mapped["drug_id"].notna()
    mapped["mapping_status"] = "atlas_inchikey_match"
    mapped.loc[~mapped["atlas_mapped"], "mapping_status"] = "unmapped_pubchem_cid"
    if keep_unmapped:
        mapped["drug_id"] = mapped["drug_id"].fillna("pubchem_cid:" + mapped["pubchem_cid"].astype(str))
        mapped["drug_name"] = mapped["drug_name"].fillna("PubChem CID " + mapped["pubchem_cid"].astype(str))
    else:
        mapped = mapped[mapped["atlas_mapped"]].copy()
    out = pd.DataFrame(
        {
            "drug_id": mapped["drug_id"],
            "drug_name": mapped["drug_name"],
            "inchikey": mapped["inchikey"],
            "smiles": mapped["smiles"],
            "target_id": mapped["target_id"],
            "accession": mapped["target_id"],
            "assay_id": int(aid),
            "activity_outcome": mapped["activity_outcome"],
            "activity_type": mapped["activity_type"],
            "activity_relation": mapped["activity_relation"],
            "activity_units": mapped["activity_units"],
            "source": "PubChem BioAssay",
            "pubchem_cid": mapped["pubchem_cid"],
            "atlas_mapped": mapped["atlas_mapped"],
            "mapping_status": mapped["mapping_status"],
            "activity_value_nM": mapped["activity_value_nM"],
        }
    ).drop_duplicates()
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, sep="\t", index=False)
    manifest = {
        "source": "PubChem BioAssay downloaded AID table",
        "aid": int(aid),
        "aid_csv_path": str(aid_csv_path),
        "summary_json_path": str(summary_json_path or ""),
        "output": str(path),
        "rows": int(len(out)),
        "atlas_mapped_rows": int(out["atlas_mapped"].sum()) if "atlas_mapped" in out.columns else 0,
        "target_accessions": target_accessions,
        "label_policy": "Active -> 1, inactive -> 0, inconclusive/unspecified -> -1 unless separately configured.",
    }
    path.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return out


def _atlas_cid_map(mapping_path: str | Path) -> dict[int, list[dict[str, Any]]]:
    mapping = pd.read_csv(mapping_path, low_memory=False)
    cid_col = "pubchem_cid_resolved" if "pubchem_cid_resolved" in mapping.columns else "pubchem_cid"
    rows: dict[int, list[dict[str, Any]]] = {}
    if cid_col not in mapping.columns:
        return rows
    for row in mapping.to_dict("records"):
        cid = pd.to_numeric(pd.Series([row.get(cid_col)]), errors="coerce").iloc[0]
        if pd.isna(cid):
            continue
        scheme = str(row.get("scheme") or "drug").strip()
        file_num = row.get("file_num")
        try:
            drug_id = f"{scheme}_{int(file_num):07d}"
        except Exception:
            drug_id = str(row.get("drug_id") or row.get("display_name") or cid)
        rows.setdefault(int(cid), []).append(
            {
                "drug_id": drug_id,
                "drug_name": row.get("display_name") or row.get("generic_name") or row.get("pubchem_name"),
                "inchikey": row.get("inchikey") or row.get("remark_inchikey"),
                "smiles": row.get("smiles") or row.get("remark_smiles"),
            }
        )
    return rows


def stage_pubchem_bioassay_for_targets(
    pair_table_path: str | Path,
    mapping_path: str | Path,
    out_path: str | Path,
    *,
    max_aids_per_target: int = 25,
    sleep_sec: float = 0.2,
    resolve_cid_inchikey: bool = False,
    cid_properties_cache: str | Path | None = None,
) -> pd.DataFrame:
    pairs = pd.read_csv(pair_table_path, low_memory=False)
    target_col = "target_uniprot" if "target_uniprot" in pairs.columns else "target_id"
    target_map = pairs[["target_id", target_col]].dropna().drop_duplicates()
    cid_map = _atlas_cid_map(mapping_path)
    atlas_index = load_atlas_compound_index(mapping_path)
    cid_property_cache_path = Path(cid_properties_cache) if cid_properties_cache is not None else Path(out_path).with_suffix(".cid_properties.csv")
    cid_property_cache = (
        pd.read_csv(cid_property_cache_path, low_memory=False)
        if resolve_cid_inchikey and cid_property_cache_path.exists()
        else pd.DataFrame(columns=["pubchem_cid", "inchikey", "smiles"])
    )
    rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    for target in target_map.to_dict("records"):
        accession = str(target.get(target_col) or "").strip()
        if not accession:
            continue
        aids = aids_for_uniprot(accession)
        used_aids = aids[: max(0, int(max_aids_per_target))]
        audit_rows.append(
            {
                "target_id": target.get("target_id"),
                "target_uniprot": accession,
                "n_aids_total": len(aids),
                "n_aids_used": len(used_aids),
            }
        )
        for aid in used_aids:
            active = cids_for_aid(aid, "active")
            time.sleep(max(0.0, sleep_sec))
            inactive = cids_for_aid(aid, "inactive")
            time.sleep(max(0.0, sleep_sec))
            for state, cids in [(1, active), (0, inactive)]:
                mapped_cids = set(cid_map).intersection(cids)
                cid_props = pd.DataFrame()
                if resolve_cid_inchikey:
                    unmapped_cids = sorted(set(cids) - mapped_cids)
                    cached = cid_property_cache[
                        pd.to_numeric(cid_property_cache.get("pubchem_cid"), errors="coerce").isin(unmapped_cids)
                    ].copy()
                    missing_cids = sorted(set(unmapped_cids) - set(pd.to_numeric(cached.get("pubchem_cid"), errors="coerce").dropna().astype(int)))
                    fetched = properties_for_cids(missing_cids)
                    if not fetched.empty:
                        cid_property_cache = pd.concat([cid_property_cache, fetched], ignore_index=True).drop_duplicates("pubchem_cid")
                    cid_props = pd.concat([cached, fetched], ignore_index=True).drop_duplicates("pubchem_cid")
                    if not cid_props.empty:
                        # Exact CID matching was already handled above. Fill salt/parent mismatches through InChIKey.
                        cid_props = map_by_inchikey(cid_props, atlas_index, source_inchikey_col="inchikey")
                for cid in cids:
                    for mapped in cid_map.get(int(cid), []):
                        rows.append(
                            {
                                "drug_id": mapped["drug_id"],
                                "drug_name": mapped.get("drug_name"),
                                "inchikey": mapped.get("inchikey"),
                                "smiles": mapped.get("smiles"),
                                "target_id": target.get("target_id"),
                                "accession": accession,
                                "assay_id": aid,
                                "activity_outcome": state,
                                "activity_type": "PubChem BioAssay summary outcome",
                                "activity_relation": "",
                                "activity_units": "",
                                "source": "PubChem BioAssay",
                                "pubchem_cid": cid,
                            }
                        )
                if not cid_props.empty:
                    for mapped in cid_props[cid_props["drug_id"].notna()].to_dict("records"):
                        rows.append(
                            {
                                "drug_id": mapped["drug_id"],
                                "drug_name": mapped.get("drug_name"),
                                "inchikey": mapped.get("inchikey"),
                                "smiles": mapped.get("smiles"),
                                "target_id": target.get("target_id"),
                                "accession": accession,
                                "assay_id": aid,
                                "activity_outcome": state,
                                "activity_type": "PubChem BioAssay summary outcome",
                                "activity_relation": "",
                                "activity_units": "",
                                "source": "PubChem BioAssay",
                                "pubchem_cid": mapped.get("pubchem_cid"),
                            }
                        )
    out = pd.DataFrame(rows, columns=PUBCHEM_BIOASSAY_COLUMNS).drop_duplicates()
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, sep="\t", index=False)
    audit = pd.DataFrame(audit_rows)
    audit.to_csv(path.with_suffix(".target_audit.csv"), index=False)
    if resolve_cid_inchikey:
        cid_property_cache_path.parent.mkdir(parents=True, exist_ok=True)
        cid_property_cache.drop_duplicates("pubchem_cid").to_csv(cid_property_cache_path, index=False)
    manifest = {
        "source": "PubChem BioAssay PUG-REST",
        "pair_table": str(pair_table_path),
        "mapping": str(mapping_path),
        "output": str(path),
        "rows": int(len(out)),
        "target_rows": int(len(audit)),
        "max_aids_per_target": int(max_aids_per_target),
        "resolve_cid_inchikey": bool(resolve_cid_inchikey),
        "cid_properties_cache": str(cid_property_cache_path) if resolve_cid_inchikey else "",
        "label_policy": "Active -> 1, inactive -> 0, inconclusive/absent not downloaded and remains unknown.",
    }
    path.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return out

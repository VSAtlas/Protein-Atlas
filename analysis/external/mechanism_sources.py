from __future__ import annotations

import gzip
import json
import re
import time
import urllib.error
import urllib.request
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

try:
    from src.path_router.path_router import _ensure_router_roots
except Exception:  # pragma: no cover - path_router is not available in minimal package installs.
    _ensure_router_roots = None  # type: ignore[assignment]


SIDER_SE_URL = "http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz"
SIDER_NAMES_URL = "http://sideeffects.embl.de/media/download/drug_names.tsv"
REACTOME_UNIPROT_URL = "https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt"
CTD_CHEM_GENE_URLS = (
    "https://ctdbase.org/reports/CTD_chem_gene_ixns.tsv.gz",
    "https://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz",
)
CTD_CHEM_DISEASE_URLS = (
    "https://ctdbase.org/reports/CTD_chemicals_diseases.tsv.gz",
    "https://ctdbase.org/reports/CTD_chemicals_diseases.csv.gz",
)
CTD_GENE_DISEASE_URLS = (
    "https://ctdbase.org/reports/CTD_genes_diseases.tsv.gz",
    "https://ctdbase.org/reports/CTD_genes_diseases.csv.gz",
)
OPEN_TARGETS_GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"

NORMALIZED_RELAXED = re.compile(r"[^a-z0-9]+")


@dataclass(frozen=True)
class MechanismSourcePaths:
    sider: Path
    ctd: Path
    opentargets_safety: Path
    reactome: Path
    manifest: Path


def repo_root() -> Path:
    if _ensure_router_roots is not None:
        try:
            return Path(_ensure_router_roots().overall)
        except Exception:
            pass
    return Path(__file__).resolve().parents[3]


def mechanism_source_paths(base_dir: str | Path | None = None) -> MechanismSourcePaths:
    base = Path(base_dir) if base_dir is not None else repo_root() / "data" / "external"
    return MechanismSourcePaths(
        sider=base / "sider" / "drug_adr.tsv",
        ctd=base / "ctd" / "ctd_edges.tsv",
        opentargets_safety=base / "opentargets" / "safety.tsv",
        reactome=base / "reactome" / "pathways.tsv",
        manifest=base / "mechanism_sources_manifest.json",
    )


def download_file(urls: str | Sequence[str], out_path: str | Path, *, overwrite: bool = False, retries: int = 3) -> Path:
    out = Path(out_path)
    if out.exists() and not overwrite:
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    candidates = [urls] if isinstance(urls, str) else list(urls)
    last_error: Exception | None = None
    for url in candidates:
        request = urllib.request.Request(str(url), headers={"User-Agent": "AtlasMechanismSources/1.0"})
        for attempt in range(max(1, retries)):
            try:
                with urllib.request.urlopen(request, timeout=180) as response:
                    out.write_bytes(response.read())
                return out
            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
                last_error = exc
                if attempt + 1 < retries:
                    time.sleep(2.0 * float(attempt + 1))
                continue
    raise RuntimeError(f"failed to download {candidates} -> {out}: {last_error}")


def _norm_key(value: Any) -> str:
    return NORMALIZED_RELAXED.sub("", str(value or "").lower())


def _nonempty(value: Any) -> bool:
    return pd.notna(value) and bool(str(value).strip())


def _split_aliases(value: Any) -> Iterable[str]:
    if not _nonempty(value):
        return []
    return [part.strip() for part in re.split(r"[;|]", str(value)) if part.strip()]


def _read_table(path: str | Path, *, names: list[str] | None = None, comment: str | None = "#") -> pd.DataFrame:
    p = Path(path)
    sep = "\t" if ".tsv" in "".join(p.suffixes).lower() or p.suffix.lower() in {".tab", ".txt"} else ","
    return pd.read_csv(p, sep=sep, comment=comment, names=names, compression="infer", low_memory=False)


def _read_pair_table_ids(pair_table: str | Path | None) -> tuple[set[str], set[str]]:
    if pair_table is None or not Path(pair_table).exists():
        return set(), set()
    df = pd.read_csv(pair_table, usecols=lambda col: col in {"drug_id", "ligand_base", "target_id", "target_gene", "target_uniprot", "pdb_id"})
    drugs = {str(v).strip() for v in df.get("drug_id", pd.Series(dtype=str)).dropna().unique() if str(v).strip()}
    if "ligand_base" in df.columns:
        drugs.update(str(v).strip() for v in df["ligand_base"].dropna().unique() if str(v).strip())
    targets: set[str] = set()
    for column in ("target_id", "target_gene", "target_uniprot", "pdb_id"):
        if column in df.columns:
            targets.update(str(v).strip() for v in df[column].dropna().unique() if str(v).strip())
    return drugs, targets


def build_fda_alias_map(mapping_path: str | Path | None) -> dict[str, set[str]]:
    if mapping_path is None or not Path(mapping_path).exists():
        return {}
    mapping = pd.read_csv(mapping_path, low_memory=False)
    out: dict[str, set[str]] = {}
    alias_columns = [
        "display_name",
        "generic_name",
        "pubchem_name",
        "pubchem_record_title",
        "rxnorm_generic_name",
        "drugcentral_generic_name",
        "remark_name",
        "sdf_title",
        "brand_names",
        "rxnorm_brand_names",
        "drugcentral_brand_names",
        "pubchem_synonyms",
        "cas",
        "inchikey",
        "remark_inchikey",
    ]
    for _idx, row in mapping.iterrows():
        scheme = str(row.get("scheme") or "").strip()
        file_num = row.get("file_num")
        if not scheme or pd.isna(file_num):
            continue
        try:
            drug_id = f"{scheme}_{int(file_num):07d}"
        except Exception:
            drug_id = f"{scheme}_{str(file_num).strip()}"
        aliases: set[str] = {drug_id, drug_id.replace("_", " ")}
        for column in alias_columns:
            value = row.get(column)
            aliases.update(str(part) for part in _split_aliases(value))
        for alias in aliases:
            key = _norm_key(alias)
            if key:
                out.setdefault(key, set()).add(drug_id)
    return out


def _map_drug_alias(value: Any, alias_map: Mapping[str, set[str]]) -> set[str]:
    return set(alias_map.get(_norm_key(value), set()))


def normalize_sider(
    side_effects_path: str | Path,
    names_path: str | Path | None,
    out_path: str | Path,
    *,
    fda_mapping_path: str | Path | None = None,
) -> pd.DataFrame:
    columns = ["stitch_flat", "stitch_stereo", "umls_label_id", "meddra_type", "umls_id", "adr"]
    se = _read_table(side_effects_path, names=columns)
    names = pd.DataFrame(columns=["stitch_stereo", "drug_name"])
    if names_path is not None and Path(names_path).exists():
        names = _read_table(names_path, names=["stitch_stereo", "drug_name"])
    alias_map = build_fda_alias_map(fda_mapping_path)
    names = names.drop_duplicates("stitch_stereo")
    merged = se.merge(names.rename(columns={"drug_name": "drug_name_stereo"}), on="stitch_stereo", how="left")
    merged = merged.merge(
        names.rename(columns={"stitch_stereo": "stitch_flat", "drug_name": "drug_name_flat"}),
        on="stitch_flat",
        how="left",
    )
    merged["drug_name"] = merged["drug_name_stereo"].fillna(merged["drug_name_flat"])
    rows: list[dict[str, Any]] = []
    for row in merged.itertuples(index=False):
        mapped = _map_drug_alias(getattr(row, "drug_name", ""), alias_map)
        drug_ids = sorted(mapped) if mapped else [getattr(row, "stitch_stereo")]
        for drug_id in drug_ids:
            rows.append(
                {
                    "drug_id": drug_id,
                    "adr": getattr(row, "adr"),
                    "source": "SIDER",
                    "_source_stitch_id": getattr(row, "stitch_stereo"),
                    "_source_drug_name": getattr(row, "drug_name", pd.NA),
                    "_source_meddra_type": getattr(row, "meddra_type"),
                    "_source_umls_id": getattr(row, "umls_id"),
                    "_mapping_status": "mapped_to_atlas_fda" if mapped else "unmapped_stitch_id",
                }
            )
    out = pd.DataFrame(rows).drop_duplicates()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    return out


def _ctd_column(df: pd.DataFrame, names: Sequence[str]) -> pd.Series:
    lower = {str(col).lower(): col for col in df.columns}
    for name in names:
        col = lower.get(name.lower())
        if col is not None:
            return df[col]
    return pd.Series(pd.NA, index=df.index)


def _ctd_read(path: str | Path) -> pd.DataFrame:
    return _read_table(path)


def _ctd_field_names(path: str | Path) -> list[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:  # type: ignore[arg-type]
        for line in handle:
            if line.startswith("# Fields:"):
                fields = next(handle, "").strip()
                return [item.strip() for item in fields.lstrip("#").split("\t") if item.strip()]
    return []


def _ctd_chunks(path: str | Path, *, chunksize: int = 250_000) -> Iterable[pd.DataFrame]:
    fields = _ctd_field_names(path)
    if not fields:
        yield _ctd_read(path)
        return
    yield from pd.read_csv(
        path,
        sep="\t",
        comment="#",
        names=fields,
        compression="infer",
        low_memory=False,
        chunksize=chunksize,
    )


def _mapped_chemical_ids(row: Mapping[str, Any], alias_map: Mapping[str, set[str]]) -> set[str]:
    hits: set[str] = set()
    for key in ("chemical_id", "_source_chemical_id", "_source_cas"):
        hits.update(_map_drug_alias(row.get(key), alias_map))
    return hits


def _alias_keys_for_drug_ids(alias_map: Mapping[str, set[str]], drug_ids: set[str]) -> set[str]:
    if not drug_ids:
        return set(alias_map)
    return {key for key, mapped_drugs in alias_map.items() if mapped_drugs & drug_ids}


def _chemical_filter_mask(df: pd.DataFrame, allowed_alias_keys: set[str]) -> pd.Series:
    if not allowed_alias_keys:
        return pd.Series(True, index=df.index)
    mask = pd.Series(False, index=df.index)
    for names in (["ChemicalName", "chemical_name", "chemical_id"], ["ChemicalID", "chemical_id"], ["CasRN", "cas"]):
        values = _ctd_column(df, names)
        mask = mask | values.map(lambda value: _norm_key(value) in allowed_alias_keys if _nonempty(value) else False)
    return mask


def _append_ctd_rows(
    rows: list[dict[str, Any]],
    frame: pd.DataFrame,
    *,
    alias_map: Mapping[str, set[str]],
    drug_ids: set[str],
) -> None:
    for row in frame.to_dict(orient="records"):
        mapped = _mapped_chemical_ids(row, alias_map)
        if mapped:
            keep = sorted(mapped & drug_ids) if drug_ids else sorted(mapped)
            for drug_id in keep:
                new_row = dict(row)
                new_row["chemical_id"] = drug_id
                new_row["_mapping_status"] = "mapped_to_atlas_fda"
                rows.append(new_row)
        elif not drug_ids:
            row["_mapping_status"] = "unmapped_ctd_chemical"
            rows.append(row)


def normalize_ctd(
    chem_gene_path: str | Path | None,
    chem_disease_path: str | Path | None,
    gene_disease_path: str | Path | None,
    out_path: str | Path,
    *,
    fda_mapping_path: str | Path | None = None,
    target_filter: set[str] | None = None,
    drug_filter: set[str] | None = None,
) -> pd.DataFrame:
    alias_map = build_fda_alias_map(fda_mapping_path)
    allowed_alias_keys = _alias_keys_for_drug_ids(alias_map, set(drug_filter or set()))
    target_keys = {_norm_key(v) for v in (target_filter or set()) if _norm_key(v)}
    drug_ids = set(drug_filter or set())
    rows: list[dict[str, Any]] = []
    if chem_gene_path and Path(chem_gene_path).exists():
        for df in _ctd_chunks(chem_gene_path):
            gene = _ctd_column(df, ["GeneSymbol", "gene_symbol", "GeneID", "target_id"])
            if target_keys:
                df = df[gene.map(lambda value: _norm_key(value) in target_keys)]
                if df.empty:
                    continue
            df = df[_chemical_filter_mask(df, allowed_alias_keys)]
            if df.empty:
                continue
            out = pd.DataFrame(
                {
                    "chemical_id": _ctd_column(df, ["ChemicalName", "chemical_name", "chemical_id"]),
                    "gene_id": _ctd_column(df, ["GeneSymbol", "gene_symbol", "GeneID", "target_id"]),
                    "disease_id": pd.NA,
                    "edge_type": "chemical_gene_interaction",
                    "pubmed_ids": _ctd_column(df, ["PubMedIDs", "pubmed_ids"]),
                    "source": "CTD",
                    "_source_chemical_id": _ctd_column(df, ["ChemicalID", "chemical_id"]),
                    "_source_cas": _ctd_column(df, ["CasRN", "cas"]),
                    "_source_interaction": _ctd_column(df, ["Interaction", "interaction"]),
                }
            )
            _append_ctd_rows(rows, out, alias_map=alias_map, drug_ids=drug_ids)
    if chem_disease_path and Path(chem_disease_path).exists():
        for df in _ctd_chunks(chem_disease_path):
            df = df[_chemical_filter_mask(df, allowed_alias_keys)]
            if df.empty:
                continue
            gene = _ctd_column(df, ["InferenceGeneSymbol", "gene_symbol", "target_id"])
            if target_keys:
                target_hit = gene.map(lambda value: _norm_key(value) in target_keys if _nonempty(value) else False)
                df = df[target_hit]
                gene = gene[target_hit]
                if df.empty:
                    continue
            out = pd.DataFrame(
                {
                    "chemical_id": _ctd_column(df, ["ChemicalName", "chemical_name", "chemical_id"]),
                    "gene_id": gene,
                    "disease_id": _ctd_column(df, ["DiseaseName", "disease_name", "disease_id"]),
                    "edge_type": "chemical_disease_association",
                    "pubmed_ids": _ctd_column(df, ["PubMedIDs", "pubmed_ids"]),
                    "source": "CTD",
                    "_source_chemical_id": _ctd_column(df, ["ChemicalID", "chemical_id"]),
                    "_source_disease_id": _ctd_column(df, ["DiseaseID", "disease_id"]),
                    "_source_cas": _ctd_column(df, ["CasRN", "cas"]),
                    "_source_direct_evidence": _ctd_column(df, ["DirectEvidence", "direct_evidence"]),
                }
            )
            _append_ctd_rows(rows, out, alias_map=alias_map, drug_ids=drug_ids)
    if gene_disease_path and Path(gene_disease_path).exists():
        for df in _ctd_chunks(gene_disease_path):
            gene = _ctd_column(df, ["GeneSymbol", "gene_symbol", "GeneID", "target_id"])
            if target_keys:
                df = df[gene.map(lambda value: _norm_key(value) in target_keys)]
                if df.empty:
                    continue
            out = pd.DataFrame(
                {
                    "chemical_id": pd.NA,
                    "gene_id": _ctd_column(df, ["GeneSymbol", "gene_symbol", "GeneID", "target_id"]),
                    "disease_id": _ctd_column(df, ["DiseaseName", "disease_name", "disease_id"]),
                    "edge_type": "gene_disease_association",
                    "pubmed_ids": _ctd_column(df, ["PubMedIDs", "pubmed_ids"]),
                    "source": "CTD",
                    "_source_gene_id": _ctd_column(df, ["GeneID", "gene_id"]),
                    "_source_disease_id": _ctd_column(df, ["DiseaseID", "disease_id"]),
                    "_source_direct_evidence": _ctd_column(df, ["DirectEvidence", "direct_evidence"]),
                    "_mapping_status": "target_filtered_gene_disease",
                }
            )
            rows.extend(out.to_dict(orient="records"))
    combined = pd.DataFrame(rows)
    if combined.empty:
        combined = pd.DataFrame(columns=["chemical_id", "gene_id", "disease_id", "edge_type", "pubmed_ids", "source"])
    combined = combined.drop_duplicates()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(out_path, sep="\t", index=False)
    return combined


def normalize_reactome(
    reactome_path: str | Path,
    out_path: str | Path,
    *,
    target_filter: set[str] | None = None,
) -> pd.DataFrame:
    columns = ["target_id", "pathway_id", "pathway_url", "pathway", "evidence_code", "species"]
    df = _read_table(reactome_path, names=columns, comment=None)
    human = df[df["species"].astype(str).str.contains("Homo sapiens", case=False, na=False)].copy()
    if target_filter:
        target_keys = {_norm_key(v) for v in target_filter if _norm_key(v)}
        human = human[human["target_id"].map(lambda value: _norm_key(value) in target_keys)]
    out = human.assign(adr=pd.NA, source="Reactome")[
        ["target_id", "pathway", "pathway_id", "adr", "source", "evidence_code", "pathway_url", "species"]
    ].drop_duplicates()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    return out


def _graphql(query: str, variables: Mapping[str, Any], *, timeout: int = 30) -> dict[str, Any]:
    payload = json.dumps({"query": query, "variables": dict(variables)}).encode("utf-8")
    request = urllib.request.Request(
        OPEN_TARGETS_GRAPHQL_URL,
        data=payload,
        headers={"Content-Type": "application/json", "User-Agent": "AtlasMechanismSources/1.0"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        data = json.loads(response.read().decode("utf-8"))
    if data.get("errors"):
        raise RuntimeError(str(data["errors"]))
    return dict(data.get("data") or {})


def _search_opentargets_target(term: str) -> dict[str, Any] | None:
    query = (
        'query Search($q:String!){ search(queryString:$q, entityNames:["target"], page:{index:0,size:3})'
        "{ hits { id name entity } } }"
    )
    data = _graphql(query, {"q": term})
    hits = ((data.get("search") or {}).get("hits") or [])
    return dict(hits[0]) if hits else None


def _fetch_opentargets_safety(ensembl_id: str) -> dict[str, Any]:
    query = (
        "query Target($id:String!){ target(ensemblId:$id){ approvedSymbol approvedName proteinIds{ id source } "
        "safetyLiabilities{ event eventId datasource effects{ direction dosing } "
        "biosamples{ tissueLabel tissueId cellLabel cellFormat } studies{ name description type } } } }"
    )
    data = _graphql(query, {"id": ensembl_id})
    return dict(data.get("target") or {})


def normalize_opentargets_safety(
    target_terms: Iterable[str],
    out_path: str | Path,
    *,
    sleep_sec: float = 0.2,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    seen_terms: set[str] = set()
    for raw_term in target_terms:
        term = str(raw_term).strip()
        key = _norm_key(term)
        if not key or key in seen_terms:
            continue
        seen_terms.add(key)
        try:
            hit = _search_opentargets_target(term)
            if not hit:
                rows.append({"target_id": term, "adr": pd.NA, "confidence": pd.NA, "pubmed_ids": "", "source": "OpenTargetsSafety", "_mapping_status": "unmapped_target"})
                continue
            target = _fetch_opentargets_safety(str(hit["id"]))
        except Exception as exc:
            rows.append({"target_id": term, "adr": pd.NA, "confidence": pd.NA, "pubmed_ids": "", "source": "OpenTargetsSafety", "_mapping_status": f"query_failed:{type(exc).__name__}"})
            continue
        time.sleep(max(0.0, sleep_sec))
        symbol = target.get("approvedSymbol") or term
        uniprots = [item.get("id") for item in target.get("proteinIds") or [] if str(item.get("source") or "").lower() == "uniprot_swissprot"]
        liabilities = target.get("safetyLiabilities") or []
        target_ids = [symbol] + [str(item) for item in uniprots if item]
        if not liabilities:
            for target_id in target_ids:
                rows.append({"target_id": target_id, "adr": pd.NA, "confidence": 0.0, "pubmed_ids": "", "source": "OpenTargetsSafety", "_source_ensembl_id": hit.get("id"), "_source_approved_symbol": symbol, "_mapping_status": "mapped_no_safety_liability"})
            continue
        for liability in liabilities:
            event = liability.get("event") or liability.get("eventId")
            studies = liability.get("studies") or []
            pubmed_ids = ";".join(str(study.get("name")) for study in studies if study.get("name"))
            metadata = json.dumps(
                {
                    "eventId": liability.get("eventId"),
                    "datasource": liability.get("datasource"),
                    "effects": liability.get("effects"),
                    "biosamples": liability.get("biosamples"),
                },
                sort_keys=True,
            )
            for target_id in target_ids:
                rows.append(
                    {
                        "target_id": target_id,
                        "adr": event,
                        "confidence": 1.0,
                        "pubmed_ids": pubmed_ids,
                        "source": "OpenTargetsSafety",
                        "_source_ensembl_id": hit.get("id"),
                        "_source_approved_symbol": symbol,
                        "_source_approved_name": target.get("approvedName"),
                        "_metadata_json": metadata,
                        "_mapping_status": "mapped_safety_liability",
                    }
                )
    out = pd.DataFrame(rows)
    if out.empty:
        out = pd.DataFrame(columns=["target_id", "adr", "confidence", "pubmed_ids", "source"])
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    out.drop_duplicates().to_csv(out_path, sep="\t", index=False)
    return out


def stage_mechanism_sources(
    *,
    base_dir: str | Path | None = None,
    pair_table: str | Path | None = None,
    fda_mapping: str | Path | None = None,
    download: bool = True,
    overwrite: bool = False,
    sources: Sequence[str] = ("sider", "ctd", "opentargets_safety", "reactome"),
) -> dict[str, Any]:
    base = Path(base_dir) if base_dir is not None else repo_root() / "data" / "external"
    paths = mechanism_source_paths(base)
    raw = base / "raw_mechanism_sources"
    drugs, targets = _read_pair_table_ids(pair_table)
    previous_manifest: dict[str, Any] = {}
    if paths.manifest.exists():
        try:
            previous_manifest = json.loads(paths.manifest.read_text(encoding="utf-8"))
        except Exception:
            previous_manifest = {}
    manifest: dict[str, Any] = {
        "source_policy": "Downloaded and normalized into the current repository. Missing evidence remains unknown, not negative.",
        "base_dir": str(base),
        "pair_table": str(pair_table) if pair_table else None,
        "fda_mapping": str(fda_mapping) if fda_mapping else None,
        "outputs": dict(previous_manifest.get("outputs") or {}),
        "row_counts": dict(previous_manifest.get("row_counts") or {}),
        "sources": sorted(set(previous_manifest.get("sources") or []) | set(sources)),
    }
    for source_name, output_path in {
        "sider": paths.sider,
        "ctd": paths.ctd,
        "opentargets_safety": paths.opentargets_safety,
        "reactome": paths.reactome,
    }.items():
        if output_path.exists():
            manifest["outputs"].setdefault(source_name, str(output_path))
            manifest["sources"] = sorted(set(manifest["sources"]) | {source_name})
            if source_name not in manifest["row_counts"]:
                try:
                    manifest["row_counts"][source_name] = int(pd.read_csv(output_path, sep="\t", usecols=[0]).shape[0])
                except Exception:
                    manifest["row_counts"][source_name] = None
    if "sider" in sources:
        se = raw / "sider" / "meddra_all_se.tsv.gz"
        names = raw / "sider" / "drug_names.tsv"
        if download:
            download_file(SIDER_SE_URL, se, overwrite=overwrite)
            download_file(SIDER_NAMES_URL, names, overwrite=overwrite)
        df = normalize_sider(se, names, paths.sider, fda_mapping_path=fda_mapping)
        manifest["outputs"]["sider"] = str(paths.sider)
        manifest["row_counts"]["sider"] = int(len(df))
    if "ctd" in sources:
        cg = raw / "ctd" / "CTD_chem_gene_ixns.tsv.gz"
        cd = raw / "ctd" / "CTD_chemicals_diseases.tsv.gz"
        gd = raw / "ctd" / "CTD_genes_diseases.tsv.gz"
        if download:
            download_file(CTD_CHEM_GENE_URLS, cg, overwrite=overwrite)
            download_file(CTD_CHEM_DISEASE_URLS, cd, overwrite=overwrite)
            download_file(CTD_GENE_DISEASE_URLS, gd, overwrite=overwrite)
        df = normalize_ctd(cg, cd, gd, paths.ctd, fda_mapping_path=fda_mapping, target_filter=targets, drug_filter=drugs)
        manifest["outputs"]["ctd"] = str(paths.ctd)
        manifest["row_counts"]["ctd"] = int(len(df))
    if "reactome" in sources:
        reactome = raw / "reactome" / "UniProt2Reactome_All_Levels.txt"
        if download:
            download_file(REACTOME_UNIPROT_URL, reactome, overwrite=overwrite)
        df = normalize_reactome(reactome, paths.reactome, target_filter=targets)
        manifest["outputs"]["reactome"] = str(paths.reactome)
        manifest["row_counts"]["reactome"] = int(len(df))
    if "opentargets_safety" in sources:
        df = normalize_opentargets_safety(targets, paths.opentargets_safety)
        manifest["outputs"]["opentargets_safety"] = str(paths.opentargets_safety)
        manifest["row_counts"]["opentargets_safety"] = int(len(df))
    paths.manifest.parent.mkdir(parents=True, exist_ok=True)
    paths.manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest

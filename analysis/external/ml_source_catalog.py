from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.external.source_tables import download_to_cache


SOURCE_CATALOG: dict[str, dict[str, Any]] = {
    "chembl": {
        "role": "curated binding and functional activity labels",
        "default_path": "data/external/chembl/bioactivity.tsv",
        "filtered_path": "data/external/chembl/bioactivity.tsv",
        "download_url": "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/",
        "size_class": "large",
        "production_use": "core direct assay evidence after assay-type and target-confidence filtering",
        "auto_downloadable": False,
    },
    "bindingdb": {
        "role": "direct binding labels",
        "default_path": "data/external/bindingdb/BindingDB_All_202605_tsv.zip",
        "filtered_path": "data/external/bindingdb/bioactivity.tsv",
        "download_url": "https://www.bindingdb.org/rwd/bind/downloads/BindingDB_All_202605_tsv.zip",
        "size_class": "large",
        "production_use": "core direct assay evidence after filtering to Atlas drugs/targets",
        "auto_downloadable": True,
    },
    "toxcast": {
        "role": "toxicology-relevant HTS activity labels",
        "default_path": "data/external/toxcast/invitrodb_v4_3.zip",
        "filtered_path": "data/external/toxcast/bioactivity.tsv",
        "download_url": "https://www.epa.gov/comptox-tools/exploring-toxcast-data",
        "size_class": "large",
        "production_use": "independent HTS activity benchmark after assay QC, cytotoxicity, and target-mapping filters",
        "auto_downloadable": False,
    },
    "tox21": {
        "role": "Tox21 toxicology HTS activity labels",
        "default_path": "data/external/tox21/bioactivity.tsv",
        "filtered_path": "data/external/tox21/bioactivity.tsv",
        "download_url": "https://www.epa.gov/comptox-tools/exploring-toxcast-data",
        "size_class": "query_or_filtered",
        "production_use": "independent HTS activity benchmark after assay QC and target-mapping filters",
        "auto_downloadable": False,
    },
    "pubchem_bioassay": {
        "role": "explicit active/inactive assay calls",
        "default_path": "data/external/pubchem_bioassay/pubchem_bioassay.tsv",
        "filtered_path": "data/external/pubchem_bioassay/pubchem_bioassay.tsv",
        "download_url": "https://pubchem.ncbi.nlm.nih.gov/docs/bioassays",
        "size_class": "query_or_filtered",
        "production_use": "source-specific measured evidence; only assay-defined active/inactive outcomes",
        "auto_downloadable": False,
    },
    "iuphar_gtopdb": {
        "role": "expert-curated ligand-target pharmacology positives",
        "default_path": "data/external/iuphar_gtopdb/interactions.tsv",
        "filtered_path": "data/external/iuphar_gtopdb/bioactivity.tsv",
        "download_url": "https://www.guidetopharmacology.org/DATA/interactions.tsv",
        "size_class": "small",
        "production_use": "high-precision positive ligand-target evidence; absence and weak affinity are not negatives",
        "auto_downloadable": True,
    },
    "ncats_inxight": {
        "role": "drug status, pharmacokinetic, and safety metadata",
        "default_path": "data/external/ncats_inxight/frdb-v2024-12-30.zip",
        "filtered_path": "data/external/ncats_inxight/pk_table.tsv",
        "download_url": "https://drugs.ncats.io/downloads-public",
        "size_class": "small",
        "production_use": "exposure/PK and safety context only; not a direct drug-target activity label source",
        "auto_downloadable": False,
    },
    "excape": {
        "role": "large chemogenomics training-only evidence",
        "default_path": "data/external/excape/excape.tsv",
        "filtered_path": "data/external/excape/bioactivity.tsv",
        "download_url": "https://jcheminf.biomedcentral.com/articles/10.1186/s13321-017-0203-5",
        "size_class": "very_large",
        "production_use": "training-only; high overlap/leakage risk with ChEMBL/PubChem evaluations",
        "auto_downloadable": False,
    },
    "lit_pcba_subset": {
        "role": "benchmark/sensitivity active-inactive panel",
        "default_path": "data/external/lit_pcba/LITPCBA_9t_subset.tar.xz",
        "filtered_path": "data/external/lit_pcba/benchmark_labels.tsv",
        "download_url": "https://zenodo.org/records/10846630/files/LITPCBA_9t_subset.tar.xz?download=1",
        "size_class": "small",
        "production_use": "benchmark-only; do not merge into production truth labels",
        "auto_downloadable": True,
    },
    "dude": {
        "role": "docking benchmark actives/decoys",
        "default_path": "data/external/dude/dude.tsv",
        "filtered_path": "data/external/dude/benchmark_labels.tsv",
        "download_url": "http://dude.docking.org/targets",
        "size_class": "benchmark",
        "production_use": "benchmark-only decoys; never production negatives",
        "auto_downloadable": False,
    },
    "muv": {
        "role": "virtual-screening benchmark",
        "default_path": "data/external/muv/muv.tsv",
        "filtered_path": "data/external/muv/benchmark_labels.tsv",
        "download_url": "https://doi.org/10.1021/ci8002649",
        "size_class": "benchmark",
        "production_use": "benchmark-only; do not merge into core Atlas truth",
        "auto_downloadable": False,
    },
}


def stage_ml_source_catalog(
    out_dir: str | Path,
    *,
    download_small: bool = False,
    download_large: bool = False,
    sources: list[str] | None = None,
    overwrite: bool = False,
) -> pd.DataFrame:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    requested = {str(source).strip().lower() for source in (sources or []) if str(source).strip()}
    for name, meta in SOURCE_CATALOG.items():
        if requested and name not in requested:
            continue
        local_path = Path(str(meta["default_path"]))
        filtered_path = Path(str(meta["filtered_path"]))
        downloaded = False
        error = ""
        should_download = bool(meta.get("auto_downloadable")) and (
            (download_small and meta["size_class"] == "small")
            or (download_large and meta["size_class"] in {"large", "very_large"})
        )
        if should_download:
            try:
                download_to_cache(str(meta["download_url"]), local_path, overwrite=overwrite)
                downloaded = True
            except Exception as exc:
                error = str(exc)
        rows.append(
            {
                "source": name,
                "role": meta["role"],
                "local_path": str(local_path),
                "exists": local_path.exists(),
                "size_bytes": local_path.stat().st_size if local_path.exists() else pd.NA,
                "filtered_path": str(filtered_path),
                "filtered_exists": filtered_path.exists(),
                "download_url": meta["download_url"],
                "size_class": meta["size_class"],
                "auto_downloadable": bool(meta.get("auto_downloadable")),
                "downloaded_this_run": downloaded,
                "production_use": meta["production_use"],
                "staging_error": error,
            }
        )
    manifest = pd.DataFrame(rows)
    manifest_path = out / "ml_source_catalog_manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    payload = {
        "manifest": str(manifest_path),
        "download_small": download_small,
        "download_large": download_large,
        "sources": sorted(requested) if requested else "all",
        "sources_present": int(manifest["exists"].sum()),
        "filtered_sources_present": int(manifest["filtered_exists"].sum()),
        "sources_total": int(len(manifest)),
        "large_source_policy": (
            "Large sources are downloadable when --download-large is set and the catalog has a direct file URL. "
            "For query/portal sources, stage filtered local exports under filtered_path, then pass them to "
            "build_four_state_evidence."
        ),
    }
    (out / "ml_source_catalog_manifest.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return manifest

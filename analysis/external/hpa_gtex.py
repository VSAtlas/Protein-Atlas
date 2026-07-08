from __future__ import annotations

import gzip
import json
import re
import zipfile
from pathlib import Path
from typing import Any

import pandas as pd
import requests

from analysis.external.mapping import normalize_columns

HPA_BASE = "https://www.proteinatlas.org"
HPA_CONSENSUS_URL = f"{HPA_BASE}/download/tsv/rna_tissue_consensus.tsv.zip"
HPA_GTEX_URL = f"{HPA_BASE}/download/tsv/rna_tissue_gtex.tsv.zip"
HPA_ATLAS_URL = f"{HPA_BASE}/download/proteinatlas.tsv.zip"
BGEE_EXPR_URL = "https://www.bgee.org/ftp/current/download/calls/expr_calls/Homo_sapiens_expr_simple.tsv.gz"
OPENTARGETS_BASELINE_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/output/baseline_expression/"
)

DEFAULT_SOURCE_URLS = {
    "hpa_rna_tissue_consensus.tsv.zip": HPA_CONSENSUS_URL,
    "hpa_rna_tissue_gtex.tsv.zip": HPA_GTEX_URL,
    "hpa_proteinatlas.tsv.zip": HPA_ATLAS_URL,
    "bgee_homo_sapiens_expr_simple.tsv.gz": BGEE_EXPR_URL,
}


def load_tissue_expression(path: str | Path) -> pd.DataFrame:
    sep = "\t" if Path(path).suffix.lower() in {".tsv", ".tab"} else ","
    return normalize_columns(
        pd.read_csv(path, sep=sep),
        {
            "target_id": ["target_id", "uniprot", "gene"],
            "tissue": ["tissue"],
            "tissue_expression": ["tpm", "expression", "tissue_expression"],
            "source": ["source"],
        },
    )


def _clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _upper(value: Any) -> str:
    return _clean(value).upper()


def _lower(value: Any) -> str:
    return _clean(value).lower()


def _read_csv(path: Path, **kwargs: Any) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False, **kwargs)


def _read_zipped_tsv(path: Path, member: str | None = None, **kwargs: Any) -> pd.DataFrame:
    with zipfile.ZipFile(path) as archive:
        selected = member or archive.namelist()[0]
        with archive.open(selected) as handle:
            return pd.read_csv(handle, sep="\t", low_memory=False, **kwargs)


def _download_file(url: str, path: Path, *, timeout: int = 120) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0:
        return {"path": str(path), "url": url, "status": "exists", "bytes": int(path.stat().st_size)}
    with requests.get(url, stream=True, timeout=timeout) as response:
        response.raise_for_status()
        with path.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
    return {"path": str(path), "url": url, "status": "downloaded", "bytes": int(path.stat().st_size)}


def download_expression_sources(raw_dir: str | Path) -> dict[str, Any]:
    raw = Path(raw_dir)
    downloads = {}
    for filename, url in DEFAULT_SOURCE_URLS.items():
        downloads[filename] = _download_file(url, raw / filename)
    return {"raw_dir": str(raw), "downloads": downloads}


def discover_opentargets_baseline_parts() -> list[str]:
    response = requests.get(OPENTARGETS_BASELINE_URL, timeout=30)
    response.raise_for_status()
    return sorted(set(re.findall(r'href="([^"]+\.parquet)"', response.text)))


def download_opentargets_baseline(raw_dir: str | Path, *, limit: int = 0) -> dict[str, Any]:
    out = Path(raw_dir) / "opentargets_baseline_expression"
    out.mkdir(parents=True, exist_ok=True)
    parts = discover_opentargets_baseline_parts()
    if limit > 0:
        parts = parts[:limit]
    downloads = []
    for part in parts:
        downloads.append(_download_file(OPENTARGETS_BASELINE_URL + part, out / part, timeout=240))
    return {"raw_dir": str(out), "n_parts": len(parts), "downloads": downloads}


def load_target_gene_set(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    sep = "\t" if p.suffix.lower() in {".tsv", ".tab"} else ","
    df = _read_csv(p, sep=sep)
    gene_col = next(
        (
            col
            for col in (
                "gene",
                "target_gene",
                "Gene",
                "EntrezGeneSymbol",
                "HumanEntrezGeneSymbol(representative)",
            )
            if col in df.columns
        ),
        None,
    )
    if gene_col is None:
        raise ValueError(f"target table has no recognized gene column: {path}")
    out = pd.DataFrame({"target_gene": df[gene_col].map(_upper)})
    if "pdb_id" in df.columns:
        out["pdb_id"] = df["pdb_id"].map(_upper)
    if "target_uniprot" in df.columns:
        out["target_uniprot"] = df["target_uniprot"].map(_upper)
    return out.loc[out["target_gene"].astype(str).str.len().gt(0)].drop_duplicates("target_gene")


def _filter_targets(df: pd.DataFrame, target_genes: set[str], ensembl_ids: set[str] | None = None) -> pd.DataFrame:
    if "Gene name" in df.columns:
        mask = df["Gene name"].map(_upper).isin(target_genes)
    elif "Gene" in df.columns:
        mask = df["Gene"].map(_upper).isin(target_genes)
    else:
        mask = pd.Series(False, index=df.index)
    if ensembl_ids and "Gene" in df.columns:
        mask = mask | df["Gene"].map(_upper).isin(ensembl_ids)
    return df.loc[mask].copy()


def _load_hpa_metadata(raw_dir: Path, target_genes: set[str]) -> pd.DataFrame:
    path = raw_dir / "hpa_proteinatlas.tsv.zip"
    if not path.exists():
        return pd.DataFrame()
    cols = [
        "Gene",
        "Ensembl",
        "Uniprot",
        "RNA tissue specificity",
        "RNA tissue distribution",
        "RNA tissue specificity score",
        "Protein tissue specificity",
        "Protein tissue distribution",
        "Protein tissue specificity score",
    ]
    df = _read_zipped_tsv(path, usecols=lambda col: col in cols)
    df = df.loc[df["Gene"].map(_upper).isin(target_genes)].copy()
    if df.empty:
        return df
    df["target_gene"] = df["Gene"].map(_upper)
    df["gene_symbol"] = df["target_gene"]
    df["ensembl_id"] = df["Ensembl"].map(_upper)
    df["target_uniprot"] = df.get("Uniprot", pd.Series("", index=df.index)).map(_upper)
    df["hpa_specificity_category"] = df.get("RNA tissue specificity")
    df["hpa_distribution_category"] = df.get("RNA tissue distribution")
    df["hpa_specificity_score"] = pd.to_numeric(df.get("RNA tissue specificity score"), errors="coerce")
    df["hpa_protein_specificity_category"] = df.get("Protein tissue specificity")
    df["hpa_protein_distribution_category"] = df.get("Protein tissue distribution")
    df["hpa_protein_specificity_score"] = pd.to_numeric(df.get("Protein tissue specificity score"), errors="coerce")
    return df[
        [
            "target_gene",
            "gene_symbol",
            "ensembl_id",
            "target_uniprot",
            "hpa_specificity_category",
            "hpa_distribution_category",
            "hpa_specificity_score",
            "hpa_protein_specificity_category",
            "hpa_protein_distribution_category",
            "hpa_protein_specificity_score",
        ]
    ].drop_duplicates("target_gene")


def _hpa_tissue_rows(raw_dir: Path, target_genes: set[str], ensembl_ids: set[str]) -> pd.DataFrame:
    frames = []
    consensus_path = raw_dir / "hpa_rna_tissue_consensus.tsv.zip"
    if consensus_path.exists():
        df = _read_zipped_tsv(consensus_path)
        df = _filter_targets(df, target_genes, ensembl_ids)
        if not df.empty:
            frames.append(
                pd.DataFrame(
                    {
                        "target_gene": df["Gene name"].map(_upper),
                        "ensembl_id": df["Gene"].map(_upper),
                        "tissue": df["Tissue"].map(_lower),
                        "source_name": "Human Protein Atlas",
                        "source_family": "HPA",
                        "expression_metric": "hpa_consensus_ntpm",
                        "expression_value": pd.to_numeric(df["nTPM"], errors="coerce"),
                    }
                )
            )
    gtex_path = raw_dir / "hpa_rna_tissue_gtex.tsv.zip"
    if gtex_path.exists():
        df = _read_zipped_tsv(gtex_path)
        df = _filter_targets(df, target_genes, ensembl_ids)
        if not df.empty:
            frames.append(
                pd.DataFrame(
                    {
                        "target_gene": df["Gene name"].map(_upper),
                        "ensembl_id": df["Gene"].map(_upper),
                        "tissue": df["Tissue"].map(_lower),
                        "source_name": "Human Protein Atlas GTEx",
                        "source_family": "GTEx",
                        "expression_metric": "gtex_tpm",
                        "expression_value": pd.to_numeric(df["TPM"], errors="coerce"),
                        "gtex_ntpm": pd.to_numeric(df["nTPM"], errors="coerce"),
                    }
                )
            )
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _bgee_rows(raw_dir: Path, target_genes: set[str], ensembl_ids: set[str]) -> pd.DataFrame:
    path = raw_dir / "bgee_homo_sapiens_expr_simple.tsv.gz"
    if not path.exists():
        return pd.DataFrame()
    frames = []
    usecols = [
        "Gene ID",
        "Gene name",
        "Anatomical entity name",
        "Expression",
        "Call quality",
        "FDR",
        "Expression score",
    ]
    with gzip.open(path, "rt") as handle:
        for chunk in pd.read_csv(handle, sep="\t", usecols=usecols, chunksize=500_000):
            mask = chunk["Gene name"].map(_upper).isin(target_genes) | chunk["Gene ID"].map(_upper).isin(ensembl_ids)
            chunk = chunk.loc[mask].copy()
            if chunk.empty:
                continue
            frames.append(
                pd.DataFrame(
                    {
                        "target_gene": chunk["Gene name"].map(_upper),
                        "ensembl_id": chunk["Gene ID"].map(_upper),
                        "tissue": chunk["Anatomical entity name"].map(_lower),
                        "source_name": "Bgee",
                        "source_family": "Bgee",
                        "expression_metric": "bgee_expression_call",
                        "expression_value": pd.to_numeric(chunk["Expression score"], errors="coerce"),
                        "bgee_expression_call": chunk["Expression"].map(_lower),
                        "bgee_call_quality": chunk["Call quality"].map(_lower),
                        "bgee_fdr": pd.to_numeric(chunk["FDR"], errors="coerce"),
                    }
                )
            )
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _opentargets_safety_context(path: Path, target_genes: set[str]) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if "target_id" not in df.columns:
        return pd.DataFrame()
    df = df.loc[df["target_id"].map(_upper).isin(target_genes)].copy()
    if df.empty:
        return pd.DataFrame()
    tissue_labels: dict[str, set[str]] = {}
    counts: dict[str, int] = {}
    for _, row in df.iterrows():
        gene = _upper(row.get("target_id"))
        counts[gene] = counts.get(gene, 0) + 1
        raw = row.get("_metadata_json")
        if not isinstance(raw, str) or not raw.strip():
            continue
        try:
            payload = json.loads(raw)
        except json.JSONDecodeError:
            continue
        for biosample in payload.get("biosamples") or []:
            label = _clean((biosample or {}).get("tissueLabel"))
            if label:
                tissue_labels.setdefault(gene, set()).add(label)
    rows = []
    for gene, count in counts.items():
        rows.append(
            {
                "target_gene": gene,
                "ot_safety_liability_count": count,
                "ot_safety_biosample_tissues": ";".join(sorted(tissue_labels.get(gene, set()))),
                "source_name": "Open Targets Safety",
                "source_family": "OpenTargets",
            }
        )
    return pd.DataFrame(rows)


def _opentargets_baseline_rows(raw_dir: Path, metadata: pd.DataFrame) -> pd.DataFrame:
    directory = raw_dir / "opentargets_baseline_expression"
    if not directory.exists() or metadata.empty:
        return pd.DataFrame()
    try:
        import pyarrow.parquet as pq
    except Exception:
        return pd.DataFrame()
    ensembl_ids = set(metadata["ensembl_id"].dropna().map(_upper))
    if not ensembl_ids:
        return pd.DataFrame()
    frames = []
    cols = [
        "targetId",
        "tissueBiosampleFromSource",
        "datasourceId",
        "median",
        "distribution_score",
        "specificity_score",
    ]
    for path in sorted(directory.glob("*.parquet")):
        table = pq.read_table(path, columns=cols)
        df = table.to_pandas()
        df = df.loc[df["targetId"].map(_upper).isin(ensembl_ids)].copy()
        if df.empty:
            continue
        frames.append(
            pd.DataFrame(
                {
                    "ensembl_id": df["targetId"].map(_upper),
                    "tissue": df["tissueBiosampleFromSource"].map(_lower),
                    "source_name": "Open Targets baseline expression",
                    "source_family": "OpenTargets",
                    "expression_metric": "ot_expression_value",
                    "expression_value": pd.to_numeric(df["median"], errors="coerce"),
                    "ot_expression_value": pd.to_numeric(df["median"], errors="coerce"),
                    "ot_expression_distribution_score": pd.to_numeric(df["distribution_score"], errors="coerce"),
                    "ot_expression_specificity_score": pd.to_numeric(df["specificity_score"], errors="coerce"),
                    "ot_expression_datasource": df["datasourceId"].map(_clean),
                }
            )
        )
    if not frames:
        return pd.DataFrame()
    rows = pd.concat(frames, ignore_index=True)
    return rows.merge(metadata[["target_gene", "ensembl_id"]].drop_duplicates("ensembl_id"), on="ensembl_id", how="left")


def _aggregate_expression(rows: pd.DataFrame, metadata: pd.DataFrame, safety: pd.DataFrame) -> pd.DataFrame:
    targets = metadata.copy() if not metadata.empty else pd.DataFrame(columns=["target_gene", "ensembl_id"])
    if targets.empty and not rows.empty:
        targets = rows[["target_gene", "ensembl_id"]].drop_duplicates("target_gene")
    agg_parts = []
    if not rows.empty:
        hpa = rows.loc[rows["expression_metric"].eq("hpa_consensus_ntpm")].copy()
        if not hpa.empty:
            idx = hpa.groupby("target_gene")["expression_value"].idxmax()
            max_rows = hpa.loc[idx, ["target_gene", "tissue", "expression_value"]].rename(
                columns={"tissue": "hpa_max_tissue", "expression_value": "hpa_consensus_ntpm"}
            )
            stats = hpa.groupby("target_gene")["expression_value"].agg(
                hpa_consensus_median_ntpm="median",
                hpa_detected_tissue_count=lambda s: int((pd.to_numeric(s, errors="coerce") > 0).sum()),
            )
            agg_parts.append(max_rows.merge(stats, on="target_gene", how="left"))
        gtex = rows.loc[rows["expression_metric"].eq("gtex_tpm")].copy()
        if not gtex.empty:
            idx = gtex.groupby("target_gene")["expression_value"].idxmax()
            max_rows = gtex.loc[idx, ["target_gene", "tissue", "expression_value"]].rename(
                columns={"tissue": "gtex_max_tissue", "expression_value": "gtex_median_tpm"}
            )
            agg_parts.append(max_rows)
        bgee = rows.loc[rows["expression_metric"].eq("bgee_expression_call")].copy()
        if not bgee.empty:
            present = bgee["bgee_expression_call"].eq("present")
            bgee_stats = bgee.assign(_present=present).groupby("target_gene").agg(
                bgee_expression_call=("_present", lambda s: "present" if bool(s.any()) else "absent"),
                bgee_expression_score=("expression_value", "max"),
                bgee_present_tissue_count=("_present", "sum"),
            ).reset_index()
            agg_parts.append(bgee_stats)
        ot = rows.loc[rows["expression_metric"].eq("ot_expression_value")].copy()
        if not ot.empty and "target_gene" in ot.columns:
            ot = ot.loc[ot["target_gene"].notna()].copy()
            if not ot.empty:
                idx = ot.groupby("target_gene")["ot_expression_value"].idxmax()
                max_rows = ot.loc[idx, ["target_gene", "tissue", "ot_expression_value", "ot_expression_specificity_score"]].rename(
                    columns={"tissue": "ot_expression_max_tissue", "ot_expression_specificity_score": "ot_expression_zscore"}
                )
                agg_parts.append(max_rows)
    out = targets.drop_duplicates("target_gene")
    for part in agg_parts:
        out = out.merge(part.drop_duplicates("target_gene"), on="target_gene", how="left")
    if not safety.empty:
        out = out.merge(safety.drop_duplicates("target_gene"), on="target_gene", how="left", suffixes=("", "_ot_safety"))
    present_sources = sorted(rows["source_family"].dropna().unique().tolist()) if not rows.empty else []
    if not safety.empty and "OpenTargets" not in present_sources:
        present_sources.append("OpenTargets")
    out["expression_source"] = ";".join(present_sources)
    return out.loc[:, ~out.columns.duplicated()].copy()


def stage_tissue_expression_sources(
    *,
    target_table_path: str | Path,
    out_dir: str | Path = "data/external/hpa_gtex",
    download: bool = True,
    download_opentargets: bool = False,
    opentargets_part_limit: int = 0,
    opentargets_safety_path: str | Path = "data/external/opentargets/safety.tsv",
) -> dict[str, Any]:
    out = Path(out_dir)
    raw = out / "raw"
    out.mkdir(parents=True, exist_ok=True)
    raw.mkdir(parents=True, exist_ok=True)
    downloads: dict[str, Any] = {}
    if download:
        downloads["expression_sources"] = download_expression_sources(raw)
    if download_opentargets:
        downloads["opentargets_baseline_expression"] = download_opentargets_baseline(
            raw, limit=opentargets_part_limit
        )

    target_table = load_target_gene_set(target_table_path)
    target_genes = set(target_table["target_gene"].dropna().astype(str))
    metadata = _load_hpa_metadata(raw, target_genes)
    ensembl_ids = set(metadata.get("ensembl_id", pd.Series(dtype="object")).dropna().map(_upper))
    source_rows = []
    hpa_rows = _hpa_tissue_rows(raw, target_genes, ensembl_ids)
    if not hpa_rows.empty:
        source_rows.append(hpa_rows)
    bgee_rows = _bgee_rows(raw, target_genes, ensembl_ids)
    if not bgee_rows.empty:
        source_rows.append(bgee_rows)
    ot_baseline = _opentargets_baseline_rows(raw, metadata)
    if not ot_baseline.empty:
        source_rows.append(ot_baseline)
    tissue_rows = pd.concat(source_rows, ignore_index=True) if source_rows else pd.DataFrame()
    safety = _opentargets_safety_context(Path(opentargets_safety_path), target_genes)
    by_target = _aggregate_expression(tissue_rows, metadata, safety)

    source_rows_path = out / "tissue_expression_source_rows.tsv"
    tissue_rows.to_csv(source_rows_path, sep="\t", index=False)
    out_path = out / "tissue_expression.tsv"
    by_target.to_csv(out_path, sep="\t", index=False)
    manifest = {
        "target_table_path": str(target_table_path),
        "out_dir": str(out),
        "raw_dir": str(raw),
        "downloads": downloads,
        "n_requested_targets": int(len(target_table)),
        "n_hpa_metadata_targets": int(metadata["target_gene"].nunique()) if not metadata.empty else 0,
        "n_source_rows": int(len(tissue_rows)),
        "n_aggregated_targets": int(len(by_target)),
        "source_rows_path": str(source_rows_path),
        "tissue_expression_path": str(out_path),
        "source_families": sorted(tissue_rows["source_family"].dropna().unique().tolist()) if not tissue_rows.empty else [],
        "opentargets_safety_targets": int(safety["target_gene"].nunique()) if not safety.empty else 0,
        "opentargets_baseline_parts_present": len(list((raw / "opentargets_baseline_expression").glob("*.parquet")))
        if (raw / "opentargets_baseline_expression").exists()
        else 0,
        "policy": "Expression context is deterministic target metadata, not a supervised label. Duplicate source rows are aggregated to one target-level row before ML joins.",
    }
    manifest_path = out / "tissue_expression.manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest

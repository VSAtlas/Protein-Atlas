from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


def load_spd_reference(spd_workbook: str | Path) -> dict[str, pd.DataFrame]:
    """Load the SPD drug, assay, and assay-result sheets with real headers."""

    path = Path(spd_workbook)
    drugs = pd.read_excel(path, sheet_name="S Data 2", header=2)
    assays = pd.read_excel(path, sheet_name="S Data 3", header=2)
    results = pd.read_excel(path, sheet_name="S Data 1", header=2)
    for frame in (drugs, assays, results):
        frame.columns = [str(col).strip() for col in frame.columns]
    return {"drugs": drugs, "assays": assays, "results": results}


def audit_pubchem_aid_spd_overlap(
    pubchem_aid_table: str | Path,
    spd_workbook: str | Path,
    out_dir: str | Path,
    *,
    target_map_path: str | Path | None = None,
) -> dict[str, Any]:
    """Compare a normalized PubChem AID table against the full local SPD workbook."""

    aid = pd.read_csv(pubchem_aid_table, sep="\t", low_memory=False)
    spd = load_spd_reference(spd_workbook)
    drugs = spd["drugs"]
    assays = spd["assays"]
    results = spd["results"]
    aid_inchikeys = set(aid.get("inchikey", pd.Series(dtype=str)).dropna().astype(str).str.upper())
    aid_targets = set(aid.get("target_id", pd.Series(dtype=str)).dropna().astype(str).str.upper())
    target_genes = set(aid_targets)
    if target_map_path is not None and Path(target_map_path).exists():
        target_map = pd.read_csv(target_map_path, low_memory=False)
        for target_col in ("target_uniprot", "target_id"):
            if target_col not in target_map.columns or "target_gene" not in target_map.columns:
                continue
            hits = target_map[target_map[target_col].astype(str).str.upper().isin(aid_targets)]
            target_genes.update(hits["target_gene"].dropna().astype(str).str.upper())
    spd_drug_col = "drugcentral inchikey"
    spd_result_drug_col = "inchi_key"
    spd_gene_col = "EntrezGeneSymbol"
    drug_overlap = drugs[drugs[spd_drug_col].astype(str).str.upper().isin(aid_inchikeys)].copy()
    result_drug_overlap = results[results[spd_result_drug_col].astype(str).str.upper().isin(aid_inchikeys)].copy()
    target_overlap = pd.DataFrame()
    result_target_overlap = pd.DataFrame()
    if spd_gene_col in assays.columns:
        # Prefer direct gene symbols when the normalized AID target can be mapped to gene-like IDs.
        target_overlap = assays[assays[spd_gene_col].astype(str).str.upper().isin(target_genes)].copy()
    if target_overlap.empty:
        # Fall back to assay-name matching for accession-only AIDs; this is an audit hint, not a canonical map.
        target_tokens = set()
        if "assay name" in assays.columns:
            for token in aid_targets:
                hits = assays[assays["assay name"].astype(str).str.upper().str.contains(token, na=False)]
                if not hits.empty:
                    target_tokens.add(token)
        if target_tokens:
            target_overlap = assays[assays["assay name"].astype(str).str.upper().apply(lambda v: any(t in v for t in target_tokens))].copy()
    if "assay_group_name" in results.columns:
        assay_names = set(target_overlap.get("assay name", pd.Series(dtype=str)).dropna().astype(str))
        result_target_overlap = results[results["assay_group_name"].astype(str).isin(assay_names)].copy()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    drug_overlap.to_csv(out / "spd_drug_overlap.csv", index=False)
    result_drug_overlap.to_csv(out / "spd_result_drug_overlap.csv", index=False)
    target_overlap.to_csv(out / "spd_target_overlap.csv", index=False)
    result_target_overlap.to_csv(out / "spd_result_target_overlap.csv", index=False)
    manifest: dict[str, Any] = {
        "pubchem_aid_table": str(pubchem_aid_table),
        "spd_workbook": str(spd_workbook),
        "n_pubchem_rows": int(len(aid)),
        "n_pubchem_unique_inchikeys": int(len(aid_inchikeys)),
        "n_pubchem_unique_targets": int(len(aid_targets)),
        "pubchem_targets": sorted(aid_targets),
        "mapped_target_genes": sorted(target_genes - aid_targets),
        "spd_drug_overlap_rows": int(len(drug_overlap)),
        "spd_result_rows_same_drug": int(len(result_drug_overlap)),
        "spd_target_overlap_rows": int(len(target_overlap)),
        "spd_result_rows_same_target": int(len(result_target_overlap)),
        "interpretation": "Drug and target overlap are reported separately; pair overlap requires both the same compound and same target context.",
    }
    (out / "pubchem_spd_overlap_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return manifest

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pandas as pd


CHEMBL_RE = re.compile(r"CHEMBL\d+", re.IGNORECASE)
UNIPROT_RE = re.compile(r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^[OPQ][0-9][A-Z0-9]{3}[0-9]$", re.IGNORECASE)


def _read_csv(path: Path, **kwargs: Any) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False, **kwargs)


def _clean_text(value: Any) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def _clean_adr(value: Any) -> str:
    text = _clean_text(value).replace("_", " ")
    return re.sub(r"\s+", " ", text).strip()


def _read_atlas_target_gene_map(pair_table_path: str | Path) -> dict[str, str]:
    path = Path(pair_table_path)
    if not path.exists():
        return {}
    usecols = [col for col in ["target_gene", "target_id", "target_uniprot"] if col in pd.read_csv(path, nrows=0).columns]
    if not usecols:
        return {}
    df = pd.read_csv(path, usecols=usecols, low_memory=False)
    mapping: dict[str, str] = {}
    for row in df.drop_duplicates().itertuples(index=False):
        data = row._asdict()
        target = _clean_text(data.get("target_uniprot")) or _clean_text(data.get("target_id"))
        gene = _clean_text(data.get("target_gene")).upper()
        if gene and target:
            mapping.setdefault(gene, target)
    return mapping


def _smiles_to_inchikey(smiles: Any) -> str | None:
    text = _clean_text(smiles)
    if not text:
        return None
    try:
        from rdkit import Chem
    except Exception:
        return None
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return None
    return Chem.MolToInchiKey(mol)


def _read_atlas_inchikey_map(mapping_path: str | Path) -> dict[str, str]:
    path = Path(mapping_path)
    if not path.exists():
        return {}
    header = pd.read_csv(path, nrows=0).columns
    usecols = [col for col in ["drug_id", "inchikey", "atlas_smiles"] if col in header]
    if "drug_id" not in usecols:
        return {}
    mapping: dict[str, str] = {}
    for chunk in pd.read_csv(path, usecols=usecols, chunksize=100_000, low_memory=False):
        if "inchikey" not in chunk.columns and "atlas_smiles" in chunk.columns:
            chunk["inchikey"] = chunk["atlas_smiles"].map(_smiles_to_inchikey)
        for drug_id, inchikey in zip(chunk["drug_id"], chunk.get("inchikey", pd.Series(dtype=object)), strict=False):
            key = _clean_text(inchikey)
            drug = _clean_text(drug_id)
            if key and drug:
                mapping.setdefault(key, drug)
    return mapping


def _chembl_to_atlas_from_bindingdb(bindingdb_bioactivity_path: str | Path) -> dict[str, str]:
    path = Path(bindingdb_bioactivity_path)
    if not path.exists():
        return {}
    mapping: dict[str, str] = {}
    usecols = ["drug_id", "bindingdb_ligand_name"]
    for chunk in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=200_000, low_memory=False):
        for drug_id, ligand_name in zip(chunk["drug_id"], chunk["bindingdb_ligand_name"], strict=False):
            if pd.isna(drug_id) or pd.isna(ligand_name):
                continue
            for chembl_id in CHEMBL_RE.findall(str(ligand_name)):
                mapping.setdefault(chembl_id.upper(), str(drug_id))
    return mapping


def _normalize_label(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    out = pd.Series(pd.NA, index=series.index, dtype="Int64")
    out = out.mask(numeric.eq(1), 1)
    out = out.mask(numeric.eq(0), 0)
    return out


def _empty_activity_table() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "drug_id",
            "target_id",
            "activity_outcome",
            "activity_value_nM",
            "activity_relation",
            "activity_type",
            "activity_units",
            "assay_id",
            "source",
            "document_ids",
            "publication_year",
            "source_drug_id",
            "source_target_id",
            "source_specific_activity_label",
            "source_label_policy",
        ]
    )


def _empty_target_adr_table() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "target_id",
            "adr_id",
            "adr_term",
            "label_state",
            "source",
            "confidence",
            "document_ids",
            "gene_symbol",
            "source_target_id",
            "source_adr_id",
            "source_label_policy",
        ]
    )


def stage_adrtarget_target_adr(
    adrtarget_root: str | Path,
    out_dir: str | Path,
    *,
    pair_table_path: str | Path = "data/pilotstudy/moe_experts/banana_four_state_with_scorch_backfill_model_ready.csv",
) -> pd.DataFrame:
    """Stage ADRtarget paper predictions as target-ADR positive mechanism evidence.

    ADRtarget rows are target-ADR associations, not drug-target activity labels.
    They are therefore kept in the target_adr_mechanism namespace and should be
    evaluated separately from direct binding/activity evidence.
    """

    root = Path(adrtarget_root)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = root / "data" / "predicted_ADR_target_final.txt"
    if not path.exists():
        staged = _empty_target_adr_table()
        staged.to_csv(out / "adrtarget_target_adr_evidence.tsv", sep="\t", index=False)
        return staged
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    raw.columns = [str(col).strip().strip('"') for col in raw.columns]
    gene_map = _read_atlas_target_gene_map(pair_table_path)
    gene = raw.get("Entrez Gene Symbol for target", pd.Series(pd.NA, index=raw.index)).astype(str).str.strip().str.upper()
    target = gene.map(gene_map).fillna(gene)
    adr = raw.get("MedDRA term name", raw.get("UMLS concept name for ADR", pd.Series(pd.NA, index=raw.index))).map(_clean_adr)
    source_adr = raw.get("MedDRA ID", raw.get("UMLS concept CUI for ADR", pd.Series(pd.NA, index=raw.index)))
    staged = pd.DataFrame(
        {
            "target_id": target,
            "adr_id": adr,
            "adr_term": adr,
            "label_state": 1,
            "source": "ADRtarget",
            "confidence": 0.75,
            "document_ids": "ADRtarget",
            "gene_symbol": gene,
            "source_target_id": raw.get("target", gene),
            "source_adr_id": source_adr,
            "source_label_policy": (
                "ADRtarget target-ADR associations are paper-derived mechanism positives. "
                "They are not biochemical activity labels and do not provide negatives."
            ),
        }
    )
    staged = staged[staged["target_id"].astype(str).str.len().gt(0) & staged["adr_id"].astype(str).str.len().gt(0)]
    staged = staged.drop_duplicates(["target_id", "adr_id", "source"])
    staged.to_csv(out / "adrtarget_target_adr_evidence.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "source": "ADRtarget",
                "path": str(path),
                "rows": int(len(raw)),
                "staged_rows": int(len(staged)),
                "mapped_to_atlas_target_rows": int(gene.isin(gene_map).sum()),
            }
        ]
    ).to_csv(out / "adrtarget_ingestion_audit.csv", index=False)
    return staged


def stage_doctor_tardis_target_adr(
    doctor_root: str | Path,
    out_dir: str | Path,
    *,
    pair_table_path: str | Path = "data/pilotstudy/moe_experts/banana_four_state_with_scorch_backfill_model_ready.csv",
) -> pd.DataFrame:
    """Stage DocTOR/T-ARDIS seed protein-ADR rows as curated target-ADR positives."""

    root = Path(doctor_root)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    node_path = root / "node_info_BIANA.txt"
    extracted = root.parent / "extracted"
    outdf_candidates = [extracted / "controlled_outdf", extracted / "community_outdf", root / "controlled_outdf"]
    outdf_path = next((path for path in outdf_candidates if path.exists()), None)
    if not node_path.exists() or outdf_path is None:
        staged = _empty_target_adr_table()
        staged.to_csv(out / "doctor_tardis_target_adr_evidence.tsv", sep="\t", index=False)
        return staged
    node = pd.read_csv(node_path, sep="\t", header=None, names=["biana_id", "uniprot"], low_memory=False)
    biana_to_uniprot = node.dropna().drop_duplicates("biana_id").set_index("biana_id")["uniprot"].astype(str).to_dict()
    raw = pd.read_csv(outdf_path, sep="\t", low_memory=False)
    seed = pd.to_numeric(raw.get("Seed", pd.Series(0, index=raw.index)), errors="coerce").fillna(0).eq(1)
    work = raw.loc[seed].copy()
    target = work["BIANA ID"].map(biana_to_uniprot)
    adr = work.get("PT", work.get("Side_effect", pd.Series(pd.NA, index=work.index))).map(_clean_adr)
    staged = pd.DataFrame(
        {
            "target_id": target,
            "adr_id": adr,
            "adr_term": adr,
            "label_state": 1,
            "source": "DocTOR/T-ARDIS",
            "confidence": 0.8,
            "document_ids": "DocTOR/T-ARDIS",
            "gene_symbol": pd.NA,
            "source_target_id": work["BIANA ID"],
            "source_adr_id": work.get("Side_effect", adr),
            "source_label_policy": (
                "Only DocTOR/T-ARDIS seed protein-ADR rows are staged as positives. "
                "Predicted non-seed rankings are not production truth labels."
            ),
        }
    )
    staged = staged[staged["target_id"].notna() & staged["adr_id"].astype(str).str.len().gt(0)]
    staged = staged.drop_duplicates(["target_id", "adr_id", "source"])
    staged.to_csv(out / "doctor_tardis_target_adr_evidence.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "source": "DocTOR/T-ARDIS",
                "path": str(outdf_path),
                "rows": int(len(raw)),
                "seed_rows": int(seed.sum()),
                "staged_rows": int(len(staged)),
                "mapped_biana_rows": int(target.notna().sum()),
            }
        ]
    ).to_csv(out / "doctor_tardis_ingestion_audit.csv", index=False)
    return staged


def stage_dtiam_sources(
    dtiam_root: str | Path,
    out_dir: str | Path,
    *,
    pair_table_path: str | Path = "data/pilotstudy/moe_experts/banana_four_state_with_scorch_backfill_model_ready.csv",
    atlas_compound_map_path: str | Path = "data/pilotstudy/source_overlap_deep_sources_v3/spd_drug_atlas_mapping_audit.csv",
) -> pd.DataFrame:
    """Stage DTIAM MoA DTI positives when SMILES and target genes map to Atlas IDs."""

    root = Path(dtiam_root)
    if (root / "DTIAM-main").exists():
        root = root / "DTIAM-main"
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    gene_map = _read_atlas_target_gene_map(pair_table_path)
    inchikey_map = _read_atlas_inchikey_map(atlas_compound_map_path)
    rows: list[pd.DataFrame] = []
    audit_rows: list[dict[str, Any]] = []
    for mode in ["activation", "inhibition"]:
        base = root / "data" / "moa" / mode
        dti_path = base / "dti.csv"
        smiles_path = base / "drug_smi.csv"
        target_path = base / "tar_gene.csv"
        if not (dti_path.exists() and smiles_path.exists() and target_path.exists()):
            audit_rows.append({"source": f"DTIAM:{mode}", "rows": 0, "staged_rows": 0, "reason": "missing_files"})
            continue
        dti = pd.read_csv(dti_path, sep=None, engine="python")
        smiles = pd.read_csv(smiles_path, sep=None, engine="python")
        targets = pd.read_csv(target_path, sep=None, engine="python")
        smiles["inchikey"] = smiles["smi"].map(_smiles_to_inchikey)
        smiles["drug_id"] = smiles["inchikey"].map(inchikey_map)
        targets["target_id"] = targets["Gene_symbol"].astype(str).str.upper().map(gene_map)
        work = dti.merge(smiles[["DrugID", "smi", "inchikey", "drug_id"]], on="DrugID", how="left").merge(
            targets[["TargetID", "Gene_symbol", "target_id"]],
            on="TargetID",
            how="left",
        )
        mapped = work["drug_id"].notna() & work["target_id"].notna()
        staged = pd.DataFrame(
            {
                "drug_id": work.loc[mapped, "drug_id"],
                "target_id": work.loc[mapped, "target_id"],
                "activity_outcome": 1,
                "activity_value_nM": pd.NA,
                "activity_relation": "=",
                "activity_type": f"DTIAM_{mode}",
                "activity_units": "derived_dti_label",
                "assay_id": f"DTIAM:{mode}",
                "source": f"DTIAM:{mode}",
                "document_ids": "DTIAM",
                "publication_year": 2020,
                "source_drug_id": work.loc[mapped, "DrugID"],
                "source_target_id": work.loc[mapped, "TargetID"],
                "source_specific_activity_label": 1,
                "source_label_policy": (
                    "DTIAM MoA rows are derived DTI positives. Use as source-aware training/sensitivity evidence, "
                    "not as independent ADR truth."
                ),
            }
        )
        if not staged.empty:
            rows.append(staged)
        audit_rows.append(
            {
                "source": f"DTIAM:{mode}",
                "rows": int(len(dti)),
                "mapped_drug_rows": int(work["drug_id"].notna().sum()),
                "mapped_target_rows": int(work["target_id"].notna().sum()),
                "staged_rows": int(len(staged)),
                "reason": "ok",
            }
        )
    staged_all = pd.concat(rows, ignore_index=True) if rows else _empty_activity_table()
    staged_all.to_csv(out / "dtiam_dti_activity.tsv", sep="\t", index=False)
    pd.DataFrame(audit_rows).to_csv(out / "dtiam_ingestion_audit.csv", index=False)
    return staged_all


def stage_ensdti_sources(
    ensdti_root: str | Path,
    out_dir: str | Path,
    *,
    bindingdb_bioactivity_path: str | Path = "data/external/bindingdb/bioactivity.tsv",
) -> pd.DataFrame:
    """Normalize EnsDTI benchmark tables into Atlas four-state-compatible rows.

    EnsDTI benchmark labels are treated as derived DTI evidence. Rows are useful
    only when a ligand identifier maps to an Atlas FDA ligand and the protein is
    already a UniProt-like target identifier.
    """

    root = Path(ensdti_root)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    chembl_to_atlas = _chembl_to_atlas_from_bindingdb(bindingdb_bioactivity_path)
    rows: list[pd.DataFrame] = []
    audit_rows: list[dict[str, Any]] = []
    for dataset_dir in sorted((root / "dataset").glob("*")):
        data_path = dataset_dir / "data.csv"
        if not data_path.exists():
            continue
        df = _read_csv(data_path)
        required = {"lig", "pro", "lab"}
        if not required.issubset(df.columns):
            audit_rows.append(
                {
                    "source": "EnsDTI",
                    "dataset": dataset_dir.name,
                    "path": str(data_path),
                    "rows": int(len(df)),
                    "staged_rows": 0,
                    "reason": "missing_required_columns",
                }
            )
            continue
        ligand = df["lig"].astype(str).str.upper()
        target = df["pro"].astype(str).str.upper()
        label = _normalize_label(df["lab"])
        mapped_drug = ligand.map(chembl_to_atlas)
        mapped = mapped_drug.notna() & target.map(lambda value: bool(UNIPROT_RE.match(value))) & label.notna()
        staged = pd.DataFrame(
            {
                "drug_id": mapped_drug[mapped],
                "target_id": target[mapped],
                "activity_outcome": label[mapped].astype(int),
                "activity_value_nM": pd.to_numeric(df.loc[mapped, "affinity"], errors="coerce")
                if "affinity" in df.columns
                else pd.NA,
                "activity_relation": "=",
                "activity_type": "derived_dti_label",
                "activity_units": "benchmark_label",
                "assay_id": "EnsDTI:" + dataset_dir.name,
                "source": "EnsDTI:" + dataset_dir.name,
                "document_ids": "EnsDTI",
                "publication_year": 2025,
                "source_drug_id": ligand[mapped],
                "source_target_id": target[mapped],
                "source_specific_activity_label": label[mapped].astype(int),
                "source_label_policy": (
                    "EnsDTI benchmark labels are derived DTI labels. Use as sensitivity/training evidence only "
                    "after source-holdout auditing; do not count as independent clinical ADR truth."
                ),
            }
        )
        if not staged.empty:
            rows.append(staged)
        audit_rows.append(
            {
                "source": "EnsDTI",
                "dataset": dataset_dir.name,
                "path": str(data_path),
                "rows": int(len(df)),
                "positive_rows": int(label.eq(1).sum()),
                "negative_rows": int(label.eq(0).sum()),
                "mapped_drug_rows": int(mapped_drug.notna().sum()),
                "uniprot_like_target_rows": int(target.map(lambda value: bool(UNIPROT_RE.match(value))).sum()),
                "staged_rows": int(len(staged)),
                "reason": "ok",
            }
        )
    staged_all = pd.concat(rows, ignore_index=True) if rows else _empty_activity_table()
    staged_path = out / "ensdti_dti_activity.tsv"
    staged_all.to_csv(staged_path, sep="\t", index=False)
    audit = pd.DataFrame(audit_rows)
    audit.to_csv(out / "ensdti_ingestion_audit.csv", index=False)
    return staged_all


def stage_mosedti_audit(
    mosedti_root: str | Path,
    out_dir: str | Path,
    *,
    bindingdb_bioactivity_path: str | Path = "data/external/bindingdb/bioactivity.tsv",
) -> pd.DataFrame:
    """Audit MoseDTI variable-data files and stage mappable rows if present."""

    root = Path(mosedti_root)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    chembl_to_atlas = _chembl_to_atlas_from_bindingdb(bindingdb_bioactivity_path)
    rows: list[pd.DataFrame] = []
    audit_rows: list[dict[str, Any]] = []
    var_data = root / "var_data"
    for path in sorted(var_data.glob("*.tsv")):
        if path.name.endswith(("_train.tsv", "_valid.tsv", "_test.tsv")):
            continue
        label = 0 if path.stem.endswith("_neg") else 1
        relation = path.stem.removesuffix("_neg")
        df = pd.read_csv(path, sep="\t", header=None, names=["compound", "relation", "gene"], low_memory=False)
        compound = df["compound"].astype(str).str.replace("Compound;;", "", regex=False).str.upper()
        target = df["gene"].astype(str).str.replace("Gene;;", "", regex=False)
        mapped_drug = compound.map(chembl_to_atlas)
        mapped = mapped_drug.notna() & target.str.upper().map(lambda value: bool(UNIPROT_RE.match(value)))
        staged = pd.DataFrame(
            {
                "drug_id": mapped_drug[mapped],
                "target_id": target[mapped].str.upper(),
                "activity_outcome": label,
                "activity_value_nM": pd.NA,
                "activity_relation": "=",
                "activity_type": f"mosedti_{relation}",
                "activity_units": "benchmark_label",
                "assay_id": "MoseDTI:" + relation,
                "source": "MoseDTI:" + relation,
                "document_ids": "MoseDTI",
                "publication_year": 2025,
                "source_drug_id": compound[mapped],
                "source_target_id": target[mapped],
                "source_specific_activity_label": label,
                "source_label_policy": (
                    "MoseDTI negatives are model/benchmark negatives. Use only as sensitivity evidence unless "
                    "source-specific label lineage is independently verified."
                ),
            }
        )
        if not staged.empty:
            rows.append(staged)
        audit_rows.append(
            {
                "source": "MoseDTI",
                "dataset": path.stem,
                "path": str(path),
                "rows": int(len(df)),
                "label": label,
                "mapped_drug_rows": int(mapped_drug.notna().sum()),
                "uniprot_like_target_rows": int(target.str.upper().map(lambda value: bool(UNIPROT_RE.match(value))).sum()),
                "staged_rows": int(len(staged)),
            }
        )
    staged_all = pd.concat(rows, ignore_index=True) if rows else _empty_activity_table()
    staged_all.to_csv(out / "mosedti_dti_activity.tsv", sep="\t", index=False)
    audit = pd.DataFrame(audit_rows)
    audit.to_csv(out / "mosedti_ingestion_audit.csv", index=False)
    return audit


def stage_paper_dti_sources(
    out_dir: str | Path,
    *,
    ensdti_root: str | Path = "data/external/paper_sources/ensdti/EnsDTI-master",
    mosedti_root: str | Path = "data/external/paper_sources/mosedti/MoseDTI-main",
    adrtarget_root: str | Path = "data/external/paper_sources/adrtarget/ADRtarget-master",
    tardis_root: str | Path = "data/external/paper_sources/tardis/T-ARDIS-master",
    doctor_root: str | Path = "data/external/paper_sources/doctor/DocTOR-master",
    dtiam_root: str | Path = "data/external/paper_sources/dtiam/github/DTIAM-main",
    bindingdb_bioactivity_path: str | Path = "data/external/bindingdb/bioactivity.tsv",
    pair_table_path: str | Path = "data/pilotstudy/moe_experts/banana_four_state_with_scorch_backfill_model_ready.csv",
    atlas_compound_map_path: str | Path = "data/pilotstudy/source_overlap_deep_sources_v3/spd_drug_atlas_mapping_audit.csv",
) -> dict[str, Any]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ens = stage_ensdti_sources(
        ensdti_root,
        out,
        bindingdb_bioactivity_path=bindingdb_bioactivity_path,
    )
    mose = stage_mosedti_audit(
        mosedti_root,
        out,
        bindingdb_bioactivity_path=bindingdb_bioactivity_path,
    )
    mosedti_rows = int(len(pd.read_csv(out / "mosedti_dti_activity.tsv", sep="\t", low_memory=False)))
    adrtarget = stage_adrtarget_target_adr(adrtarget_root, out, pair_table_path=pair_table_path)
    doctor_root_path = Path(doctor_root)
    if not doctor_root_path.exists():
        doctor_root_path = Path(tardis_root)
    doctor = stage_doctor_tardis_target_adr(doctor_root_path, out, pair_table_path=pair_table_path)
    dtiam = stage_dtiam_sources(
        dtiam_root,
        out,
        pair_table_path=pair_table_path,
        atlas_compound_map_path=atlas_compound_map_path,
    )
    source_status = pd.DataFrame(
        [
            {
                "source": "ADRtarget",
                "local_root": str(adrtarget_root),
                "ingested_rows": int(len(adrtarget)),
                "status": "staged_target_adr_rows" if len(adrtarget) else "audited_no_current_atlas_overlap",
                "note": "Paper-linked predicted_ADR_target_final.txt staged as target-ADR positives; no negatives supplied.",
            },
            {
                "source": "DocTOR/T-ARDIS",
                "local_root": str(doctor_root_path),
                "ingested_rows": int(len(doctor)),
                "status": "staged_seed_target_adr_rows" if len(doctor) else "audited_no_current_atlas_overlap",
                "note": "Only seed protein-ADR rows from DocTOR outputs are staged as curated positives; model-ranked non-seeds are not truth labels.",
            },
            {
                "source": "DTIAM",
                "local_root": str(dtiam_root),
                "ingested_rows": int(len(dtiam)),
                "status": "staged_mappable_rows" if len(dtiam) else "audited_no_current_atlas_overlap",
                "note": "DTIAM GitHub MoA DTI rows are staged when SMILES and target genes map to Atlas IDs.",
            },
            {
                "source": "MoseDTI",
                "local_root": str(mosedti_root),
                "ingested_rows": mosedti_rows,
                "status": "audited_no_current_atlas_overlap" if mose["staged_rows"].sum() == 0 else "staged_mappable_rows",
                "note": "Variable-data files use benchmark labels; current mapping produced no Atlas-overlap rows.",
            },
            {
                "source": "EnsDTI",
                "local_root": str(ensdti_root),
                "ingested_rows": int(len(ens)),
                "status": "staged_mappable_rows" if len(ens) else "audited_no_current_atlas_overlap",
                "note": "Only rows mapping through BindingDB CHEMBL mentions to Atlas ligand IDs and UniProt-like targets are staged.",
            },
        ]
    )
    status_path = out / "paper_dti_source_status.csv"
    source_status.to_csv(status_path, index=False)
    manifest = {
        "out_dir": str(out),
        "ensdti_rows": int(len(ens)),
        "mosedti_rows": mosedti_rows,
        "adrtarget_target_adr_rows": int(len(adrtarget)),
        "doctor_tardis_target_adr_rows": int(len(doctor)),
        "dtiam_rows": int(len(dtiam)),
        "status_table": str(status_path),
        "normalized_activity_files": [
            str(out / "ensdti_dti_activity.tsv"),
            str(out / "mosedti_dti_activity.tsv"),
            str(out / "dtiam_dti_activity.tsv"),
        ],
        "normalized_target_adr_files": [
            str(out / "adrtarget_target_adr_evidence.tsv"),
            str(out / "doctor_tardis_target_adr_evidence.tsv"),
        ],
        "policy": (
            "Paper DTI benchmark labels are not clinical ADR truth. Use source-holdout and report as "
            "DTI/binding sensitivity evidence unless source-specific production semantics are verified."
        ),
    }
    (out / "paper_dti_source_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return manifest

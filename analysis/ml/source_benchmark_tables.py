from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.reporting.ligand_annotation_groups import resolve_ligand_annotations
from analysis.ml.build_ml_dataset import add_standard_ml_columns, deduplicate_ml_rows
from analysis.ml.labels import binary_label_series, truthy_series
from analysis.ml.year_metadata import add_availability_year_metadata


ID_COLS = [
    "drug_id",
    "target_id",
    "pdb_id",
    "ligand_chemotype",
    "scaffold_key",
    "protein_class",
    "target_gene",
    "target_uniprot",
    "ligand_base",
    "display_name",
    "generic_name",
    "ligand_chemotype_source",
    "scaffold_source",
]
FEATURE_COLS = [
    "atlas_score",
    "consensus_score",
    "free_cmax_um",
    "free_cmax_uM",
    "cmax_um",
    "fraction_unbound_plasma",
]
MODEL_READY_FEATURE_COLS = ["consensus_z_score", "atlas_score", "consensus_score"]
MODEL_READY_PROVENANCE_COLS = [
    "atlas_score_source_for_ml",
    "consensus_z_score_source",
    "final_score_source",
    "z_selected_source",
    "label_source",
    "source_family",
    "upstream_source",
    "source_label_policy",
    "assay_type",
    "endpoint_type",
    "activity_type",
    "assay_count",
    "activity_publication_year",
    "activity_publication_year_source",
    "database_release_year",
    "database_release_year_source",
    "document_year",
    "evidence_publication_year",
    "source_available_date",
    "source_release_date",
    "availability_year",
    "availability_year_source",
    "availability_year_confidence",
]
FDA_MAPPING_PATH = Path("chemdb/data/fda_mapping_from_pdbqt.csv")
UNASSIGNED_CHEMOTYPE = "Other / Unassigned"
SCAFFOLD_ASSIGNED_CHEMOTYPE = "Other / Structural scaffold assigned"


def _truth(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(False, index=df.index)
    return truthy_series(df[col]).fillna(False)


def _source_label(df: pd.DataFrame, prefix: str) -> tuple[pd.Series, pd.Series]:
    active = _truth(df, f"{prefix}_active")
    inactive = _truth(df, f"{prefix}_inactive")
    conflict = _truth(df, f"{prefix}_conflict") | (active & inactive)
    label = pd.Series(pd.NA, index=df.index, dtype="Int64")
    label = label.mask(inactive & ~conflict, 0)
    label = label.mask(active & ~conflict, 1)
    return label, conflict


def _keep(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    return df[[col for col in cols if col in df.columns]].copy()


def _source_metadata(prefixes: list[str]) -> list[str]:
    cols: list[str] = []
    for prefix in prefixes:
        cols.extend(
            [
                f"{prefix}_activity_nM",
                f"{prefix}_activity_type",
                f"{prefix}_endpoint_type",
                f"{prefix}_assay_type",
                f"{prefix}_assay_mode",
                f"{prefix}_activity_relation",
                f"{prefix}_assay_ids",
                f"{prefix}_assay_count",
                f"{prefix}_source",
                f"{prefix}_activity_document_ids",
                f"{prefix}_activity_publication_year_min",
                f"{prefix}_activity_publication_year_max",
                f"{prefix}_database_release_year",
                f"{prefix}_label_status",
                f"{prefix}_threshold_nM",
                f"{prefix}_missing_reason",
            ]
        )
    return cols


def _write_summary(path: Path, df: pd.DataFrame, label_col: str, extra: dict[str, Any] | None = None) -> dict[str, Any]:
    labels = binary_label_series(df[label_col]) if label_col in df.columns else pd.Series(pd.NA, index=df.index)
    summary: dict[str, Any] = {
        "path": str(path),
        "rows": int(len(df)),
        "label_col": label_col,
        "labelable_rows": int(labels.notna().sum()),
        "positive_rows": int(labels.eq(1).sum()),
        "negative_rows": int(labels.eq(0).sum()),
        "unknown_rows": int(labels.isna().sum()),
        "ambiguous_or_excluded_rows": int(labels.eq(-1).sum()),
        "positive_rate": float(labels.dropna().mean()) if labels.notna().any() else None,
    }
    if extra:
        summary.update(extra)
    return summary


def _write_model_ready(df: pd.DataFrame, label_col: str, out_path: Path) -> pd.DataFrame:
    seen: set[str] = set()
    keep = []
    for col in [*ID_COLS, "source_objective", *MODEL_READY_PROVENANCE_COLS, label_col, *MODEL_READY_FEATURE_COLS]:
        if col in df.columns and col not in seen:
            keep.append(col)
            seen.add(col)
    out = df[keep].copy()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    return out


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _norm_key(value: Any) -> str:
    return " ".join(str(value or "").strip().lower().split())


def _murcko_scaffold_key(smiles: str) -> str:
    if not str(smiles or "").strip():
        return ""
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold

        mol = Chem.MolFromSmiles(str(smiles).strip())
        if mol is None:
            return ""
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
        if not scaffold:
            scaffold = Chem.MolToSmiles(mol, isomericSmiles=False)
        digest = hashlib.sha1(scaffold.encode("utf-8")).hexdigest()[:12]
        return f"murcko:{digest}"
    except Exception:
        return ""


def _load_fda_mapping_scaffolds(root: Path) -> dict[str, dict[str, str]]:
    mapping_path = root / FDA_MAPPING_PATH
    if not mapping_path.exists():
        return {}
    wanted = {
        "path",
        "generic_name",
        "display_name",
        "pubchem_name",
        "drugcentral_generic_name",
        "smiles_neutral",
        "smiles",
        "inchikey",
    }
    try:
        mapping = pd.read_csv(
            mapping_path,
            low_memory=False,
            usecols=lambda col: col in wanted,
        )
    except Exception:
        return {}
    lookup: dict[str, dict[str, str]] = {}
    for row in mapping.to_dict("records"):
        smiles = str(row.get("smiles_neutral") or row.get("smiles") or "").strip()
        scaffold_key = _murcko_scaffold_key(smiles)
        if not scaffold_key:
            continue
        record = {
            "scaffold_key": scaffold_key,
            "inchikey": str(row.get("inchikey") or "").strip(),
            "scaffold_source": "fda_mapping_from_pdbqt_murcko",
        }
        path_text = str(row.get("path") or "").strip()
        if path_text:
            lookup[f"base:{Path(path_text).stem.lower()}"] = record
        for field in ("generic_name", "display_name", "pubchem_name", "drugcentral_generic_name"):
            key = _norm_key(row.get(field))
            if key:
                lookup[f"name:{key}"] = record
    return lookup


def _lookup_mapping_scaffold(
    lookup: dict[str, dict[str, str]],
    *,
    drug_id: str,
    ligand_base: str,
) -> dict[str, str]:
    for key in (
        f"base:{ligand_base.strip().lower()}",
        f"name:{_norm_key(drug_id)}",
    ):
        found = lookup.get(key)
        if found:
            return found
    return {}


def _ensure_ligand_chemotype(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "ligand_chemotype" not in out.columns:
        out["ligand_chemotype"] = ""
    if "scaffold_key" not in out.columns:
        out["scaffold_key"] = ""
    if "ligand_chemotype_source" not in out.columns:
        out["ligand_chemotype_source"] = ""
    if "scaffold_source" not in out.columns:
        out["scaffold_source"] = ""

    missing_chemotype = out["ligand_chemotype"].fillna("").astype(str).str.strip().isin(
        {"", UNASSIGNED_CHEMOTYPE}
    )
    missing_scaffold = out["scaffold_key"].fillna("").astype(str).str.strip().isin(
        {"", "other-unassigned"}
    )
    needs_annotation = missing_chemotype | missing_scaffold
    if not needs_annotation.any():
        return out

    root = _repo_root()
    mapping_scaffolds = _load_fda_mapping_scaffolds(root)
    cache: dict[tuple[str, str], dict[str, Any]] = {}
    for idx in out.index[needs_annotation]:
        drug_id = str(out.at[idx, "drug_id"]) if "drug_id" in out.columns else ""
        ligand_base = str(out.at[idx, "ligand_base"]) if "ligand_base" in out.columns else ""
        key = (drug_id.strip().lower(), ligand_base.strip().lower())
        if key not in cache:
            cache[key] = resolve_ligand_annotations(
                root,
                ligand_name=drug_id,
                ligand_base=ligand_base,
                library="fda",
            )
        resolved = cache[key]
        mapped = _lookup_mapping_scaffold(mapping_scaffolds, drug_id=drug_id, ligand_base=ligand_base)
        if missing_chemotype.at[idx]:
            chemotype = resolved.get("chemotype_primary", "") or UNASSIGNED_CHEMOTYPE
            if chemotype == UNASSIGNED_CHEMOTYPE and mapped:
                chemotype = SCAFFOLD_ASSIGNED_CHEMOTYPE
            out.at[idx, "ligand_chemotype"] = chemotype
            out.at[idx, "ligand_chemotype_source"] = (
                "fda_mapping_murcko"
                if chemotype == SCAFFOLD_ASSIGNED_CHEMOTYPE
                else "ligand_annotation_catalog"
            )
        if missing_scaffold.at[idx]:
            scaffold = mapped.get("scaffold_key") or resolved.get("scaffold_key", "") or "other-unassigned"
            out.at[idx, "scaffold_key"] = scaffold
            out.at[idx, "scaffold_source"] = mapped.get("scaffold_source") or "ligand_annotation_catalog"
    return out


def _finalize_table(df: pd.DataFrame, label_col: str, out_path: Path) -> pd.DataFrame:
    out = add_standard_ml_columns(df)
    out = _ensure_ligand_chemotype(out)
    out = add_availability_year_metadata(out)
    out[label_col] = binary_label_series(out[label_col])
    out = out.dropna(subset=[label_col]).copy()
    out[label_col] = out[label_col].astype(int)
    out = deduplicate_ml_rows(out, label_col)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    return out


def _curated_table(source: pd.DataFrame) -> pd.DataFrame:
    chembl_label, chembl_conflict = _source_label(source, "chembl")
    papyrus_label, papyrus_conflict = _source_label(source, "papyrus")
    disagreement = chembl_label.notna() & papyrus_label.notna() & chembl_label.ne(papyrus_label)
    label = pd.Series(pd.NA, index=source.index, dtype="Int64")
    label = label.mask(papyrus_label.notna(), papyrus_label)
    label = label.mask(chembl_label.notna(), chembl_label)
    label = label.mask(chembl_conflict | papyrus_conflict | disagreement, pd.NA)
    out = _keep(source, [*ID_COLS, *FEATURE_COLS, *_source_metadata(["chembl", "papyrus"])])
    out["curated_bioactivity_label"] = label
    out["source_objective"] = "curated_bioactivity"
    out["label_source"] = out["source_objective"]
    out["assay_type"] = source.get("chembl_assay_type", source.get("papyrus_assay_type", pd.Series(pd.NA, index=source.index)))
    out["endpoint_type"] = source.get("chembl_activity_type", source.get("papyrus_activity_type", pd.Series(pd.NA, index=source.index)))
    out["activity_type"] = out["endpoint_type"]
    out["source_family"] = "curated_bioactivity"
    out["upstream_source"] = "chembl;papyrus"
    out["source_conflict"] = (chembl_conflict | papyrus_conflict | disagreement).astype(bool)
    out["source_label_policy"] = "ChEMBL/Papyrus active or inactive; source conflicts/disagreements are unknown."
    if "chembl_activity_publication_year_min" in source.columns or "papyrus_activity_publication_year_min" in source.columns:
        years = pd.concat(
            [
                pd.to_numeric(source.get("chembl_activity_publication_year_min", pd.Series(pd.NA, index=source.index)), errors="coerce"),
                pd.to_numeric(source.get("papyrus_activity_publication_year_min", pd.Series(pd.NA, index=source.index)), errors="coerce"),
            ],
            axis=1,
        )
        out["activity_publication_year"] = years.min(axis=1, skipna=True)
    return out


def _toxcast_table(source: pd.DataFrame) -> pd.DataFrame:
    label, conflict = _source_label(source, "toxcast")
    out = _keep(source, [*ID_COLS, *FEATURE_COLS, *_source_metadata(["toxcast"])])
    out["toxcast_hts_label"] = label
    out["source_objective"] = "toxcast_hts"
    out["label_source"] = out["source_objective"]
    out["assay_type"] = source.get("toxcast_assay_type", pd.Series("HTS", index=source.index))
    out["endpoint_type"] = source.get("toxcast_endpoint_type", source.get("toxcast_activity_type", pd.Series("hit_call", index=source.index)))
    out["activity_type"] = out["endpoint_type"]
    out["source_family"] = "toxcast"
    out["upstream_source"] = "epa_invitrodb"
    out["source_conflict"] = conflict.astype(bool)
    out["source_label_policy"] = (
        "ToxCast HTS active or inactive after local ToxCast staging filters; conflicting active/inactive pairs are unknown."
    )
    if "toxcast_database_release_year" in source.columns:
        out["database_release_year"] = source["toxcast_database_release_year"]
    return out


def _spd_table(source: pd.DataFrame) -> pd.DataFrame:
    label_col = "ml_binary_label" if "ml_binary_label" in source.columns else "spd_exposure_relevant"
    out = _keep(
        source,
        [
            *ID_COLS,
            *FEATURE_COLS,
            "z_selected",
            "z_selected_source",
            "consensus_z_score",
            "consensus_z_score_source",
            "final_score",
            "final_score_source",
            "atlas_score_source_for_ml",
            "z_vs_decoys_consensus",
            "z_vs_compare_run_consensus",
            "z_vs_decoys_blend",
            "SCORCH_score_used",
            "exposure_margin",
            "spd_exposure_relevant",
            "spd_exposure_weak",
            "spd_exposure_unlikely",
            "spd_label_status",
            "spd_ac50_uM",
            "spd_assay_count",
            "spd_assay_ids",
            "spd_assay_name",
            "ml_binary_label",
            "ml_supervised_eligible",
            "ml_exclude_reason",
            "dedup_drug_key",
            "dedup_target_key",
        ],
    )
    if "drug_id" not in out.columns and "dedup_drug_key" in out.columns:
        out["drug_id"] = out["dedup_drug_key"]
    if "target_id" not in out.columns and "dedup_target_key" in out.columns:
        out["target_id"] = out["dedup_target_key"]
    out["spd_exposure_label"] = binary_label_series(source[label_col])
    out["source_objective"] = "spd_exposure_relevance"
    out["label_source"] = out["source_objective"]
    out["assay_type"] = "secondary_pharmacology"
    out["endpoint_type"] = "AC50/free_cmax_margin"
    out["activity_type"] = "exposure_relevance"
    out["source_family"] = "spd"
    out["upstream_source"] = "novartis_spd"
    out["source_label_policy"] = "SPD exposure relevance: AC50/free-Cmax <= 10 where available; missing SPD assays are unknown."
    out["database_release_year"] = 2023
    out["database_release_year_source"] = "Sutherland et al. 2023 SPD publication year; coarse dataset provenance, not row-level assay year"
    out["availability_year"] = 2023
    out["availability_year_source"] = "SPD public release/publication year; use only as coarse provenance, not prospective row-level evidence timing"
    return out


def _transfer_table(curated: pd.DataFrame, toxcast: pd.DataFrame) -> pd.DataFrame:
    curated_part = curated.rename(columns={"curated_bioactivity_label": "source_specific_activity_label"}).copy()
    toxcast_part = toxcast.rename(columns={"toxcast_hts_label": "source_specific_activity_label"}).copy()
    common = sorted(set(curated_part.columns) | set(toxcast_part.columns))
    combined = pd.concat(
        [
            curated_part.reindex(columns=common),
            toxcast_part.reindex(columns=common),
        ],
        ignore_index=True,
    )
    return combined


def _matched_transfer_table(transfer: pd.DataFrame) -> pd.DataFrame:
    curated = transfer[transfer["source_objective"].eq("curated_bioactivity")].copy()
    toxcast = transfer[transfer["source_objective"].eq("toxcast_hts")].copy()
    if curated.empty or toxcast.empty:
        return transfer.iloc[0:0].copy()
    mask = pd.Series(True, index=toxcast.index)
    for col in ("target_id", "ligand_chemotype"):
        if col in curated.columns and col in toxcast.columns:
            allowed = set(curated[col].dropna().astype(str))
            mask &= toxcast[col].fillna("").astype(str).isin(allowed)
    for col in ("atlas_score", "consensus_score"):
        if col in curated.columns and col in toxcast.columns:
            values = pd.to_numeric(curated[col], errors="coerce")
            if values.notna().sum() >= 10:
                low = values.quantile(0.01)
                high = values.quantile(0.99)
                mask &= pd.to_numeric(toxcast[col], errors="coerce").between(low, high)
    matched = pd.concat([curated, toxcast[mask]], ignore_index=True)
    matched["domain_match_policy"] = "target_id + ligand_chemotype overlap, with atlas/consensus score inside curated 1-99 percentile range."
    return matched


def build_source_benchmark_tables(
    bioactivity_source_path: str | Path,
    out_dir: str | Path,
    *,
    spd_table_path: str | Path | None = None,
) -> dict[str, Any]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    source = pd.read_csv(bioactivity_source_path, low_memory=False)
    curated = _finalize_table(_curated_table(source), "curated_bioactivity_label", out / "ml_curated_bioactivity_table.csv")
    toxcast = _finalize_table(_toxcast_table(source), "toxcast_hts_label", out / "ml_toxcast_hts_table.csv")
    transfer = _finalize_table(_transfer_table(curated, toxcast), "source_specific_activity_label", out / "ml_source_transfer_table.csv")
    matched = _finalize_table(_matched_transfer_table(transfer), "source_specific_activity_label", out / "ml_matched_source_transfer_table.csv")
    _write_model_ready(curated, "curated_bioactivity_label", out / "model_ready" / "ml_curated_bioactivity_model_ready.csv")
    _write_model_ready(toxcast, "toxcast_hts_label", out / "model_ready" / "ml_toxcast_hts_model_ready.csv")
    _write_model_ready(transfer, "source_specific_activity_label", out / "model_ready" / "ml_source_transfer_model_ready.csv")
    _write_model_ready(matched, "source_specific_activity_label", out / "model_ready" / "ml_matched_source_transfer_model_ready.csv")
    summaries = {
        "purpose": "Separate curated bioactivity, ToxCast HTS, and SPD exposure objectives so label definitions are not pooled silently.",
        "curated_bioactivity": _write_summary(out / "ml_curated_bioactivity_table.csv", curated, "curated_bioactivity_label"),
        "toxcast_hts": _write_summary(out / "ml_toxcast_hts_table.csv", toxcast, "toxcast_hts_label"),
        "source_transfer": _write_summary(out / "ml_source_transfer_table.csv", transfer, "source_specific_activity_label"),
        "matched_source_transfer": _write_summary(
            out / "ml_matched_source_transfer_table.csv",
            matched,
            "source_specific_activity_label",
            {"match_policy": "target_id + ligand_chemotype + score-range overlap"},
        ),
    }
    if spd_table_path is not None and Path(spd_table_path).exists():
        spd_source = pd.read_csv(spd_table_path, low_memory=False)
        spd = _finalize_table(_spd_table(spd_source), "spd_exposure_label", out / "ml_spd_exposure_table.csv")
        _write_model_ready(spd, "spd_exposure_label", out / "model_ready" / "ml_spd_exposure_model_ready.csv")
        summaries["spd_exposure"] = _write_summary(out / "ml_spd_exposure_table.csv", spd, "spd_exposure_label")
    else:
        summaries["spd_exposure"] = {"status": "missing", "reason": "SPD table path was not provided or does not exist."}
    (out / "source_benchmark_tables_manifest.json").write_text(
        json.dumps(summaries, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return summaries

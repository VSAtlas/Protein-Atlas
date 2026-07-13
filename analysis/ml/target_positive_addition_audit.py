from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.audit_utils import load_table


DEFAULT_FAMILIES = ("Protease", "Ion Channel", "Kinase", "Enzyme")
DEFAULT_LABEL = "combined_activity_ml_label"
DEFAULT_ACTIVE_NM = 1000.0
DEFAULT_FAMILY_POSITIVE_FLOOR = 50
DEFAULT_TARGET_POSITIVE_FLOOR = 5

DEFAULT_SOURCE_PATHS: dict[str, str] = {
    "BindingDB": "data/external/bindingdb/bioactivity.tsv",
    "ChEMBL": "data/external/chembl/bioactivity.tsv",
    "Papyrus": "data/external/papyrus/bioactivity.tsv",
    "TTD": "data/external/ttd/ttd_activity.tsv",
    "DrugCentral": "data/external/drugcentral/drugcentral_activity.tsv",
    "IUPHAR_GtoPdb": "data/external/iuphar_gtopdb/interactions.tsv",
    "CardiacSafety": "data/external/cardiac_safety/cardiac_ion_channel_bioactivity.tsv",
}

DIRECT_ENDPOINTS = {"ki", "kd", "ic50", "ec50", "ac50"}
HIGH_CONFIDENCE_ENDPOINTS = {"ki", "kd", "ic50"}
UPPER_BOUND_RELATIONS = {">", ">=", ">>"}


def audit_target_positive_additions(
    *,
    dataset_path: str | Path,
    out_dir: str | Path,
    label_col: str = DEFAULT_LABEL,
    families: Sequence[str] = DEFAULT_FAMILIES,
    source_paths: dict[str, str | Path] | None = None,
    active_nm: float = DEFAULT_ACTIVE_NM,
    family_positive_floor: int = DEFAULT_FAMILY_POSITIVE_FLOOR,
    target_positive_floor: int = DEFAULT_TARGET_POSITIVE_FLOOR,
    top_per_target: int = 25,
    include_existing_pairs: bool = False,
) -> dict[str, Any]:
    """Rank measured-positive additions for sparse Atlas target families.

    This audit does not create labels and does not mutate the model-ready table.
    It identifies source-backed measured positive drug-target candidates that
    overlap current Atlas/SPD targets, with enough provenance to decide whether
    a candidate should be docked, used as FDA-only evidence, or kept as a
    non-FDA sensitivity/probe row.
    """

    dataset = Path(dataset_path)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    base = load_table(dataset)
    if label_col not in base.columns:
        raise ValueError(f"Missing label column {label_col!r} in {dataset}")

    source_paths = {**DEFAULT_SOURCE_PATHS, **{k: str(v) for k, v in (source_paths or {}).items()}}
    focus_families = tuple(families or DEFAULT_FAMILIES)

    target_meta = _target_metadata(base)
    if focus_families:
        target_meta = target_meta[target_meta["target_family"].isin(focus_families)].copy()

    family_balance = _balance_table(base, label_col, ["target_family"], family_positive_floor)
    target_gap = _balance_table(
        base[base["target_uniprot"].isin(set(target_meta["target_uniprot"].dropna().astype(str)))].copy(),
        label_col,
        ["target_family", "target_gene", "target_uniprot", "pdb_id"],
        target_positive_floor,
    )
    target_gap["family_positive_floor"] = int(family_positive_floor)
    target_gap["target_positive_floor"] = int(target_positive_floor)
    target_gap["candidate_priority"] = target_gap.apply(_target_priority, axis=1)

    current_keys = _current_keys(base)
    source_frames = _load_sources(source_paths, active_nm=active_nm)
    candidates = _match_sources_to_targets(
        source_frames,
        target_meta=target_meta,
        current_keys=current_keys,
        active_nm=active_nm,
        include_existing_pairs=include_existing_pairs,
    )
    candidates = _attach_need_scores(candidates, target_gap, family_balance)
    candidates = _rank_candidates(candidates)

    recommended_pairs = (
        candidates.groupby(["target_gene", "target_uniprot", "pdb_id"], dropna=False, group_keys=False)
        .head(top_per_target)
        .reset_index(drop=True)
    )
    recommended_targets = _recommended_targets(target_gap, candidates)

    outputs = {
        "current_family_balance": out / "current_family_balance.csv",
        "target_positive_gap_summary": out / "target_positive_gap_summary.csv",
        "measured_positive_candidates": out / "measured_positive_candidates.csv",
        "recommended_target_positive_additions": out / "recommended_target_positive_additions.csv",
        "recommended_pair_additions": out / "recommended_pair_additions.csv",
    }
    family_balance.to_csv(outputs["current_family_balance"], index=False)
    target_gap.to_csv(outputs["target_positive_gap_summary"], index=False)
    candidates.to_csv(outputs["measured_positive_candidates"], index=False)
    recommended_targets.to_csv(outputs["recommended_target_positive_additions"], index=False)
    recommended_pairs.to_csv(outputs["recommended_pair_additions"], index=False)

    manifest = {
        "dataset": str(dataset),
        "out_dir": str(out),
        "label_col": label_col,
        "families": list(focus_families),
        "active_nm": float(active_nm),
        "family_positive_floor": int(family_positive_floor),
        "target_positive_floor": int(target_positive_floor),
        "top_per_target": int(top_per_target),
        "include_existing_pairs": bool(include_existing_pairs),
        "n_rows": int(len(base)),
        "n_focus_targets": int(target_meta["target_uniprot"].nunique()),
        "n_source_candidate_rows": int(len(candidates)),
        "n_recommended_pair_rows": int(len(recommended_pairs)),
        "source_paths": {key: str(path) for key, path in source_paths.items()},
        "outputs": {key: str(path) for key, path in outputs.items()},
        "policy": {
            "measured_positive_definition": f"activity_nM <= {active_nm:g} with non-upper-bound relation, or explicit source active label",
            "candidate_rows_are_not_labels": True,
            "non_fda_or_unverified_candidates_are_sensitivity_only": True,
            "existing_pairs_are_excluded_by_default": True,
            "absence_of_candidate_source_rows_is_not_negative_evidence": True,
        },
    }
    (out / "target_positive_addition_audit_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True, default=str),
        encoding="utf-8",
    )
    _write_markdown_summary(out / "target_positive_addition_audit.md", manifest, family_balance, recommended_targets)
    return manifest


def _read_tsv(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False)


def _target_metadata(df: pd.DataFrame) -> pd.DataFrame:
    cols = [col for col in ["target_gene", "target_uniprot", "target_id", "pdb_id", "target_family", "protein_class"] if col in df.columns]
    meta = df[cols].copy()
    if "target_uniprot" not in meta.columns:
        meta["target_uniprot"] = meta.get("target_id", "")
    for col in ["target_gene", "target_uniprot", "pdb_id", "target_family", "protein_class"]:
        if col not in meta.columns:
            meta[col] = ""
        meta[col] = meta[col].astype(str).str.strip()
    meta = meta[meta["target_uniprot"].ne("") | meta["target_gene"].ne("")]
    return meta.drop_duplicates(["target_gene", "target_uniprot", "pdb_id"]).reset_index(drop=True)


def _balance_table(df: pd.DataFrame, label_col: str, group_cols: Sequence[str], floor: int) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    label = pd.to_numeric(df[label_col], errors="coerce")
    work = df.copy()
    work["_label_num"] = label
    for keys, group in work.groupby(list(group_cols), dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        positives = int(group["_label_num"].eq(1).sum())
        negatives = int(group["_label_num"].eq(0).sum())
        unknown = int(group["_label_num"].isna().sum())
        labelable = positives + negatives
        row = {col: key for col, key in zip(group_cols, keys, strict=False)}
        row.update(
            {
                "n_rows": int(len(group)),
                "n_labelable": labelable,
                "n_positive": positives,
                "n_negative": negatives,
                "n_unknown": unknown,
                "positive_rate": positives / labelable if labelable else math.nan,
                "positive_floor": int(floor),
                "positive_gap_to_floor": max(0, int(floor) - positives),
            }
        )
        rows.append(row)
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    return result.sort_values(["positive_gap_to_floor", "n_positive", "n_labelable"], ascending=[False, True, False]).reset_index(drop=True)


def _target_priority(row: pd.Series) -> str:
    positives = int(row.get("n_positive", 0) or 0)
    negatives = int(row.get("n_negative", 0) or 0)
    if positives == 0 and negatives >= 10:
        return "highest_no_positives_many_negatives"
    if positives < 3 and negatives >= 10:
        return "high_sparse_positives_many_negatives"
    if positives < int(row.get("target_positive_floor", DEFAULT_TARGET_POSITIVE_FLOOR)):
        return "moderate_below_target_floor"
    return "lower_already_has_target_floor"


def _current_keys(df: pd.DataFrame) -> dict[str, set[Any]]:
    def values(col: str) -> set[str]:
        if col not in df.columns:
            return set()
        return set(df[col].dropna().astype(str).str.strip().loc[lambda s: s.ne("")])

    target = values("target_uniprot") | values("target_id")
    drug = values("drug_id") | values("ligand_base")
    pair_drug = set()
    pair_inchikey = set()
    pair_smiles = set()
    if target and ("target_uniprot" in df.columns or "target_id" in df.columns):
        target_col = "target_uniprot" if "target_uniprot" in df.columns else "target_id"
        for drug_col in [col for col in ["drug_id", "ligand_base"] if col in df.columns]:
            pair_drug |= set(zip(df[target_col].astype(str), df[drug_col].astype(str)))
        if "inchikey" in df.columns:
            pair_inchikey = set(zip(df[target_col].astype(str), df["inchikey"].astype(str)))
        if "canonical_smiles" in df.columns:
            pair_smiles = set(zip(df[target_col].astype(str), df["canonical_smiles"].astype(str)))
        elif "smiles" in df.columns:
            pair_smiles = set(zip(df[target_col].astype(str), df["smiles"].astype(str)))
    return {
        "targets": target,
        "drugs": drug,
        "inchikeys": values("inchikey"),
        "smiles": values("canonical_smiles") | values("smiles"),
        "pair_drug": pair_drug,
        "pair_inchikey": pair_inchikey,
        "pair_smiles": pair_smiles,
    }


def _load_sources(source_paths: dict[str, str | Path], *, active_nm: float) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for source_name, source_path in source_paths.items():
        path = Path(source_path)
        if not path.exists():
            continue
        try:
            if source_name == "IUPHAR_GtoPdb":
                frame = _normalize_iuphar(path, active_nm=active_nm)
            else:
                frame = _normalize_generic_activity(path, source_name=source_name, active_nm=active_nm)
        except Exception as exc:  # pragma: no cover - manifest captures bad local source shape.
            frames.append(
                pd.DataFrame(
                    [
                        {
                            "source": source_name,
                            "source_path": str(path),
                            "source_load_error": str(exc),
                        }
                    ]
                )
            )
            continue
        frames.append(frame)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True, sort=False)


def _normalize_generic_activity(path: Path, *, source_name: str, active_nm: float) -> pd.DataFrame:
    raw = _read_tsv(path)
    work = pd.DataFrame(index=raw.index)
    work["source"] = source_name
    work["source_path"] = str(path)
    work["drug_id"] = _first_existing(raw, ["drug_id", "ligand_id"])
    work["drug_name"] = _first_existing(raw, ["drug_name", "bindingdb_ligand_name", "ligand_name"])
    work["target_gene"] = _first_existing(raw, ["gene_symbol", "target_gene", "cardiac_gene_symbol"])
    work["target_uniprot"] = _first_existing(raw, ["uniprot", "target_uniprot", "target_id"])
    work["activity_nM"] = pd.to_numeric(_first_existing(raw, ["activity_nM", "standard_value_nm"]), errors="coerce")
    work["activity_relation"] = _first_existing(raw, ["activity_relation", "standard_relation"]).astype(str).str.strip()
    work["activity_type"] = _first_existing(raw, ["activity_type", "standard_type", "endpoint_type"]).astype(str).str.upper()
    work["assay_id"] = _first_existing(raw, ["assay_id", "bindingdb_reactant_set_id"])
    work["publication_year"] = _first_existing(raw, ["publication_year", "activity_publication_year", "database_release_year"])
    work["pubchem_cid"] = pd.to_numeric(
        _first_existing(raw, ["pubchem_cid", "pubchem_cid_resolved"]),
        errors="coerce",
    ).astype("Int64")
    work["inchikey"] = _first_existing(raw, ["inchikey", "standard_inchi_key"])
    work["smiles"] = _first_existing(raw, ["smiles", "canonical_smiles"])
    work["source_specific_label"] = _first_existing(raw, ["source_specific_activity_label", "activity_class"])
    work["source_row_id"] = raw.index.astype(str)
    work = _measured_positive_rows(work, active_nm=active_nm)
    return work


def _normalize_iuphar(path: Path, *, active_nm: float) -> pd.DataFrame:
    raw = _read_tsv(path)
    work = pd.DataFrame(index=raw.index)
    work["source"] = "IUPHAR_GtoPdb"
    work["source_path"] = str(path)
    work["drug_id"] = ""
    work["drug_name"] = _first_existing(raw, ["Ligand"])
    work["target_gene"] = _first_existing(raw, ["Target Gene Symbol"])
    work["target_uniprot"] = _first_existing(raw, ["Target UniProt ID"])
    work["activity_nM"] = pd.to_numeric(_first_existing(raw, ["Original Affinity Median nm", "Affinity Median"]), errors="coerce")
    work["activity_relation"] = _first_existing(raw, ["Original Affinity Relation"])
    work["activity_type"] = _first_existing(raw, ["Type"]).astype(str).str.upper()
    work["assay_id"] = _first_existing(raw, ["PubMed ID"])
    work["publication_year"] = ""
    work["pubchem_cid"] = pd.Series(pd.NA, index=raw.index, dtype="Int64")
    work["inchikey"] = ""
    work["smiles"] = ""
    work["source_specific_label"] = ""
    work["source_row_id"] = raw.index.astype(str)
    return _measured_positive_rows(work, active_nm=active_nm)


def _first_existing(df: pd.DataFrame, cols: Iterable[str]) -> pd.Series:
    for col in cols:
        if col in df.columns:
            return df[col]
    return pd.Series("", index=df.index)


def _measured_positive_rows(df: pd.DataFrame, *, active_nm: float) -> pd.DataFrame:
    work = df.copy()
    relation = work["activity_relation"].astype(str).str.strip()
    endpoint = work["activity_type"].astype(str).str.lower().str.strip()
    label = work["source_specific_label"].astype(str).str.lower().str.strip()
    value = pd.to_numeric(work["activity_nM"], errors="coerce")
    explicit_active = label.isin({"1", "1.0", "active", "positive", "strict_positive"})
    bound_safe = value.le(active_nm) & ~relation.isin(UPPER_BOUND_RELATIONS)
    direct_endpoint = endpoint.isin(DIRECT_ENDPOINTS) | endpoint.eq("")
    work = work[(explicit_active | bound_safe) & direct_endpoint].copy()
    work["measured_positive_policy"] = f"activity_nM <= {active_nm:g} and relation not upper-bound, or explicit active"
    work["endpoint_confidence"] = endpoint.map(lambda value: "high" if value in HIGH_CONFIDENCE_ENDPOINTS else ("medium" if value in DIRECT_ENDPOINTS else "unknown"))
    work["activity_uM"] = pd.to_numeric(work["activity_nM"], errors="coerce") / 1000.0
    return work


def _match_sources_to_targets(
    sources: pd.DataFrame,
    *,
    target_meta: pd.DataFrame,
    current_keys: dict[str, set[Any]],
    active_nm: float,
    include_existing_pairs: bool,
) -> pd.DataFrame:
    if sources.empty or target_meta.empty:
        return pd.DataFrame()

    target_by_uniprot = target_meta[target_meta["target_uniprot"].ne("")].copy()
    target_by_gene = target_meta[target_meta["target_gene"].ne("")].copy()
    left = sources.copy()
    for col in ["target_uniprot", "target_gene", "drug_id", "inchikey", "smiles"]:
        if col not in left.columns:
            left[col] = ""
        left[col] = left[col].astype(str).str.strip()

    matches: list[pd.DataFrame] = []
    if not target_by_uniprot.empty:
        by_uniprot = left.merge(
            target_by_uniprot,
            on="target_uniprot",
            how="inner",
            suffixes=("", "_atlas"),
        )
        by_uniprot["target_match_type"] = "uniprot"
        matches.append(by_uniprot)
    if not target_by_gene.empty:
        by_gene = left.merge(
            target_by_gene,
            on="target_gene",
            how="inner",
            suffixes=("", "_atlas"),
        )
        by_gene["target_match_type"] = "gene_symbol"
        matches.append(by_gene)
    if not matches:
        return pd.DataFrame()
    candidates = pd.concat(matches, ignore_index=True, sort=False)
    candidates = _canonicalize_matched_target_columns(candidates)
    candidates = candidates.drop_duplicates(
        ["source", "source_row_id", "target_gene", "target_uniprot", "pdb_id", "drug_id", "inchikey", "smiles"]
    )

    candidates["present_ligand_in_phase1"] = candidates.apply(lambda row: _present_ligand(row, current_keys), axis=1)
    candidates["present_pair_in_phase1"] = candidates.apply(lambda row: _present_pair(row, current_keys), axis=1)
    if not include_existing_pairs:
        candidates = candidates[~candidates["present_pair_in_phase1"]].copy()

    candidates["ligand_domain"] = candidates.apply(_ligand_domain, axis=1)
    candidates["training_domain"] = candidates["ligand_domain"].map(
        {
            "atlas_mapped_rdk": "fda_atlas_current_ligand",
            "external_mapped_current_ligand": "current_ligand_needs_mapping_review",
            "external_non_fda_or_unverified": "sensitivity_or_external_probe_only",
        }
    )
    candidates["production_truth_allowed"] = candidates["ligand_domain"].isin({"atlas_mapped_rdk", "external_mapped_current_ligand"})
    candidates["training_allowed"] = True
    candidates["active_nm_threshold"] = float(active_nm)
    candidates["ligand_base"] = candidates.apply(_stable_ligand_base, axis=1)
    candidates["display_name"] = candidates["drug_name"].where(
        candidates["drug_name"].astype(str).str.strip().ne(""),
        candidates["ligand_base"],
    )
    candidates["generic_name"] = candidates["display_name"]
    candidates["sources"] = candidates["source"]
    candidates["raw_sources"] = candidates["source"] + ":" + candidates["source_row_id"].astype(str)
    candidates["projected_label"] = 1
    candidates["projected_label_status"] = "strict_measured_positive"
    candidates["is_binding_positive_le1uM"] = True
    candidates["addon_label_policy"] = "measured_activity_positive_le_1uM"
    candidates["benchmark_only"] = False
    return candidates


def _stable_ligand_base(row: pd.Series) -> str:
    drug_id = str(row.get("drug_id", "") or "").strip()
    if drug_id:
        return drug_id
    source = re.sub(
        r"[^a-z0-9]+",
        "_",
        str(row.get("source", "external")).lower(),
    ).strip("_")
    identity = "|".join(
        str(row.get(field, "") or "").strip()
        for field in ("inchikey", "smiles", "drug_name", "source_row_id")
    )
    digest = hashlib.sha256(identity.encode("utf-8")).hexdigest()[:12]
    return f"{source or 'external'}:{digest}"


def _canonicalize_matched_target_columns(df: pd.DataFrame) -> pd.DataFrame:
    work = df.copy()
    for col in ["target_gene", "target_uniprot", "pdb_id", "target_family", "protein_class"]:
        atlas_col = f"{col}_atlas"
        if col not in work.columns:
            work[col] = ""
        if atlas_col in work.columns:
            atlas_values = work[atlas_col].astype(str).str.strip()
            current_values = work[col].astype(str).str.strip()
            work[col] = atlas_values.where(atlas_values.ne("") & atlas_values.ne("nan"), current_values)
    return work


def _present_ligand(row: pd.Series, current_keys: dict[str, set[Any]]) -> bool:
    drug_id = str(row.get("drug_id", "")).strip()
    inchikey = str(row.get("inchikey", "")).strip()
    smiles = str(row.get("smiles", "")).strip()
    return bool(
        (drug_id and drug_id in current_keys["drugs"])
        or (inchikey and inchikey in current_keys["inchikeys"])
        or (smiles and smiles in current_keys["smiles"])
    )


def _present_pair(row: pd.Series, current_keys: dict[str, set[Any]]) -> bool:
    target = str(row.get("target_uniprot", "")).strip()
    drug_id = str(row.get("drug_id", "")).strip()
    inchikey = str(row.get("inchikey", "")).strip()
    smiles = str(row.get("smiles", "")).strip()
    return bool(
        (target and drug_id and (target, drug_id) in current_keys["pair_drug"])
        or (target and inchikey and (target, inchikey) in current_keys["pair_inchikey"])
        or (target and smiles and (target, smiles) in current_keys["pair_smiles"])
    )


def _ligand_domain(row: pd.Series) -> str:
    drug_id = str(row.get("drug_id", "")).strip()
    if drug_id.startswith("rdk_"):
        return "atlas_mapped_rdk"
    if bool(row.get("present_ligand_in_phase1", False)):
        return "external_mapped_current_ligand"
    return "external_non_fda_or_unverified"


def _attach_need_scores(candidates: pd.DataFrame, target_gap: pd.DataFrame, family_balance: pd.DataFrame) -> pd.DataFrame:
    if candidates.empty:
        return candidates
    target_need_cols = [
        "target_family",
        "target_gene",
        "target_uniprot",
        "pdb_id",
        "n_positive",
        "n_negative",
        "positive_rate",
        "positive_gap_to_floor",
        "candidate_priority",
    ]
    merged = candidates.merge(
        target_gap[[col for col in target_need_cols if col in target_gap.columns]].rename(
            columns={
                "n_positive": "target_current_positives",
                "n_negative": "target_current_negatives",
                "positive_rate": "target_positive_rate",
                "positive_gap_to_floor": "target_positive_gap",
            }
        ),
        on=["target_family", "target_gene", "target_uniprot", "pdb_id"],
        how="left",
    )
    family_cols = ["target_family", "n_positive", "n_negative", "positive_rate", "positive_gap_to_floor"]
    merged = merged.merge(
        family_balance[[col for col in family_cols if col in family_balance.columns]].rename(
            columns={
                "n_positive": "family_current_positives",
                "n_negative": "family_current_negatives",
                "positive_rate": "family_positive_rate",
                "positive_gap_to_floor": "family_positive_gap",
            }
        ),
        on="target_family",
        how="left",
    )
    merged["audit_priority_score"] = merged.apply(_candidate_priority_score, axis=1)
    merged["selection_reason"] = merged.apply(_selection_reason, axis=1)
    return merged


def _candidate_priority_score(row: pd.Series) -> float:
    score = 0.0
    score += min(float(row.get("family_positive_gap", 0) or 0), 50.0) * 2.0
    score += min(float(row.get("target_positive_gap", 0) or 0), 10.0) * 3.0
    if row.get("ligand_domain") == "atlas_mapped_rdk":
        score += 20.0
    elif row.get("present_ligand_in_phase1"):
        score += 12.0
    if row.get("endpoint_confidence") == "high":
        score += 8.0
    elif row.get("endpoint_confidence") == "medium":
        score += 4.0
    value = pd.to_numeric(pd.Series([row.get("activity_nM")]), errors="coerce").iloc[0]
    if pd.notna(value):
        score += max(0.0, 10.0 - math.log10(max(float(value), 1e-6)))
    return score


def _selection_reason(row: pd.Series) -> str:
    bits: list[str] = []
    if float(row.get("family_positive_gap", 0) or 0) > 0:
        bits.append("family_positive_gap")
    if float(row.get("target_positive_gap", 0) or 0) > 0:
        bits.append("target_positive_gap")
    if row.get("ligand_domain") == "atlas_mapped_rdk":
        bits.append("current_fda_rdk_ligand")
    elif row.get("present_ligand_in_phase1"):
        bits.append("current_ligand_mapping_review")
    else:
        bits.append("non_fda_or_unverified_probe")
    if row.get("endpoint_confidence") == "high":
        bits.append("high_confidence_endpoint")
    return ";".join(bits)


def _rank_candidates(candidates: pd.DataFrame) -> pd.DataFrame:
    if candidates.empty:
        return candidates
    work = candidates.copy()
    work["_domain_rank"] = work["ligand_domain"].map(
        {
            "atlas_mapped_rdk": 0,
            "external_mapped_current_ligand": 1,
            "external_non_fda_or_unverified": 2,
        }
    ).fillna(3)
    work["_endpoint_rank"] = work["endpoint_confidence"].map({"high": 0, "medium": 1, "unknown": 2}).fillna(3)
    work["_activity_sort"] = pd.to_numeric(work["activity_nM"], errors="coerce")
    work = work.sort_values(
        [
            "audit_priority_score",
            "target_positive_gap",
            "_domain_rank",
            "_endpoint_rank",
            "_activity_sort",
            "source",
        ],
        ascending=[False, False, True, True, True, True],
        na_position="last",
    )
    dedupe_cols = ["target_uniprot", "pdb_id", "drug_id", "inchikey", "smiles"]
    work = work.drop_duplicates([col for col in dedupe_cols if col in work.columns], keep="first")
    return work.drop(columns=["_domain_rank", "_endpoint_rank", "_activity_sort"], errors="ignore").reset_index(drop=True)


def _recommended_targets(target_gap: pd.DataFrame, candidates: pd.DataFrame) -> pd.DataFrame:
    if target_gap.empty:
        return target_gap
    grouped = pd.DataFrame()
    if not candidates.empty:
        grouped = candidates.groupby(["target_gene", "target_uniprot", "pdb_id", "target_family"], dropna=False).agg(
            n_candidate_rows=("source", "size"),
            n_atlas_mapped_rdk_candidates=("ligand_domain", lambda s: int(s.eq("atlas_mapped_rdk").sum())),
            n_current_ligand_candidates=("present_ligand_in_phase1", lambda s: int(pd.Series(s).fillna(False).sum())),
            n_non_fda_or_unverified_candidates=("ligand_domain", lambda s: int(s.eq("external_non_fda_or_unverified").sum())),
            best_activity_nM=("activity_nM", "min"),
            candidate_sources=("source", lambda s: ";".join(sorted(set(map(str, s))))),
        ).reset_index()
    result = target_gap.merge(
        grouped,
        on=["target_gene", "target_uniprot", "pdb_id", "target_family"],
        how="left",
    )
    fill_zero = [
        "n_candidate_rows",
        "n_atlas_mapped_rdk_candidates",
        "n_current_ligand_candidates",
        "n_non_fda_or_unverified_candidates",
    ]
    for col in fill_zero:
        if col in result.columns:
            result[col] = result[col].fillna(0).astype(int)
    result["candidate_sources"] = result.get("candidate_sources", "").fillna("")
    result["recommended_action"] = result.apply(_recommended_action, axis=1)
    return result.sort_values(
        ["positive_gap_to_floor", "n_atlas_mapped_rdk_candidates", "n_current_ligand_candidates", "n_candidate_rows"],
        ascending=[False, False, False, False],
    ).reset_index(drop=True)


def _recommended_action(row: pd.Series) -> str:
    if int(row.get("positive_gap_to_floor", 0) or 0) <= 0:
        return "lower_priority_target_floor_met"
    if int(row.get("n_atlas_mapped_rdk_candidates", 0) or 0) > 0:
        return "dock_current_fda_rdk_ligands_first"
    if int(row.get("n_current_ligand_candidates", 0) or 0) > 0:
        return "review_mapping_then_dock_current_ligands"
    if int(row.get("n_candidate_rows", 0) or 0) > 0:
        return "non_fda_or_unverified_probe_sensitivity_only"
    return "no_local_measured_positive_candidates_found"


def _write_markdown_summary(
    path: Path,
    manifest: dict[str, Any],
    family_balance: pd.DataFrame,
    recommended_targets: pd.DataFrame,
) -> None:
    lines = [
        "# Target Positive Addition Audit",
        "",
        "This audit ranks measured-positive candidates for sparse Atlas/SPD target families.",
        "It does not relabel the model table and does not treat missing source evidence as negative.",
        "",
        f"- Dataset: `{manifest['dataset']}`",
        f"- Label: `{manifest['label_col']}`",
        f"- Candidate rows: {manifest['n_source_candidate_rows']}",
        f"- Recommended pair rows: {manifest['n_recommended_pair_rows']}",
        "",
        "## Family Balance",
        "",
    ]
    if not family_balance.empty:
        preview = family_balance.head(20)
        lines.extend(_markdown_table(preview))
    lines.extend(["", "## Top Target Gaps", ""])
    if not recommended_targets.empty:
        cols = [
            col
            for col in [
                "target_family",
                "target_gene",
                "pdb_id",
                "n_positive",
                "n_negative",
                "positive_gap_to_floor",
                "n_atlas_mapped_rdk_candidates",
                "n_current_ligand_candidates",
                "n_non_fda_or_unverified_candidates",
                "recommended_action",
            ]
            if col in recommended_targets.columns
        ]
        lines.extend(_markdown_table(recommended_targets[cols].head(30)))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _markdown_table(df: pd.DataFrame) -> list[str]:
    if df.empty:
        return ["No rows."]
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(_format_cell(row[col]) for col in cols) + " |")
    return lines


def _format_cell(value: Any) -> str:
    if pd.isna(value):
        return ""
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value).replace("|", "/")

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.external.spd import aggregate_spd_assays, normalize_spd_table


def _as_key(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def _normalize_target_gene(raw: pd.Series) -> pd.Series:
    def _extract_compact_gene(value: str) -> str:
        value = re.sub(r"\(.*?\)", " ", str(value))
        value = value.replace(",", " ").replace(";", " ").replace("/", " ")
        if not value:
            return ""
        for token in str(value).upper().split():
            token = token.strip(" ,;/()")
            if len(token) < 2 or _looks_like_assay_word(token):
                continue
            if not re.match(r"^[A-Z0-9][A-Z0-9_\-]{1,15}$", token):
                continue
            return token
        # legacy fallback: allow single-word tokens extracted from compact text
        token = re.sub(r"[^A-Z0-9_\-]", " ", str(value).upper())
        compact = [t for t in token.split() if re.match(r"^[A-Z0-9][A-Z0-9_\-]{1,15}$", t)]
        return compact[0] if compact else ""

    raw = raw.fillna("").astype(str).str.strip()
    out = []
    for value in raw:
        candidate = _extract_compact_gene(value)
        out.append(candidate)
    return pd.Series(out, index=raw.index, dtype=object)


def _looks_like_assay_word(token: str) -> bool:
    token_upper = token.upper()
    stop_words = {
        "INHIBITION",
        "BINDING",
        "AGONIST",
        "ANTAGONIST",
        "MODULATOR",
        "INCREASE",
        "DECREASE",
        "ACTIVATION",
        "SUPPRESSION",
        "POTENCY",
        "CONCENTRATION",
        "ACTIVITY",
    }
    return token_upper in stop_words


def _coerce_target_from_assay_text(raw: pd.Series) -> pd.Series:
    words = (
        raw.fillna("")
        .astype(str)
        .str.replace("/", " ", regex=False)
        .str.replace("\\", " ", regex=False)
        .str.replace(",", " ", regex=False)
        .str.replace(";", " ", regex=False)
        .str.split()
    )
    out = pd.Series("", index=raw.index, dtype="object")

    for idx, pieces in words.items():
        if not isinstance(pieces, list):
            out.iloc[idx] = ""
            continue
        candidate: str = ""
        for token in pieces:
            token = str(token).strip().upper()
            if not token or _looks_like_assay_word(token):
                continue
            if not re.match(r"^[A-Z0-9][A-Z0-9_\-]{1,15}$", token):
                continue
            # Prefer compact gene-like tokens; avoid long assay words such as
            # 'ACETYLCHOLINESTERASE' when this is only assay text.
            if len(token) > 12 and token.isalpha() and token.isupper():
                continue
            candidate = token
            break
        out.iloc[idx] = candidate
    return out


def _candidate_columns(df: pd.DataFrame, columns: list[str] | tuple[str, ...]) -> list[str]:
    canonical: dict[str, str] = {}
    normalized: dict[str, str] = {}
    for col in df.columns:
        raw_col = str(col)
        canonical[raw_col.strip().lower()] = raw_col
        normalized[re.sub(r"[^a-z0-9]+", "", raw_col.strip().lower())] = raw_col
    picked: list[str] = []
    for raw in columns:
        token = str(raw).strip()
        if token in canonical:
            picked.append(canonical[token])
            continue
        normalized_token = re.sub(r"[^a-z0-9]+", "", token.lower())
        if normalized_token in normalized:
            picked.append(normalized[normalized_token])
    return picked


def _first_nonempty(df: pd.DataFrame, columns: list[str] | tuple[str, ...]) -> pd.Series:
    out = pd.Series("", index=df.index, dtype="object")
    for col in _candidate_columns(df, columns):
        if col not in df.columns:
            continue
        value = _as_key(df[col])
        out = out.where(out.astype(str).str.len() > 0, value)
    return out


def _extract_target_id(df: pd.DataFrame, use_assay_fallback: bool = True) -> pd.Series:
    gene_cols = _first_nonempty(
        df,
        [
            "HumanEntrezGeneSymbol(representative)",
            "EntrezGeneSymbol",
            "gene_symbol",
            "target_id",
        ],
    )
    gene_id = _normalize_target_gene(gene_cols)
    if not use_assay_fallback:
        return gene_id

    missing = (gene_id == "") | gene_id.isna()
    if not missing.any():
        return gene_id
    assay_fallback = _coerce_target_from_assay_text(
        _first_nonempty(
            df,
            [
                "assay_group_name",
                "assay name",
                "assay_name",
            ],
        )
    )
    return gene_id.mask(missing, assay_fallback)


def _join_unique(values: pd.Series) -> str:
    normalized = [str(v).strip() for v in values if str(v).strip() and str(v).strip().lower() != "nan"]
    if not normalized:
        return ""
    return ";".join(sorted(set(normalized)))


def build_spd_target_assay_manifest(
    spd_path: str | Path,
    mapping_path: str | Path | None = None,
) -> pd.DataFrame:
    """Build per-target assay provenance rows from raw SPD assay tables."""
    raw = normalize_spd_table(spd_path, mapping_path)
    if raw.empty:
        return pd.DataFrame(
            columns=[
                "target_id",
                "assay_ids",
                "assay_names",
                "assay_count",
                "source_rows",
            ]
        )
    target_col = _extract_target_id(raw)
    raw["target_id"] = target_col
    raw["assay_id"] = _as_key(_first_nonempty(raw, ["assay_id", "assay group ID", "preferred assay ID", "assay_group_id"]))
    raw["assay_name"] = _as_key(
        _first_nonempty(raw, ["assay_name", "assay_group_name", "assay name", "name"])
    )
    raw["source"] = _as_key(raw.get("source", pd.Series("SPD", index=raw.index)))

    work = raw[raw["target_id"] != ""].copy()
    if work.empty:
        return pd.DataFrame(
            columns=[
                "target_id",
                "assay_ids",
                "assay_names",
                "assay_count",
                "source_rows",
            ]
        )
    rows = []
    for target_id, subset in work.groupby("target_id", dropna=False):
        rows.append(
            {
                "target_id": target_id,
                "assay_ids": _join_unique(subset["assay_id"]),
                "assay_names": _join_unique(subset["assay_name"]),
                "assay_count": int(len(subset)),
                "source_rows": int(len(subset)),
            }
        )
    return pd.DataFrame(rows).sort_values("target_id").reset_index(drop=True)


def _prepare_spd_for_panel(
    spd_path: str | Path,
    mapping_path: str | Path | None,
    margin_strong: float,
    margin_weak: float,
) -> pd.DataFrame:
    raw = normalize_spd_table(spd_path, mapping_path)
    raw["_spd_original_drug_id"] = raw.get("drug_id", pd.Series("", index=raw.index))
    raw["_spd_original_target_id"] = raw.get("target_id", pd.Series("", index=raw.index))
    raw["drug_id"] = _first_nonempty(raw, ["name", "drug_id"]).str.lower()
    raw["target_id"] = _extract_target_id(raw)
    raw["target_id"] = raw["target_id"].astype(str).str.upper()
    return aggregate_spd_assays(raw, margin_strong=margin_strong, margin_weak=margin_weak)


def stage_spd_full_panel(
    spd_path: str | Path,
    pair_table_path: str | Path,
    out_dir: str | Path,
    *,
    mapping_path: str | Path | None = None,
    margin_strong: float = 10.0,
    margin_weak: float = 100.0,
    include_assay_manifest: bool = False,
) -> dict[str, Any]:
    """Stage SPD assay coverage against the current Atlas pair table.

    This does not pretend to run missing dockings. It writes the complete SPD
    assayed drug-target panel, identifies the subset already present in Atlas,
    and emits run-list tables for missing SPD drugs/targets.
    """

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    spd = _prepare_spd_for_panel(spd_path, mapping_path, margin_strong, margin_weak)
    assay_manifest = (
        build_spd_target_assay_manifest(spd_path, mapping_path) if include_assay_manifest else pd.DataFrame()
    )
    pair = pd.read_csv(pair_table_path, low_memory=False)
    for col in ("drug_id", "target_id"):
        if col not in pair.columns:
            raise ValueError(f"pair table lacks required column: {col}")
        if col not in spd.columns:
            raise ValueError(f"SPD table lacks required normalized column: {col}")

    pair["_spd_panel_drug_key"] = _first_nonempty(
        pair,
        ["generic_name", "display_name", "mapped_drug_name", "ligand_display", "drug_id"],
    ).str.lower()
    pair["_spd_panel_target_key"] = _first_nonempty(
        pair,
        ["target_gene", "target_id", "target_uniprot"],
    ).str.upper()
    spd["_spd_panel_drug_key"] = _as_key(spd["drug_id"]).str.lower()
    spd["_spd_panel_target_key"] = _as_key(spd["target_id"]).str.upper()

    atlas_pairs = (
        pair[["_spd_panel_drug_key", "_spd_panel_target_key"]]
        .rename(columns={"_spd_panel_drug_key": "drug_match_key", "_spd_panel_target_key": "target_match_key"})
        .drop_duplicates()
    )
    atlas_drugs = set(atlas_pairs["drug_match_key"])
    atlas_targets = set(atlas_pairs["target_match_key"])
    spd["drug_match_key"] = spd["_spd_panel_drug_key"]
    spd["target_match_key"] = spd["_spd_panel_target_key"]
    coverage = spd.merge(
        atlas_pairs.assign(atlas_pair_available=True),
        on=["drug_match_key", "target_match_key"],
        how="left",
    )
    coverage["atlas_pair_available"] = coverage["atlas_pair_available"].eq(True)
    coverage["atlas_drug_available"] = coverage["drug_match_key"].isin(atlas_drugs)
    coverage["atlas_target_available"] = coverage["target_match_key"].isin(atlas_targets)
    coverage["atlas_missing_reason"] = ""
    coverage.loc[~coverage["atlas_drug_available"], "atlas_missing_reason"] = "missing_drug"
    coverage.loc[~coverage["atlas_target_available"], "atlas_missing_reason"] = coverage[
        "atlas_missing_reason"
    ].where(
        coverage["atlas_missing_reason"].eq(""),
        coverage["atlas_missing_reason"] + ";",
    ) + "missing_target"
    coverage.loc[
        coverage["atlas_drug_available"] & coverage["atlas_target_available"] & ~coverage["atlas_pair_available"],
        "atlas_missing_reason",
    ] = "missing_pair_docking"
    coverage.loc[coverage["atlas_pair_available"], "atlas_missing_reason"] = ""

    missing = coverage[~coverage["atlas_pair_available"]].copy()
    missing_targets = (
        missing[["target_id", "atlas_target_available"]]
        .drop_duplicates()
        .sort_values(["atlas_target_available", "target_id"])
    )
    missing_drugs = (
        missing[["drug_id", "atlas_drug_available"]]
        .drop_duplicates()
        .sort_values(["atlas_drug_available", "drug_id"])
    )

    panel_path = out / "spd_full_panel_assay_pairs.csv"
    coverage_path = out / "spd_full_panel_atlas_coverage.csv"
    missing_pairs_path = out / "spd_full_panel_missing_atlas_pairs.csv"
    missing_targets_path = out / "spd_full_panel_missing_targets.csv"
    missing_drugs_path = out / "spd_full_panel_missing_drugs.csv"
    assay_manifest_path = out / "spd_target_assay_manifest.csv"
    spd.to_csv(panel_path, index=False)
    coverage.to_csv(coverage_path, index=False)
    missing.to_csv(missing_pairs_path, index=False)
    missing_targets.to_csv(missing_targets_path, index=False)
    missing_drugs.to_csv(missing_drugs_path, index=False)
    if not assay_manifest.empty:
        assay_manifest.to_csv(assay_manifest_path, index=False)

    manifest: dict[str, Any] = {
        "status": "written",
        "role": "SPD full-panel coverage and Atlas run-list staging",
        "spd_assay_pairs": int(len(spd)),
        "spd_unique_drugs": int(spd["drug_id"].nunique()),
        "spd_unique_targets": int(spd["target_id"].nunique()),
        "atlas_pair_rows": int(len(pair)),
        "atlas_unique_drugs": int(pair["drug_id"].nunique()),
        "atlas_unique_targets": int(pair["target_id"].nunique()),
        "spd_pairs_already_in_atlas": int(coverage["atlas_pair_available"].sum()),
        "spd_pairs_missing_from_atlas": int((~coverage["atlas_pair_available"]).sum()),
        "spd_targets_missing_from_atlas": int((~missing_targets["atlas_target_available"]).sum()),
        "spd_drugs_missing_from_atlas": int((~missing_drugs["atlas_drug_available"]).sum()),
        "outputs": {
            "panel": str(panel_path),
            "coverage": str(coverage_path),
            "missing_pairs": str(missing_pairs_path),
            "missing_targets": str(missing_targets_path),
            "missing_drugs": str(missing_drugs_path),
            "assay_manifest": str(assay_manifest_path) if not assay_manifest.empty else "",
        },
        "next": (
            "Run Atlas docking/rescoring for the missing SPD targets and drugs before "
            "using SPD as a manuscript-scale exposure-relevance ML training set."
        ),
    }
    (out / "spd_full_panel_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return manifest

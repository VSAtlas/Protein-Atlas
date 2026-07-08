from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.adr_site_mapping import add_adr_site_mapping, adr_site_match_score

TPM_HIGH = 10.0
TPM_LOW = 1.0
SENSITIVE_SITE_TOKENS = ("heart", "liver", "kidney", "brain", "cns", "lung", "reproductive", "blood", "immune")

SPECIFICITY_SCORES = {
    "tissue enriched": 1.0,
    "tissue_enriched": 1.0,
    "group enriched": 0.75,
    "group_enriched": 0.75,
    "tissue enhanced": 0.5,
    "tissue_enhanced": 0.5,
    "low tissue specificity": 0.1,
    "low specificity": 0.1,
    "not detected": 0.0,
}
DISTRIBUTION_PENALTIES = {
    "detected in single": 0.0,
    "detected in some": 0.1,
    "detected in many": 0.35,
    "detected in all": 0.7,
    "single": 0.0,
    "some": 0.1,
    "many": 0.35,
    "all": 0.7,
    "not detected": 0.0,
}


def _read(path: str | Path | None) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame()
    p = Path(path)
    if not p.exists():
        return pd.DataFrame()
    sep = "\t" if p.suffix.lower() in {".tsv", ".tab"} else ","
    return pd.read_csv(p, sep=sep, low_memory=False)


def _clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _lower(value: Any) -> str:
    return _clean(value).lower()


def _numeric(df: pd.DataFrame, cols: list[str]) -> pd.Series:
    out = pd.Series(pd.NA, index=df.index, dtype="Float64")
    for col in cols:
        if col not in df.columns:
            continue
        vals = pd.to_numeric(df[col], errors="coerce").astype("Float64")
        out = out.where(out.notna(), vals)
    return out


def _category_score(df: pd.DataFrame, cols: list[str], mapping: dict[str, float]) -> pd.Series:
    out = pd.Series(pd.NA, index=df.index, dtype="Float64")
    for col in cols:
        if col not in df.columns:
            continue
        vals = df[col].map(lambda value: mapping.get(_lower(value), pd.NA)).astype("Float64")
        out = out.where(out.notna(), vals)
    return out


def _presence_from_tpm(values: pd.Series) -> pd.Series:
    out = pd.Series(pd.NA, index=values.index, dtype="Float64")
    out = out.mask(values.ge(TPM_HIGH).fillna(False), 1.0)
    out = out.mask((values.ge(TPM_LOW) & values.lt(TPM_HIGH)).fillna(False), 0.6)
    out = out.mask((values.gt(0) & values.lt(TPM_LOW)).fillna(False), 0.25)
    out = out.mask(values.eq(0).fillna(False), 0.0)
    return out


def _bgee_floor(df: pd.DataFrame) -> pd.Series:
    if "bgee_expression_call" not in df.columns:
        return pd.Series(pd.NA, index=df.index, dtype="Float64")
    calls = df["bgee_expression_call"].map(_lower)
    out = pd.Series(pd.NA, index=df.index, dtype="Float64")
    out = out.mask(calls.eq("present"), 0.6)
    out = out.mask(calls.eq("absent"), 0.25)
    return out


def _concordance_score(df: pd.DataFrame) -> pd.Series:
    present_votes = pd.DataFrame(index=df.index)
    hpa = _numeric(df, ["hpa_consensus_ntpm", "hpa_ntpm", "hpa_tpm"])
    gtex = _numeric(df, ["gtex_median_tpm", "gtex_tpm"])
    bgee_call = df.get("bgee_expression_call", pd.Series("", index=df.index)).map(_lower)
    if hpa.notna().any():
        present_votes["hpa"] = hpa.gt(0)
    if gtex.notna().any():
        present_votes["gtex"] = gtex.gt(0)
    if bgee_call.astype(bool).any():
        present_votes["bgee"] = bgee_call.eq("present")
    if present_votes.empty:
        return pd.Series(pd.NA, index=df.index, dtype="Float64")
    n_sources = present_votes.notna().sum(axis=1)
    n_present = present_votes.fillna(False).sum(axis=1)
    score = pd.Series(pd.NA, index=df.index, dtype="Float64")
    score = score.mask(n_sources.ge(3) & n_present.ge(3), 1.0)
    score = score.mask(n_sources.ge(2) & n_present.ge(2) & score.isna(), 0.7)
    score = score.mask(n_present.eq(1) & score.isna(), 0.4)
    score = score.mask((n_sources.ge(2) & n_present.eq(0)) & score.isna(), 0.2)
    return score


def _merge_expression(pair_table: pd.DataFrame, expression: pd.DataFrame) -> pd.DataFrame:
    if expression.empty:
        return pair_table.copy()
    left = pair_table.copy()
    right = expression.copy()
    key_sets = [
        ["target_id", "adr_site_group"],
        ["target_gene", "adr_site_group"],
        ["gene_symbol", "adr_site_group"],
        ["target_id"],
        ["target_gene"],
        ["gene_symbol"],
    ]
    for keys in key_sets:
        if set(keys).issubset(left.columns) and set(keys).issubset(right.columns):
            right = right.drop_duplicates(keys, keep="first")
            return left.merge(right, on=keys, how="left", suffixes=("", "_expression"))
    return left


def add_tissue_site_scores(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    out = add_adr_site_mapping(df)
    tpm = _numeric(
        out,
        [
            "hpa_consensus_ntpm",
            "hpa_ntpm",
            "hpa_tpm",
            "gtex_median_tpm",
            "gtex_tpm",
            "ot_expression_value",
        ],
    )
    presence = _presence_from_tpm(tpm)
    bgee = _bgee_floor(out)
    presence = presence.combine(bgee, lambda a, b: max([v for v in (a, b) if pd.notna(v)], default=pd.NA))
    if "expression_presence_score" not in out.columns:
        out["expression_presence_score"] = presence

    specificity = _category_score(
        out,
        ["hpa_specificity_category", "ot_specificity_category", "tissue_specificity_category"],
        SPECIFICITY_SCORES,
    )
    if "hpa_specificity_score" not in out.columns:
        out["hpa_specificity_score"] = specificity
    if "site_specificity_score" not in out.columns:
        out["site_specificity_score"] = specificity

    distribution = _category_score(
        out,
        ["hpa_distribution_category", "ot_distribution_category", "tissue_distribution_category"],
        DISTRIBUTION_PENALTIES,
    )
    if "distribution_penalty" not in out.columns:
        out["distribution_penalty"] = distribution

    concordance = _concordance_score(out)
    if "expression_concordance_score" not in out.columns:
        out["expression_concordance_score"] = concordance

    components = pd.DataFrame(
        {
            "presence": pd.to_numeric(out["expression_presence_score"], errors="coerce"),
            "specificity": pd.to_numeric(out["site_specificity_score"], errors="coerce"),
            "concordance": pd.to_numeric(out["expression_concordance_score"], errors="coerce"),
        }
    )
    weights = {"presence": 0.45, "specificity": 0.30, "concordance": 0.25}
    numer = sum(components[col].fillna(0.0) * weight for col, weight in weights.items())
    denom = sum(components[col].notna().astype(float) * weight for col, weight in weights.items())
    out["site_relevance_score"] = (numer / denom).where(denom.gt(0), pd.NA)

    if "adr_site_match_score" not in out.columns:
        out["adr_site_match_score"] = out.apply(adr_site_match_score, axis=1)
    adr_match = pd.to_numeric(out["adr_site_match_score"], errors="coerce")
    adr_confidence = pd.to_numeric(
        out.get("adr_site_mapping_confidence", pd.Series(pd.NA, index=out.index)),
        errors="coerce",
    )
    adr_adjustment = (0.7 * adr_match) + (0.3 * adr_confidence)
    adr_components = pd.DataFrame(
        {
            "target_expression": pd.to_numeric(out["site_relevance_score"], errors="coerce"),
            "adr_site_context": adr_adjustment,
        }
    )
    adr_weights = {"target_expression": 0.65, "adr_site_context": 0.35}
    adr_numer = sum(adr_components[col].fillna(0.0) * weight for col, weight in adr_weights.items())
    adr_denom = sum(adr_components[col].notna().astype(float) * weight for col, weight in adr_weights.items())
    out["site_relevance_score_by_adr"] = (adr_numer / adr_denom).where(adr_denom.gt(0), pd.NA)

    site_text = out.get("site_name", out.get("adr_site_group", pd.Series("", index=out.index))).map(_lower)
    sensitive = site_text.map(
        lambda text: any(token in text for token in SENSITIVE_SITE_TOKENS) if text else pd.NA
    ).astype("Float64")
    if "sensitive_offsite_expression_score" not in out.columns:
        out["sensitive_offsite_expression_score"] = sensitive
    offsite = pd.to_numeric(
        out.get("sensitive_offsite_expression_score", sensitive), errors="coerce"
    ).astype("Float64")
    if "ot_safety_liability_count" in out.columns:
        safety = pd.to_numeric(out["ot_safety_liability_count"], errors="coerce")
        safety_flag = safety.gt(0).where(safety.notna(), pd.NA).astype("Float64")
    else:
        safety_flag = pd.Series(pd.NA, index=out.index, dtype="Float64")
    adr_overlap = pd.to_numeric(out.get("adr_site_overlap_score", pd.Series(pd.NA, index=out.index)), errors="coerce")
    safety_components = pd.DataFrame(
        {
            "distribution": pd.to_numeric(out["distribution_penalty"], errors="coerce"),
            "offsite": offsite,
            "safety": safety_flag,
            "adr": adr_overlap,
        }
    )
    safety_weights = {"distribution": 0.35, "offsite": 0.30, "safety": 0.20, "adr": 0.15}
    safety_numer = sum(safety_components[col].fillna(0.0) * weight for col, weight in safety_weights.items())
    safety_denom = sum(safety_components[col].notna().astype(float) * weight for col, weight in safety_weights.items())
    out["site_safety_risk_score"] = (safety_numer / safety_denom).where(safety_denom.gt(0), pd.NA)

    summary = {
        "rows": int(len(out)),
        "site_relevance_nonmissing": int(out["site_relevance_score"].notna().sum()),
        "site_relevance_by_adr_nonmissing": int(out["site_relevance_score_by_adr"].notna().sum()),
        "adr_site_group_nonmissing": int(out["adr_site_group"].notna().sum()) if "adr_site_group" in out.columns else 0,
        "site_safety_risk_nonmissing": int(out["site_safety_risk_score"].notna().sum()),
        "policy": "Deterministic tissue/site scoring from normal expression and curated ADR-site mappings; missing evidence remains missing.",
    }
    return out, summary

def build_tissue_site_scores(
    pair_table_path: str | Path,
    out_path: str | Path,
    *,
    target_expression_path: str | Path | None = None,
) -> pd.DataFrame:
    pairs = _read(pair_table_path)
    if pairs.empty:
        raise ValueError(f"pair table is empty or missing: {pair_table_path}")
    expression = _read(target_expression_path)
    merged = _merge_expression(pairs, expression)
    scored, summary = add_tissue_site_scores(merged)
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    scored.to_csv(path, index=False)
    summary.update(
        {
            "pair_table": str(pair_table_path),
            "target_expression_table": str(target_expression_path) if target_expression_path else None,
            "out_path": str(path),
        }
    )
    path.with_suffix(".manifest.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    return scored

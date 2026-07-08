from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd


@dataclass(frozen=True)
class SiteRule:
    site_group: str
    site_name: str
    terms: tuple[str, ...]
    source: str
    confidence: float


ADR_SITE_RULES: tuple[SiteRule, ...] = (
    SiteRule(
        "heart",
        "heart",
        (
            "qt",
            "qtc",
            "torsade",
            "arrhythm",
            "ventricular",
            "cardiac",
            "heart",
            "myocard",
            "tachycard",
            "bradycard",
            "palpitation",
            "electrocardiogram",
        ),
        "adr_term_keyword",
        0.75,
    ),
    SiteRule(
        "liver",
        "liver",
        ("hepat", "liver", "bilirubin", "jaundice", "cholest", "transaminase", "alanine aminotransferase", "aspartate aminotransferase"),
        "adr_term_keyword",
        0.75,
    ),
    SiteRule(
        "brain_cns",
        "brain",
        (
            "brain",
            "cns",
            "central nervous",
            "sedation",
            "somnolence",
            "drowsiness",
            "hallucination",
            "psychosis",
            "seizure",
            "extrapyramidal",
            "akathisia",
            "dyskinesia",
            "cognitive",
            "depression",
            "anxiety",
            "dizziness",
        ),
        "adr_term_keyword",
        0.75,
    ),
    SiteRule(
        "gi",
        "gastrointestinal tract",
        (
            "gastrointestinal",
            "gastric",
            "stomach",
            "intestinal",
            "bowel",
            "nausea",
            "vomit",
            "diarrhea",
            "diarrhoea",
            "constipation",
            "ulcer",
            "gi bleeding",
            "gastrointestinal bleeding",
        ),
        "adr_term_keyword",
        0.75,
    ),
    SiteRule(
        "kidney",
        "kidney",
        ("renal", "kidney", "nephro", "creatinine", "urinary", "proteinuria"),
        "adr_term_keyword",
        0.75,
    ),
    SiteRule(
        "endocrine",
        "endocrine tissue",
        ("endocrine", "thyroid", "glucose", "insulin", "adrenal", "hormone", "steroid", "prolactin", "reproductive", "fertility"),
        "adr_term_keyword",
        0.70,
    ),
    SiteRule(
        "immune_blood",
        "immune and blood",
        (
            "immune",
            "hypersensitivity",
            "allergic",
            "anaphyl",
            "blood",
            "neutrop",
            "thrombocyt",
            "leukopen",
            "anaemia",
            "anemia",
            "lymphocyte",
        ),
        "adr_term_keyword",
        0.70,
    ),
    SiteRule(
        "lung",
        "lung",
        ("lung", "pulmonary", "respiratory", "bronch", "dyspnea", "dyspnoea", "pneumon"),
        "adr_term_keyword",
        0.70,
    ),
    SiteRule(
        "skin",
        "skin",
        ("skin", "rash", "dermat", "alopecia", "urticaria", "pruritus", "photosensitivity"),
        "adr_term_keyword",
        0.65,
    ),
)

SOC_SITE_HINTS: tuple[SiteRule, ...] = (
    SiteRule("heart", "heart", ("cardiac disorders",), "meddra_soc_keyword", 0.90),
    SiteRule("liver", "liver", ("hepatobiliary disorders",), "meddra_soc_keyword", 0.90),
    SiteRule("brain_cns", "brain", ("nervous system disorders", "psychiatric disorders"), "meddra_soc_keyword", 0.85),
    SiteRule("gi", "gastrointestinal tract", ("gastrointestinal disorders",), "meddra_soc_keyword", 0.90),
    SiteRule("kidney", "kidney", ("renal and urinary disorders",), "meddra_soc_keyword", 0.90),
    SiteRule("endocrine", "endocrine tissue", ("endocrine disorders", "metabolism and nutrition disorders"), "meddra_soc_keyword", 0.85),
    SiteRule("immune_blood", "immune and blood", ("immune system disorders", "blood and lymphatic system disorders"), "meddra_soc_keyword", 0.85),
    SiteRule("lung", "lung", ("respiratory, thoracic and mediastinal disorders",), "meddra_soc_keyword", 0.85),
    SiteRule("skin", "skin", ("skin and subcutaneous tissue disorders",), "meddra_soc_keyword", 0.85),
)

SITE_TISSUE_ALIASES: dict[str, tuple[str, ...]] = {
    "heart": ("heart", "cardiac", "myocard", "ventricle", "atrium"),
    "liver": ("liver", "hepatic", "hepatocyte"),
    "brain_cns": ("brain", "cerebral", "cortex", "cerebell", "hippocampus", "amygdala", "spinal cord", "cns"),
    "gi": ("colon", "intestin", "stomach", "duodenum", "ileum", "rectum", "gastro", "esophagus", "oesophagus"),
    "kidney": ("kidney", "renal", "glomer", "tubule"),
    "endocrine": ("thyroid", "adrenal", "pituitary", "pancreas", "testis", "ovary", "endocrine"),
    "immune_blood": ("blood", "spleen", "lymph", "immune", "bone marrow", "leukocyte", "monocyte"),
    "lung": ("lung", "pulmonary", "bronch"),
    "skin": ("skin", "epiderm", "dermis"),
}


def _clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _lower(value: Any) -> str:
    return _clean(value).lower()


def _first_nonempty(row: pd.Series, cols: tuple[str, ...]) -> str:
    for col in cols:
        if col in row.index:
            text = _clean(row[col])
            if text:
                return text
    return ""


def _match_rule(text: str, rules: tuple[SiteRule, ...]) -> tuple[str, str, str, float, str] | None:
    if not text:
        return None
    lowered = text.lower()
    for rule in rules:
        matched = [term for term in rule.terms if term in lowered]
        if matched:
            return rule.site_group, rule.site_name, rule.source, rule.confidence, ";".join(matched[:8])
    return None


def infer_adr_site_mapping(row: pd.Series) -> dict[str, object]:
    existing = _first_nonempty(row, ("adr_site_group", "site_group"))
    if existing:
        site_group = existing.lower().replace(" ", "_")
        return {
            "adr_site_group": site_group,
            "site_name": _first_nonempty(row, ("site_name",)) or site_group.replace("_", " "),
            "adr_site_mapping_source": _first_nonempty(row, ("adr_site_mapping_source",)) or "existing_column",
            "adr_site_mapping_confidence": pd.to_numeric(
                pd.Series([row.get("adr_site_mapping_confidence", 0.95)]), errors="coerce"
            ).iloc[0],
            "adr_site_terms_matched": _first_nonempty(row, ("adr_site_terms_matched",)),
        }

    soc_text = " | ".join(_lower(row.get(col, "")) for col in ("adr_soc", "meddra_soc", "soc_name", "system_organ_class"))
    existing_mapping_source = _first_nonempty(row, ("adr_site_mapping_source",))
    soc_match = _match_rule(soc_text, SOC_SITE_HINTS)
    if soc_match is not None:
        site_group, site_name, source, confidence, terms = soc_match
        return {
            "adr_site_group": site_group,
            "site_name": _first_nonempty(row, ("site_name",)) or site_name,
            "adr_site_mapping_source": existing_mapping_source or source,
            "adr_site_mapping_confidence": confidence,
            "adr_site_terms_matched": terms,
        }

    term_text = " | ".join(
        _lower(row.get(col, ""))
        for col in ("adr_term", "adr_pt", "meddra_pt", "preferred_term", "side_effect_name", "toxicity_term")
    )
    term_match = _match_rule(term_text, ADR_SITE_RULES)
    if term_match is not None:
        site_group, site_name, source, confidence, terms = term_match
        return {
            "adr_site_group": site_group,
            "site_name": _first_nonempty(row, ("site_name",)) or site_name,
            "adr_site_mapping_source": existing_mapping_source or source,
            "adr_site_mapping_confidence": confidence,
            "adr_site_terms_matched": terms,
        }

    return {
        "adr_site_group": pd.NA,
        "site_name": _first_nonempty(row, ("site_name",)) or pd.NA,
        "adr_site_mapping_source": "unmapped",
        "adr_site_mapping_confidence": pd.NA,
        "adr_site_terms_matched": "",
    }


def add_adr_site_mapping(df: pd.DataFrame) -> pd.DataFrame:
    """Infer broad ADR-site context from existing ADR fields.

    This creates context/provenance columns only. It does not create supervised
    labels and must not be interpreted as evidence that an ADR is caused by a
    target.
    """

    if df.empty:
        return df.copy()
    mapped = df.apply(infer_adr_site_mapping, axis=1, result_type="expand")
    out = df.copy()
    for col in mapped.columns:
        if col in out.columns:
            out[col] = out[col].where(out[col].notna() & out[col].astype(str).str.strip().ne(""), mapped[col])
        else:
            out[col] = mapped[col]
    return out


def adr_site_match_score(row: pd.Series) -> float | pd.NA:
    site_group = _lower(row.get("adr_site_group", ""))
    if not site_group:
        return pd.NA
    aliases = SITE_TISSUE_ALIASES.get(site_group, (site_group.replace("_", " "),))
    tissue_text = " | ".join(
        _lower(row.get(col, ""))
        for col in (
            "hpa_max_tissue",
            "gtex_max_tissue",
            "bgee_max_tissue",
            "ot_expression_max_tissue",
            "ot_safety_biosample_tissues",
            "site_name",
        )
    )
    if not tissue_text:
        return pd.NA
    if any(alias in tissue_text for alias in aliases):
        return 1.0
    if site_group in {"brain_cns", "immune_blood"} and any(token in tissue_text for token in site_group.split("_")):
        return 0.7
    return 0.0

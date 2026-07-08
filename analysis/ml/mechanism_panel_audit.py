from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pandas as pd


PANEL_SPECS: dict[str, dict[str, object]] = {
    "cardiac_qt": {
        "site_groups": ["heart"],
        "terms": [
            "qt", "torsade", "arrhythm", "ventricular", "cardiac", "cardiotox",
            "herg", "kcnh2", "long qt", "repolarization",
        ],
    },
    "cns_sedation": {
        "site_groups": ["brain_cns"],
        "terms": [
            "sedation", "somnolence", "drows", "cns", "central nervous", "seizure",
            "extrapyramidal", "parkinson", "hallucination", "psychosis", "brain",
        ],
    },
    "liver": {
        "site_groups": ["liver"],
        "terms": ["liver", "hepatic", "hepatotox", "cholestasis", "bile", "bilirubin", "transaminase"],
    },
    "gi": {
        "site_groups": ["gi"],
        "terms": ["gastro", "intestinal", "stomach", "ulcer", "bleed", "nausea", "vomit", "diarrhea"],
    },
    "kidney": {
        "site_groups": ["kidney"],
        "terms": ["kidney", "renal", "nephro", "creatinine"],
    },
    "endocrine": {
        "site_groups": ["endocrine"],
        "terms": ["endocrine", "thyroid", "adrenal", "hormone", "glucose", "insulin", "prolactin"],
    },
    "immune_blood": {
        "site_groups": ["immune_blood"],
        "terms": ["immune", "blood", "neutro", "platelet", "anemia", "leukopenia", "thrombocyt"],
    },
}

TEXT_COLS = [
    "adr_term",
    "adr_id",
    "meddra_pt",
    "target_adr_terms_aggregated",
    "mechanism_pu_mechanism_label_source_all",
]
SOURCE_COLS = ["label_source", "mechanism_pu_label_source_all", "source_family", "upstream_source"]


def _tokens(value: Any) -> list[str]:
    if pd.isna(value):
        return []
    out: list[str] = []
    for token in str(value).split(";"):
        cleaned = " ".join(token.strip().split())
        if cleaned:
            out.append(cleaned)
    return out


def _combined_text(df: pd.DataFrame) -> pd.Series:
    text = pd.Series("", index=df.index, dtype="object")
    for col in TEXT_COLS:
        if col not in df.columns:
            continue
        text = text.str.cat(df[col].fillna("").astype(str), sep=" ")
    return text.str.lower()


def _source_count(frame: pd.DataFrame) -> int:
    values: set[str] = set()
    for col in SOURCE_COLS:
        if col not in frame.columns:
            continue
        for value in frame[col].dropna().astype(str):
            values.update(_tokens(value))
    return len(values)


def panel_component_masks(df: pd.DataFrame, panel: str) -> tuple[pd.Series, pd.Series]:
    if panel == "all":
        all_rows = pd.Series(True, index=df.index)
        return all_rows, pd.Series(False, index=df.index)
    spec = PANEL_SPECS.get(panel)
    if spec is None:
        raise ValueError(f"unknown ADR panel: {panel}")
    site_groups = {str(value).lower() for value in spec.get("site_groups", [])}  # type: ignore[union-attr]
    site_mask = pd.Series(False, index=df.index)
    if "adr_site_group" in df.columns and site_groups:
        site_mask = df["adr_site_group"].fillna("").astype(str).str.lower().isin(site_groups)
    terms = [str(term).lower() for term in spec.get("terms", [])]  # type: ignore[union-attr]
    term_mask = pd.Series(False, index=df.index)
    if terms:
        pattern = re.compile("|".join(re.escape(term) for term in terms))
        term_mask = _combined_text(df).map(lambda text: bool(pattern.search(text)))
    return site_mask.fillna(False).astype(bool), term_mask.fillna(False).astype(bool)


def panel_mask(df: pd.DataFrame, panel: str) -> pd.Series:
    site_mask, term_mask = panel_component_masks(df, panel)
    return (site_mask | term_mask).fillna(False).astype(bool)


def filter_mechanism_panel(df: pd.DataFrame, panel: str) -> pd.DataFrame:
    return df.loc[panel_mask(df, panel)].copy()


def _label_counts(frame: pd.DataFrame, label_col: str) -> dict[str, int]:
    labels = pd.to_numeric(frame.get(label_col, pd.Series(pd.NA, index=frame.index)), errors="coerce")
    return {
        "rows": int(len(frame)),
        "labelable_rows": int(labels.notna().sum()),
        "positive_rows": int(labels.eq(1).sum()),
        "negative_rows": int(labels.eq(0).sum()),
        "unknown_rows": int(labels.isna().sum()),
    }


def _unique_count(frame: pd.DataFrame, col: str) -> int:
    return int(frame[col].dropna().astype(str).nunique()) if col in frame.columns else 0


def _source_holdout_ready(frame: pd.DataFrame, label_col: str, min_pos: int, min_neg: int) -> int:
    source_col = "label_source" if "label_source" in frame.columns else None
    if source_col is None:
        return 0
    ready = 0
    labels = pd.to_numeric(frame.get(label_col, pd.Series(pd.NA, index=frame.index)), errors="coerce")
    for source in sorted({token for value in frame[source_col].dropna().astype(str) for token in _tokens(value)}):
        source_mask = frame[source_col].fillna("").astype(str).map(lambda value: source in _tokens(value))
        source_labels = labels.loc[source_mask]
        if int(source_labels.eq(1).sum()) >= min_pos and int(source_labels.eq(0).sum()) >= min_neg:
            ready += 1
    return ready


def _recommend(row: dict[str, Any]) -> str:
    pos = int(row["positive_rows"])
    neg = int(row["negative_rows"])
    labelable = int(row["labelable_rows"])
    targets = int(row["unique_targets"])
    ready_sources = int(row["source_holdout_ready_sources"])
    site_pos = int(row.get("positive_site_group_rows", 0))
    site_neg = int(row.get("negative_site_group_rows", 0))
    term_only_pos = int(row.get("positive_term_only_rows", 0))
    if pos >= 50 and neg >= 50 and targets >= 5 and ready_sources >= 1 and site_pos >= 10 and site_neg >= 10:
        return "best_current_panel_candidate"
    if pos >= 20 and neg >= 20 and targets >= 3 and site_pos > 0:
        return "trainable_sensitivity_panel"
    if pos >= 20 and neg >= 20 and term_only_pos >= pos * 0.5:
        return "keyword_matched_trainable_needs_manual_review"
    if pos >= 5 and neg >= 5:
        return "small_exploratory_panel"
    if pos > 0:
        return "positive_only_or_negative_sparse_evidence_panel"
    if neg > 0:
        return "negative_only_panel_not_trainable"
    if labelable == 0:
        return "evidence_context_only_no_supervised_labels"
    return "not_ready"


def _readiness_score(row: dict[str, Any]) -> float:
    pos = min(int(row["positive_rows"]) / 50.0, 1.0)
    neg = min(int(row["negative_rows"]) / 50.0, 1.0)
    targets = min(int(row["unique_targets"]) / 10.0, 1.0)
    drugs = min(int(row["unique_drugs"]) / 50.0, 1.0)
    sources = min(int(row["unique_sources"]) / 3.0, 1.0)
    heldout = min(int(row["source_holdout_ready_sources"]) / 2.0, 1.0)
    return round(0.30 * pos + 0.25 * neg + 0.15 * targets + 0.10 * drugs + 0.10 * sources + 0.10 * heldout, 4)


def panel_readiness_audit(
    df: pd.DataFrame,
    *,
    label_col: str = "mechanism_pu_label",
    min_source_pos: int = 5,
    min_source_neg: int = 5,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, Any]] = []
    source_rows: list[dict[str, Any]] = []
    panels = list(PANEL_SPECS)
    for panel in panels:
        site_mask, term_mask = panel_component_masks(df, panel)
        mask = site_mask | term_mask
        frame = df.loc[mask].copy()
        site_frame = df.loc[site_mask].copy()
        term_only_frame = df.loc[term_mask & ~site_mask].copy()
        counts = _label_counts(frame, label_col)
        site_counts = _label_counts(site_frame, label_col)
        term_only_counts = _label_counts(term_only_frame, label_col)
        labels = pd.to_numeric(frame.get(label_col, pd.Series(pd.NA, index=frame.index)), errors="coerce")
        row: dict[str, Any] = {
            "panel": panel,
            **counts,
            "site_group_match_rows": site_counts["rows"],
            "labelable_site_group_rows": site_counts["labelable_rows"],
            "positive_site_group_rows": site_counts["positive_rows"],
            "negative_site_group_rows": site_counts["negative_rows"],
            "term_only_match_rows": term_only_counts["rows"],
            "labelable_term_only_rows": term_only_counts["labelable_rows"],
            "positive_term_only_rows": term_only_counts["positive_rows"],
            "negative_term_only_rows": term_only_counts["negative_rows"],
            "positive_rate_labelable": float(labels.mean()) if labels.notna().any() else None,
            "unique_drugs": _unique_count(frame, "drug_id"),
            "unique_targets": _unique_count(frame, "target_id"),
            "unique_pdbs": _unique_count(frame, "pdb_id"),
            "unique_sources": _source_count(frame),
            "source_holdout_ready_sources": _source_holdout_ready(frame, label_col, min_source_pos, min_source_neg),
        }
        row["readiness_score"] = _readiness_score(row)
        row["recommendation"] = _recommend(row)
        rows.append(row)
        if "label_source" in frame.columns:
            for source, group in frame.groupby("label_source", dropna=False):
                source_counts = _label_counts(group, label_col)
                source_rows.append({"panel": panel, "label_source": source, **source_counts})
    readiness = pd.DataFrame(rows).sort_values(
        ["readiness_score", "positive_rows", "negative_rows"], ascending=[False, False, False]
    )
    source_balance = pd.DataFrame(source_rows)
    return readiness, source_balance


def write_panel_readiness_audit(
    table_path: str | Path,
    out_path: str | Path,
    *,
    label_col: str = "mechanism_pu_label",
    source_balance_path: str | Path | None = None,
) -> pd.DataFrame:
    df = pd.read_csv(table_path, low_memory=False)
    readiness, source_balance = panel_readiness_audit(df, label_col=label_col)
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    readiness.to_csv(out, index=False)
    if source_balance_path is not None:
        source_out = Path(source_balance_path)
        source_out.parent.mkdir(parents=True, exist_ok=True)
        source_balance.to_csv(source_out, index=False)
    manifest = {
        "source_table": str(table_path),
        "out_path": str(out),
        "source_balance_path": str(source_balance_path) if source_balance_path else None,
        "label_col": label_col,
        "panels": list(PANEL_SPECS),
        "policy": "ADR panel audit filters existing mechanism labels by ADR/site context; it does not create new labels or change model architecture.",
    }
    out.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return readiness

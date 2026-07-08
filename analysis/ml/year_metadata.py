from __future__ import annotations

from typing import Any

import pandas as pd


YEAR_PRIORITY = [
    ("source_available_date", "source_available_date", "high"),
    ("source_release_date", "source_release_date", "high"),
    ("activity_publication_year", "activity_publication_year", "high"),
    ("document_year", "document_year", "medium"),
    ("evidence_publication_year", "evidence_publication_year", "medium"),
    ("database_release_year", "database_release_year", "low"),
]


def _year_from_value(value: Any) -> float | None:
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    if not text:
        return None
    parsed = pd.to_datetime(text, errors="coerce")
    if pd.notna(parsed):
        return float(parsed.year)
    numeric = pd.to_numeric(pd.Series([text]), errors="coerce").iloc[0]
    if pd.notna(numeric) and 1800 <= float(numeric) <= 2200:
        return float(int(numeric))
    return None


def add_availability_year_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """Add conservative row-level availability year fields from source metadata.

    This does not infer assay years from a database dump unless no better source is
    available. The confidence/source columns are intended to keep temporal claims
    honest in downstream audit reports.
    """
    out = df.copy()
    availability: list[float | None] = []
    basis: list[str | None] = []
    confidence: list[str | None] = []
    for _idx, row in out.iterrows():
        chosen_year = None
        chosen_basis = None
        chosen_confidence = None
        for col, source_name, conf in YEAR_PRIORITY:
            if col not in out.columns:
                continue
            year = _year_from_value(row.get(col))
            if year is None:
                continue
            chosen_year = year
            chosen_basis = source_name
            chosen_confidence = conf
            break
        availability.append(chosen_year)
        basis.append(chosen_basis)
        confidence.append(chosen_confidence)
    out["availability_year"] = pd.Series(availability, index=out.index)
    out["availability_year_source"] = pd.Series(basis, index=out.index)
    out["availability_year_confidence"] = pd.Series(confidence, index=out.index)
    return out

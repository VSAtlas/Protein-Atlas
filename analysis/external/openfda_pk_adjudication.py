from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd


@dataclass(frozen=True)
class CuratedCmaxContext:
    decision: str
    dose_value: float | None
    dose_unit: str
    route: str
    regimen: str
    context: str


def _normalize_unit(value: Any) -> str:
    unit = str(value or "").strip().replace("μ", "u").replace("µ", "u").casefold()
    aliases = {
        "ug/ml": "mcg/ml",
        "mcg/ml": "mcg/ml",
        "ng/ml": "ng/ml",
        "nm": "nm",
        "um": "um",
        "nmol/l": "nmol/l",
    }
    return aliases.get(unit, unit)


def _key(drug: str, value: float, unit: str) -> tuple[str, float, str]:
    return (drug.casefold(), float(value), _normalize_unit(unit))


def _context(
    decision: str,
    dose_value: float | None,
    dose_unit: str,
    route: str,
    regimen: str,
    context: str,
) -> CuratedCmaxContext:
    return CuratedCmaxContext(
        decision=decision,
        dose_value=dose_value,
        dose_unit=dose_unit,
        route=route,
        regimen=regimen,
        context=context,
    )


# These entries were checked against the cached SPL source excerpts. They are
# deliberately context-specific: none is a universal per-drug Cmax.
_CURATED_CONTEXTS = {
    _key("alectinib", 665, "ng/mL"): _context(
        "accept_context_only", None, "", "oral", "steady state", "ALK-positive NSCLC patients"
    ),
    _key("azacitidine", 10.6, "uM"): _context(
        "accept_context_only", None, "", "intravenous", "clinical Cmax", "IV azacitidine"
    ),
    _key("azacitidine", 145, "ng/mL"): _context(
        "accept_model_context", 300, "mg", "oral", "single dose", "ONUREG"
    ),
    _key("bromocriptine", 0.82, "nM"): _context(
        "accept_context_only", None, "", "oral", "therapeutic use", "patients; dose not bound to value"
    ),
    _key("budesonide", 3.3, "nmol/L"): _context(
        "accept_model_context", 960, "mcg", "inhalation", "single dose", "COPD patients"
    ),
    _key("carbinoxamine", 24, "ng/mL"): _context(
        "accept_model_context", 8, "mg", "oral", "single dose", "healthy volunteers"
    ),
    _key("clindamycin", 69.2, "ng/mL"): _context(
        "accept_model_context", 100, "mg", "vaginal", "single dose", "healthy volunteers"
    ),
    _key("clozapine", 319, "ng/mL"): _context(
        "accept_model_context", 100, "mg", "oral", "twice daily steady state", "tablet"
    ),
    _key("clozapine", 413, "ng/mL"): _context(
        "accept_model_context", 100, "mg", "oral", "twice daily steady state", "ODT"
    ),
    _key("clozapine", 275, "ng/mL"): _context(
        "accept_context_only", None, "", "oral", "steady state", "oral suspension; dose range reported"
    ),
    _key("dacomitinib", 108, "ng/mL"): _context(
        "accept_model_context", 45, "mg", "oral", "once daily steady state", "solid-tumor patients"
    ),
    _key("dasatinib", 82.2, "ng/mL"): _context(
        "accept_model_context", 100, "mg", "oral", "once daily steady state", "adult subjects"
    ),
    _key("enasidenib", 1.4, "mcg/mL"): _context(
        "accept_model_context", 100, "mg", "oral", "single dose", "parent drug"
    ),
    _key("lofexidine", 0.82, "ng/mL"): _context(
        "accept_model_context", 0.36, "mg", "oral", "single dose", "solution"
    ),
    _key("mercaptopurine", 69, "ng/mL"): _context(
        "accept_model_context", 50, "mg", "oral", "single fasted dose", "adult healthy subjects"
    ),
    _key("mercaptopurine", 93, "ng/mL"): _context(
        "accept_model_context", 50, "mg", "oral", "single fasted dose", "adult healthy subjects"
    ),
    _key("nabilone", 2, "ng/mL"): _context(
        "accept_model_context", 2, "mg", "oral", "single dose", "radiolabeled nabilone"
    ),
    _key("nimodipine", 69.9, "ng/mL"): _context(
        "accept_model_context", 60, "mg", "oral", "single dose", "NYMALIZE"
    ),
    _key("paroxetine", 61.7, "ng/mL"): _context(
        "accept_model_context", 30, "mg", "oral", "once daily steady state", "immediate-release tablet"
    ),
    _key("paroxetine", 30, "ng/mL"): _context(
        "accept_model_context", 25, "mg", "oral", "once daily steady state", "controlled release"
    ),
    _key("pitolisant", 73, "ng/mL"): _context(
        "accept_model_context", 35.6, "mg", "oral", "once daily steady state", "WAKIX"
    ),
    _key("pravastatin", 26.5, "ng/mL"): _context(
        "accept_model_context", 20, "mg", "oral", "single fasted dose", "parent drug"
    ),
    _key("ranolazine", 2600, "ng/mL"): _context(
        "accept_model_context", 1000, "mg", "oral", "twice daily steady state", "extended release"
    ),
    _key("saxagliptin", 24, "ng/mL"): _context(
        "accept_model_context", 5, "mg", "oral", "single dose", "healthy subjects; parent drug"
    ),
    _key("selegiline", 1, "ng/mL"): _context(
        "accept_model_context", 10, "mg", "oral", "single dose", "parent drug"
    ),
    _key("sirolimus", 2590, "ng/mL"): _context(
        "accept_context_only", None, "", "intravenous", "recommended regimen", "FYARRO formulation"
    ),
    _key("sitagliptin", 950, "nM"): _context(
        "accept_model_context", 100, "mg", "oral", "single dose", "healthy volunteers"
    ),
    _key("sonidegib", 1030, "ng/mL"): _context(
        "accept_model_context", 200, "mg", "oral", "once daily steady state", "cancer patients"
    ),
    _key("temozolomide", 7.5, "mcg/mL"): _context(
        "accept_model_context", 150, "mg/m2", "oral", "single dose", "parent drug"
    ),
    _key("temsirolimus", 585, "ng/mL"): _context(
        "accept_model_context", 25, "mg", "intravenous", "single dose", "cancer patients; whole blood"
    ),
    _key("timolol", 0.46, "ng/mL"): _context(
        "accept_context_only", None, "", "ophthalmic", "0.5% twice daily", "systemic plasma exposure"
    ),
    _key("tolcapone", 3, "mcg/mL"): _context(
        "accept_model_context", 100, "mg", "oral", "three times daily", "steady-state approximation"
    ),
    _key("zanamivir", 43, "ng/mL"): _context(
        "accept_model_context", 10, "mg", "inhalation", "single dose", "children aged 6-12"
    ),
}

_KNOWN_EXTRACTION_FAILURES = {
    "abiraterone": "standard_deviation_extracted_instead_of_cmax",
    "aminocaproic acid": "standard_deviation_extracted_instead_of_peak_mean",
    "deoxycholic acid": "standard_deviation_extracted_instead_of_cmax",
    "dorzolamide": "combination_product_other_analyte_timolol",
    "lomitapide": "auc_value_extracted_instead_of_cmax",
    "repaglinide": "auc_values_extracted_instead_of_cmax",
    "rimantadine": "standard_deviation_extracted_instead_of_peak_mean",
    "tacrolimus": "standard_deviation_extracted_instead_of_cmax",
}


def _single_candidate(row: pd.Series) -> tuple[float, str] | None:
    values = [part for part in str(row.get("cmax_values_raw") or "").split(";") if part]
    units = [part for part in str(row.get("cmax_units_raw") or "").split(";") if part]
    if len(values) != 1 or len(set(map(_normalize_unit, units))) != 1:
        return None
    try:
        return float(values[0]), _normalize_unit(units[0])
    except ValueError:
        return None


def adjudicate_openfda_cmax(review: pd.DataFrame) -> pd.DataFrame:
    """Apply source-text adjudication to every cached SPL Cmax candidate row."""

    decisions: list[dict[str, Any]] = []
    for _, row in review.iterrows():
        drug = str(row.get("drug_id") or "").strip().casefold()
        candidate = _single_candidate(row)
        curated = _CURATED_CONTEXTS.get((drug, *candidate)) if candidate else None
        failure = _KNOWN_EXTRACTION_FAILURES.get(drug)
        if int(row.get("cmax_candidate_count") or 0) <= 0:
            decision = "not_applicable_no_cmax"
            reason = "no_cmax_candidate"
        elif failure:
            decision = "reject_extraction"
            reason = failure
        elif curated:
            decision = curated.decision
            reason = "source_text_value_and_analyte_verified"
        elif candidate is None:
            decision = "reject_current_row_shape"
            reason = "multiple_values_require_context_splitting"
        else:
            decision = "reject_unresolved_context"
            reason = "value_not_safely_bound_to_analyte_dose_route_regimen"

        decisions.append(
            {
                "adjudication_status": decision,
                "adjudication_reason": reason,
                "adjudicated_dose_value": curated.dose_value if curated else pd.NA,
                "adjudicated_dose_unit": curated.dose_unit if curated else "",
                "adjudicated_route": curated.route if curated else "",
                "adjudicated_regimen": curated.regimen if curated else "",
                "adjudicated_context": curated.context if curated else "",
                "acceptable_for_contextual_pk": bool(curated),
                "acceptable_for_model_training": bool(
                    curated and curated.decision == "accept_model_context"
                ),
                "acceptable_as_universal_drug_cmax": False,
                "manual_review_version": "atlas_spl_cmax_adjudication_v1",
            }
        )
    return pd.concat(
        [review.reset_index(drop=True), pd.DataFrame(decisions)],
        axis=1,
    )


def cmax_adjudication_summary(adjudicated: pd.DataFrame) -> pd.DataFrame:
    candidates = adjudicated.loc[
        pd.to_numeric(
            adjudicated.get("cmax_candidate_count"),
            errors="coerce",
        ).fillna(0).gt(0)
    ]
    return (
        candidates.groupby("adjudication_status", dropna=False)
        .agg(
            records=("source_record_id", "size"),
            drugs=("drug_id", "nunique"),
            unique_contexts=(
                "adjudicated_context",
                lambda values: sum(bool(str(value)) for value in set(values)),
            ),
        )
        .reset_index()
    )

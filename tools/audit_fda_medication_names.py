"""Local audit for FDA mapping display names.

This intentionally does not query the network.  It checks whether the mapped
display names are populated, non-identifier-like, and worth human review as
non-medication FDA substances such as color additives.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Sequence

from chemdb.tools.fill_fda_mapping_names import row_display_name_quality_reasons
from tools.fda_mapping_index import extract_rdk_id


DEFAULT_MAPPING_CSV = Path("chemdb/data/fda_mapping_from_pdbqt.csv")
DEFAULT_CLASSIFIED_NON_MEDICATION_CSV = Path("docs/fda_non_medication_substances.csv")
DEFAULT_DETAIL_CSV = Path("docs/fda_medication_name_audit.csv")
DEFAULT_UNIQUE_CSV = Path("docs/fda_medication_name_review_unique.csv")
DEFAULT_SUMMARY_JSON = Path("docs/fda_medication_name_audit_summary.json")

NON_MEDICATION_PATTERNS: tuple[tuple[str, re.Pattern[str]], ...] = (
    ("color_additive_name", re.compile(r"\b(?:fd&c|d&c|d and c|lake|tartrazine)\b", re.I)),
    (
        "color_additive_name",
        re.compile(
            r"\b(?:carminic acid|carmine|chocolate brown|brilliant blue|allura red|"
            r"sunset yellow|erythrosine|indigo carmine|fast green|quinoline yellow|"
            r"ponceau|amaranth|annatto)\b",
            re.I,
        ),
    ),
    (
        "excipient_like_name",
        re.compile(
            r"\b(?:magnesium stearate|cellulose|microcrystalline cellulose|"
            r"lactose|starch|sucrose|sorbitol|xanthan gum|povidone|crospovidone|"
            r"polyethylene glycol|propylene glycol|polysorbate|hypromellose|"
            r"talc|silicon dioxide)\b",
            re.I,
        ),
    ),
    (
        "inorganic_excipient_like_name",
        re.compile(r"\b(?:titanium dioxide|iron oxide|ferric oxide)\b", re.I),
    ),
)

DRUG_NAME_SUFFIX_RE = re.compile(
    r"(?:"
    r"afil|azole|barbital|caine|cillin|cycline|dipine|floxacin|gliflozin|"
    r"gliptin|glitazone|lukast|mab|mycin|navir|olol|olone|pam|parin|"
    r"platin|prazole|pril|sartan|setron|statin|terol|thiazide|triptan|"
    r"vir|xaban|zepam|zolam"
    r")$",
    re.I,
)

THERAPEUTIC_WORD_RE = re.compile(
    r"\b(?:"
    r"acetate|alafenamide|benzoate|bromide|carbonate|carbamate|chloride|"
    r"citrate|dipropionate|fumarate|hydrochloride|maleate|mesylate|nitrate|"
    r"phosphate|propionate|sodium|succinate|sulfate|tartrate|valerate"
    r")\b",
    re.I,
)


def _clean(value: object) -> str:
    return str(value or "").strip()


def _norm(value: object) -> str:
    text = _clean(value).lower()
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def _split_names(value: object) -> list[str]:
    return [part.strip() for part in re.split(r"[|;]", _clean(value)) if part.strip()]


def _is_same_or_contained(left: str, right: str) -> bool:
    left_key = _norm(left)
    right_key = _norm(right)
    if not left_key or not right_key:
        return False
    return left_key == right_key or left_key in right_key or right_key in left_key


def _non_medication_reasons(name: str) -> list[str]:
    reasons: list[str] = []
    for reason, pattern in NON_MEDICATION_PATTERNS:
        if pattern.search(name):
            reasons.append(reason)
    return list(dict.fromkeys(reasons))


def _has_medication_name_shape(name: str) -> bool:
    words = [word for word in re.split(r"[^A-Za-z]+", name) if word]
    if not words:
        return False
    return any(DRUG_NAME_SUFFIX_RE.search(word) for word in words) or bool(
        THERAPEUTIC_WORD_RE.search(name)
    )


def _trusted_source_reasons(row: dict[str, str], display_name: str) -> list[str]:
    reasons: list[str] = []
    for field in ("rxnorm_generic_name", "drugcentral_generic_name"):
        value = _clean(row.get(field))
        if value and _is_same_or_contained(display_name, value):
            reasons.append(f"{field}_supports_display")
    for field in ("rxnorm_brand_names", "drugcentral_brand_names", "brand_names"):
        for value in _split_names(row.get(field)):
            if _is_same_or_contained(display_name, value):
                reasons.append(f"{field}_supports_display")
                break
    return reasons


def _conflict_reasons(row: dict[str, str], display_name: str) -> list[str]:
    reasons: list[str] = []
    display_key = _norm(display_name)
    if not display_key:
        return reasons
    for field in ("rxnorm_generic_name", "drugcentral_generic_name"):
        value = _clean(row.get(field))
        value_key = _norm(value)
        if value_key and not _is_same_or_contained(display_name, value):
            reasons.append(f"{field}_differs")
    return reasons


def _row_audit(
    row: dict[str, str],
    *,
    include_cross_field_conflicts: bool,
    include_weak_signals: bool,
) -> dict[str, str]:
    display_name = _clean(row.get("display_name"))
    quality_reasons = row_display_name_quality_reasons(row, strict=True)
    non_med_reasons = _non_medication_reasons(display_name)
    support_reasons = _trusted_source_reasons(row, display_name)
    conflict_reasons = _conflict_reasons(row, display_name)

    review_reasons: list[str] = []
    blocking_reasons = [f"name_quality:{reason}" for reason in quality_reasons]
    review_reasons.extend(f"non_medication:{reason}" for reason in non_med_reasons)
    if include_cross_field_conflicts:
        review_reasons.extend(f"cross_field:{reason}" for reason in conflict_reasons)
    if (
        include_weak_signals
        and display_name
        and not support_reasons
        and not _has_medication_name_shape(display_name)
    ):
        review_reasons.append("weak_local_medication_signal")

    if blocking_reasons:
        severity = "blocking"
    elif review_reasons:
        severity = "review"
    else:
        severity = "ok"

    return {
        "severity": severity,
        "rdk_id": extract_rdk_id(row.get("path")) or "",
        "file_num": _clean(row.get("file_num")),
        "display_name": display_name,
        "reason_codes": ";".join(blocking_reasons + review_reasons),
        "support_codes": ";".join(support_reasons),
        "pubchem_record_title": _clean(row.get("pubchem_record_title")),
        "generic_name": _clean(row.get("generic_name")),
        "rxnorm_generic_name": _clean(row.get("rxnorm_generic_name")),
        "drugcentral_generic_name": _clean(row.get("drugcentral_generic_name")),
        "sdf_title": _clean(row.get("sdf_title")),
        "inchikey": _clean(row.get("inchikey")),
        "path": _clean(row.get("path")),
    }


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _load_classified_non_medication(path: Path | None) -> dict[str, dict[str, str]]:
    if not path or not path.exists():
        return {}
    rows = _read_rows(path)
    classified: dict[str, dict[str, str]] = {}
    for row in rows:
        rdk_id = _clean(row.get("rdk_id"))
        display_name = _clean(row.get("display_name"))
        if not rdk_id or not display_name:
            raise ValueError(f"incomplete non-medication classification row: {row}")
        classified[rdk_id] = {
            "display_name": display_name,
            "substance_class": _clean(row.get("substance_class")),
            "status": _clean(row.get("status")),
            "note": _clean(row.get("note")),
        }
    return classified


def _write_csv(path: Path, rows: Sequence[dict[str, str]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _unique_review_rows(rows: Iterable[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row["severity"] not in {"blocking", "review"}:
            continue
        grouped[_norm(row["display_name"]) or row["rdk_id"]].append(row)

    unique_rows: list[dict[str, str]] = []
    for grouped_rows in grouped.values():
        first = grouped_rows[0]
        reason_counter = Counter(
            reason
            for row in grouped_rows
            for reason in row["reason_codes"].split(";")
            if reason
        )
        unique_rows.append(
            {
                "severity": "blocking"
                if any(row["severity"] == "blocking" for row in grouped_rows)
                else "review",
                "display_name": first["display_name"],
                "row_count": str(len(grouped_rows)),
                "example_rdk_id": first["rdk_id"],
                "reason_codes": ";".join(sorted(reason_counter)),
                "example_pubchem_record_title": first["pubchem_record_title"],
                "example_rxnorm_generic_name": first["rxnorm_generic_name"],
                "example_drugcentral_generic_name": first["drugcentral_generic_name"],
            }
        )
    unique_rows.sort(
        key=lambda row: (
            row["severity"] != "blocking",
            "non_medication:" not in row["reason_codes"],
            row["display_name"].lower(),
        )
    )
    return unique_rows


def run_audit(
    *,
    mapping_csv: Path,
    classified_non_medication_csv: Path | None,
    detail_csv: Path,
    unique_csv: Path,
    summary_json: Path,
    include_cross_field_conflicts: bool = False,
    include_weak_signals: bool = False,
) -> dict[str, object]:
    rows = _read_rows(mapping_csv)
    classified_non_medication = _load_classified_non_medication(
        classified_non_medication_csv
    )
    audited_rows = [
        _row_audit(
            row,
            include_cross_field_conflicts=include_cross_field_conflicts,
            include_weak_signals=include_weak_signals,
        )
        for row in rows
    ]
    classified_count = 0
    classified_unique_names: set[str] = set()
    for row in audited_rows:
        classification = classified_non_medication.get(row["rdk_id"])
        if not classification:
            continue
        expected_name = classification["display_name"]
        if _norm(row["display_name"]) != _norm(expected_name):
            row["reason_codes"] = ";".join(
                [
                    reason
                    for reason in [
                        row["reason_codes"],
                        "classification_mismatch:display_name",
                    ]
                    if reason
                ]
            )
            row["severity"] = "blocking"
            continue
        classified_count += 1
        classified_unique_names.add(_norm(row["display_name"]))
        if row["severity"] == "review" and all(
            reason.startswith("non_medication:")
            for reason in row["reason_codes"].split(";")
            if reason
        ):
            row["severity"] = "classified_non_medication"
    review_rows = [
        row for row in audited_rows if row["severity"] in {"blocking", "review"}
    ]
    unique_rows = _unique_review_rows(audited_rows)

    fieldnames = [
        "severity",
        "rdk_id",
        "file_num",
        "display_name",
        "reason_codes",
        "support_codes",
        "pubchem_record_title",
        "generic_name",
        "rxnorm_generic_name",
        "drugcentral_generic_name",
        "sdf_title",
        "inchikey",
        "path",
    ]
    _write_csv(detail_csv, review_rows, fieldnames)
    _write_csv(
        unique_csv,
        unique_rows,
        [
            "severity",
            "display_name",
            "row_count",
            "example_rdk_id",
            "reason_codes",
            "example_pubchem_record_title",
            "example_rxnorm_generic_name",
            "example_drugcentral_generic_name",
        ],
    )

    severity_counts = Counter(row["severity"] for row in audited_rows)
    reason_counts = Counter(
        reason
        for row in review_rows
        for reason in row["reason_codes"].split(";")
        if reason
    )
    unique_reason_counts = Counter(
        reason
        for row in unique_rows
        for reason in row["reason_codes"].split(";")
        if reason
    )
    summary: dict[str, object] = {
        "mapping_csv": str(mapping_csv),
        "classified_non_medication_csv": str(classified_non_medication_csv or ""),
        "total_rows": len(audited_rows),
        "classified_non_medication_rows": classified_count,
        "classified_non_medication_unique_names": len(classified_unique_names),
        "unique_display_names_reviewed": len(unique_rows),
        "row_severity_counts": dict(sorted(severity_counts.items())),
        "review_row_count": len(review_rows),
        "blocking_row_count": severity_counts.get("blocking", 0),
        "reason_counts": dict(sorted(reason_counts.items())),
        "unique_reason_counts": dict(sorted(unique_reason_counts.items())),
        "detail_csv": str(detail_csv),
        "unique_review_csv": str(unique_csv),
        "include_cross_field_conflicts": include_cross_field_conflicts,
        "include_weak_signals": include_weak_signals,
    }
    summary_json.parent.mkdir(parents=True, exist_ok=True)
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping-csv", type=Path, default=DEFAULT_MAPPING_CSV)
    parser.add_argument(
        "--classified-non-medication-csv",
        type=Path,
        default=DEFAULT_CLASSIFIED_NON_MEDICATION_CSV,
        help=(
            "Curated local table of FDA substances that are valid mappings but "
            "not medication generic names."
        ),
    )
    parser.add_argument("--detail-csv", type=Path, default=DEFAULT_DETAIL_CSV)
    parser.add_argument("--unique-csv", type=Path, default=DEFAULT_UNIQUE_CSV)
    parser.add_argument("--summary-json", type=Path, default=DEFAULT_SUMMARY_JSON)
    parser.add_argument(
        "--include-cross-field-conflicts",
        action="store_true",
        help=(
            "Also review rows where auxiliary RxNorm/DrugCentral fields disagree "
            "with display_name. This is noisy for the current generated CSV."
        ),
    )
    parser.add_argument(
        "--include-weak-signals",
        action="store_true",
        help=(
            "Also review names without local medication-shape/support signals. "
            "This is intentionally broad and can produce thousands of rows."
        ),
    )
    parser.add_argument(
        "--fail-on-blocking",
        action="store_true",
        help="Exit with code 2 if populated-name quality failures remain.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    summary = run_audit(
        mapping_csv=args.mapping_csv,
        classified_non_medication_csv=args.classified_non_medication_csv,
        detail_csv=args.detail_csv,
        unique_csv=args.unique_csv,
        summary_json=args.summary_json,
        include_cross_field_conflicts=bool(args.include_cross_field_conflicts),
        include_weak_signals=bool(args.include_weak_signals),
    )
    print(f"rows audited: {summary['total_rows']}")
    print(f"blocking rows: {summary['blocking_row_count']}")
    print(f"review rows: {summary['review_row_count']}")
    print(f"unique review names: {summary['unique_display_names_reviewed']}")
    print(f"detail_csv: {summary['detail_csv']}")
    print(f"unique_review_csv: {summary['unique_review_csv']}")
    print(f"summary_json: {args.summary_json}")
    blocking_row_count = int(str(summary["blocking_row_count"]))
    if args.fail_on_blocking and blocking_row_count > 0:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

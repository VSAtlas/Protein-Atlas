"""Apply source-backed display-name overrides for final FDA mapping gaps."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import re
import shutil
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from chemdb.tools.fill_fda_mapping_names import (
    display_name_quality_reasons,
    row_display_name_quality_reasons,
)
from tools.fda_mapping_index import MappingIndex, extract_rdk_id, iter_audited_rows


DEFAULT_MAPPING_CSV = Path("chemdb/data/fda_mapping_from_pdbqt.csv")
DEFAULT_OVERRIDES_CSV = Path("docs/fda_remaining_name_overrides.csv")


@dataclass(frozen=True)
class Override:
    rdk_id: str
    inchikey_prefix: str
    display_name: str
    source: str
    source_url: str
    note: str


def _clean(value: object) -> str:
    return str(value or "").strip()


def _rdk_id_from_row(row: dict[str, str]) -> str:
    return extract_rdk_id(row.get("path")) or ""


def _load_target_ids(path: Path | None) -> set[str]:
    if not path:
        return set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {_clean(row.get("rdk_id")) for row in reader if _clean(row.get("rdk_id"))}


def _load_overrides(path: Path) -> dict[str, Override]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        overrides: dict[str, Override] = {}
        for row in reader:
            override = Override(
                rdk_id=_clean(row.get("rdk_id")),
                inchikey_prefix=_clean(row.get("inchikey_prefix")).upper(),
                display_name=_clean(row.get("display_name")),
                source=_clean(row.get("source")),
                source_url=_clean(row.get("source_url")),
                note=_clean(row.get("note")),
            )
            if not override.rdk_id or not override.display_name:
                raise ValueError(f"incomplete override row: {row}")
            if override.rdk_id in overrides:
                raise ValueError(f"duplicate override for {override.rdk_id}")
            reasons = display_name_quality_reasons(override.display_name, strict=True)
            if reasons:
                raise ValueError(
                    f"override {override.rdk_id} has low-quality name "
                    f"{override.display_name!r}: {','.join(reasons)}"
                )
            overrides[override.rdk_id] = override
    return overrides


def _write_rows(
    path: Path,
    rows: Sequence[dict[str, str]],
    fieldnames: Sequence[str],
    *,
    delimiter: str = ",",
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def _audit_rows(mapping_csv: Path) -> list[dict[str, str]]:
    index = MappingIndex(mapping_csv)
    audit_rows: list[dict[str, str]] = []
    for row, reasons in iter_audited_rows(index.rows):
        rdk_id = extract_rdk_id(row.path) or ""
        audit_rows.append(
            {
                "rdk_id": rdk_id,
                "display_name": row.display_name,
                "suggested_name": row.preferred_name()
                if row.preferred_name() != row.display_name
                else "",
                "reason_codes": ";".join(reasons),
                "path": row.path,
            }
        )
    return audit_rows


def _write_name_audit(path: Path, mapping_csv: Path) -> int:
    audit_rows = _audit_rows(mapping_csv)
    _write_rows(
        path,
        audit_rows,
        ["rdk_id", "display_name", "suggested_name", "reason_codes", "path"],
        delimiter="\t",
    )
    return len(audit_rows)


def _write_final_audit(path: Path, summary_path: Path, mapping_csv: Path) -> None:
    index = MappingIndex(mapping_csv)
    rows: list[dict[str, str]] = []
    review_reasons: Counter[str] = Counter()
    blocking_reasons: Counter[str] = Counter()
    for row in index.rows:
        rdk_id = extract_rdk_id(row.path) or ""
        reasons = [f"row_quality_{reason}" for reason in row.audit_reasons()]
        raw_reasons = display_name_quality_reasons(
            row.pubchem_record_title,
            row.pubchem_iupac_name,
            strict=True,
        )
        reasons.extend(f"raw_quality_{reason}" for reason in raw_reasons)
        if (
            row.pubchem_record_title
            and row.display_name
            and row.pubchem_record_title.lower() != row.display_name.lower()
            and not display_name_quality_reasons(row.pubchem_record_title, strict=True)
        ):
            reasons.append("brand_display_with_clean_pubchem_title")
        reasons = list(dict.fromkeys(reason for reason in reasons if reason))
        if not reasons:
            continue
        severity = "blocking" if any(reason.startswith("row_quality_") for reason in reasons) else "review"
        for reason in reasons:
            review_reasons[reason] += 1
            if severity == "blocking":
                blocking_reasons[reason] += 1
        rows.append(
            {
                "severity": severity,
                "file_num": re.sub(r"\D+", "", rdk_id).lstrip("0") or "",
                "rdk_id": rdk_id,
                "display_name": row.display_name,
                "reason_codes": ";".join(reasons),
                "pubchem_record_title": row.pubchem_record_title,
                "rxnorm_generic_name": row.rxnorm_generic_name,
                "drugcentral_generic_name": row.drugcentral_generic_name,
                "sdf_title": row.sdf_title,
                "inchikey": row.inchikey,
                "path": row.path,
            }
        )
    rows.sort(key=lambda item: (item["severity"] != "blocking", int(item["file_num"] or 0)))
    _write_rows(
        path,
        rows,
        [
            "severity",
            "file_num",
            "rdk_id",
            "display_name",
            "reason_codes",
            "pubchem_record_title",
            "rxnorm_generic_name",
            "drugcentral_generic_name",
            "sdf_title",
            "inchikey",
            "path",
        ],
    )
    summary = {
        "mapping_csv": str(mapping_csv),
        "total_rows": len(index.rows),
        "blocking_rows": sum(1 for row in rows if row["severity"] == "blocking"),
        "review_rows": len(rows),
        "benchmark_audit_flagged_rows": sum(1 for row in rows if row["severity"] == "blocking"),
        "blocking_reason_counts": dict(sorted(blocking_reasons.items())),
        "review_reason_counts": dict(sorted(review_reasons.items())),
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_flagged_report(path: Path, mapping_csv: Path) -> int:
    rows, _fieldnames = _iter_rows(mapping_csv)
    flagged_rows: list[dict[str, str]] = []
    for row in rows:
        reasons = row_display_name_quality_reasons(row, strict=True)
        if not reasons:
            continue
        flagged_rows.append(
            {
                "sdf_title": _clean(row.get("sdf_title")),
                "inchikey": _clean(row.get("inchikey")),
                "current_display_name": _clean(row.get("display_name")),
                "reason_codes": ";".join(reasons),
                "pubchem_cid_resolved": _clean(row.get("pubchem_cid_resolved")),
                "pubchem_record_title": _clean(row.get("pubchem_record_title")),
                "pubchem_unii_list": _clean(row.get("pubchem_unii_list")),
            }
        )
    _write_rows(
        path,
        flagged_rows,
        [
            "sdf_title",
            "inchikey",
            "current_display_name",
            "reason_codes",
            "pubchem_cid_resolved",
            "pubchem_record_title",
            "pubchem_unii_list",
        ],
    )
    return len(flagged_rows)


def _iter_rows(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader), list(reader.fieldnames or [])


def apply_overrides(
    *,
    mapping_csv: Path,
    overrides_csv: Path,
    output_csv: Path,
    audit_tsv: Path | None,
    report_tsv: Path,
    unresolved_path: Path,
) -> int:
    overrides = _load_overrides(overrides_csv)
    rows, fieldnames = _iter_rows(mapping_csv)
    target_ids = _load_target_ids(audit_tsv) or set(overrides)
    report_rows: list[dict[str, str]] = []
    unresolved: list[str] = []
    updated = 0

    for row in rows:
        rdk_id = _rdk_id_from_row(row)
        if rdk_id not in target_ids:
            continue
        override = overrides.get(rdk_id)
        if override is None:
            unresolved.append(f"{rdk_id}\tmissing_override")
            continue
        inchikey = _clean(row.get("inchikey") or row.get("remark_inchikey")).upper()
        if override.inchikey_prefix and not inchikey.startswith(override.inchikey_prefix):
            unresolved.append(f"{rdk_id}\tinchikey_mismatch:{inchikey}")
            continue
        old_name = _clean(row.get("display_name"))
        row["display_name"] = override.display_name
        row["generic_name"] = override.display_name
        report_rows.append(
            {
                "rdk_id": rdk_id,
                "file_num": _clean(row.get("file_num")),
                "old_display_name": old_name,
                "new_display_name": override.display_name,
                "inchikey": inchikey,
                "source": override.source,
                "source_url": override.source_url,
                "note": override.note,
            }
        )
        if old_name != override.display_name:
            updated += 1

    if unresolved:
        unresolved_path.write_text("\n".join(unresolved) + "\n", encoding="utf-8")
    else:
        unresolved_path.write_text("", encoding="utf-8")
    _write_rows(output_csv, rows, fieldnames)
    _write_rows(
        report_tsv,
        report_rows,
        [
            "rdk_id",
            "file_num",
            "old_display_name",
            "new_display_name",
            "inchikey",
            "source",
            "source_url",
            "note",
        ],
    )
    return updated


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping-csv", type=Path, default=DEFAULT_MAPPING_CSV)
    parser.add_argument("--overrides-csv", type=Path, default=DEFAULT_OVERRIDES_CSV)
    parser.add_argument("--out-csv", type=Path)
    parser.add_argument("--inplace", action="store_true")
    parser.add_argument("--audit-tsv", type=Path, default=Path("docs/fda_mapping_name_audit.tsv"))
    parser.add_argument("--report-tsv", type=Path, default=Path("docs/fda_remaining_name_overrides.applied.tsv"))
    parser.add_argument("--unresolved-path", type=Path, default=Path("docs/fda_mapping_repair.unresolved.txt"))
    parser.add_argument("--refresh-audits", action="store_true")
    parser.add_argument("--name-audit-tsv", type=Path, default=Path("docs/fda_mapping_name_audit.tsv"))
    parser.add_argument(
        "--flagged-report-csv",
        type=Path,
        default=Path("docs/fda_mapping_repair.flagged_name_quality.csv"),
    )
    parser.add_argument("--final-audit-tsv", type=Path, default=Path("docs/fda_mapping_final_audit.tsv"))
    parser.add_argument(
        "--final-audit-summary-json",
        type=Path,
        default=Path("docs/fda_mapping_final_audit_summary.json"),
    )
    parser.add_argument("--fail-if-unresolved", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    mapping_csv = Path(args.mapping_csv)
    if args.inplace:
        output_csv = mapping_csv
        timestamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_path = mapping_csv.with_name(f"{mapping_csv.name}.{timestamp}.bak")
        shutil.copy2(mapping_csv, backup_path)
    elif args.out_csv:
        output_csv = Path(args.out_csv)
    else:
        parser.error("--out-csv is required unless --inplace is set")

    updated = apply_overrides(
        mapping_csv=mapping_csv,
        overrides_csv=Path(args.overrides_csv),
        output_csv=output_csv,
        audit_tsv=Path(args.audit_tsv) if args.audit_tsv else None,
        report_tsv=Path(args.report_tsv),
        unresolved_path=Path(args.unresolved_path),
    )
    unresolved_count = sum(
        1
        for line in Path(args.unresolved_path).read_text(encoding="utf-8").splitlines()
        if line.strip()
    )
    if args.refresh_audits:
        flagged = _write_name_audit(Path(args.name_audit_tsv), output_csv)
        quality_flagged = _write_flagged_report(Path(args.flagged_report_csv), output_csv)
        _write_final_audit(
            Path(args.final_audit_tsv),
            Path(args.final_audit_summary_json),
            output_csv,
        )
    else:
        flagged = -1
        quality_flagged = -1
    print(f"rows updated: {updated}")
    print(f"unresolved rows: {unresolved_count}")
    if flagged >= 0:
        print(f"name-audit flagged rows: {flagged}")
        print(f"quality-report flagged rows: {quality_flagged}")
    if args.fail_if_unresolved and unresolved_count:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

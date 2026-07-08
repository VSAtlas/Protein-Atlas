from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path
from typing import Any, Iterable

from analysis.reporting.ligand_side_effect_cache import (
    build_ligand_query,
    load_spd_exposure_rows,
    resolve_mapping_context,
    spd_exposure_for_query,
)

LOG = logging.getLogger("fill-master-exposure")


def _clean_text(value: Any) -> str:
    return " ".join(str(value or "").split()).strip()


def _is_decoy(row: dict[str, Any]) -> bool:
    return _clean_text(row.get("is_decoy")).lower() in {"1", "true", "yes", "y"}


def _has_exposure(row: dict[str, Any]) -> bool:
    return any(
        _clean_text(row.get(col))
        for col in (
            "free_cmax_um",
            "cmax_um",
            "fraction_unbound_plasma",
            "exposure_source",
        )
    )


def fill_master_exposure(
    *,
    repo_root: Path,
    input_csv: Path,
    output_csv: Path,
    spd_supplement_xlsx: Path,
    fda_mapping_csv: str = "",
    preserve_existing: bool = True,
) -> int:
    spd_rows = load_spd_exposure_rows(spd_supplement_xlsx)
    if not spd_rows:
        raise ValueError(f"no SPD exposure rows loaded from {spd_supplement_xlsx}")
    mapping_csv, fda_index, mapping_rows = resolve_mapping_context(
        repo_root, "", fda_mapping_csv
    )
    LOG.info("mapping_csv=%s spd_exposure_rows=%d", mapping_csv, len(spd_rows))

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    filled = 0
    rows = 0
    with input_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        exposure_cols = [
            "free_cmax_um",
            "cmax_um",
            "fraction_unbound_plasma",
            "exposure_source",
        ]
        for col in exposure_cols:
            if col not in fieldnames:
                fieldnames.append(col)
        with output_csv.open("w", encoding="utf-8", newline="") as out_handle:
            writer = csv.DictWriter(out_handle, fieldnames=fieldnames)
            writer.writeheader()
            for row in reader:
                rows += 1
                if not _is_decoy(row) and (not preserve_existing or not _has_exposure(row)):
                    query = build_ligand_query(
                        repo_root,
                        row,
                        mapping_csv=mapping_csv,
                        fda_index=fda_index,
                        mapping_rows=mapping_rows,
                    )
                    exposure = spd_exposure_for_query(spd_rows, query=query)
                    if exposure:
                        for col in exposure_cols:
                            value = exposure.get(col)
                            if value not in (None, "") and (
                                not preserve_existing or not _clean_text(row.get(col))
                            ):
                                row[col] = value
                        if _has_exposure(row):
                            filled += 1
                writer.writerow(row)
    LOG.info("wrote %s rows=%d exposure_rows=%d", output_csv, rows, filled)
    return filled


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Fill master_rows exposure fields from local SPD supplementary Cmax/PPB data."
    )
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--spd-supplement-xlsx", required=True, type=Path)
    parser.add_argument("--fda-mapping-csv", default="")
    parser.add_argument("--overwrite-existing", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s [%(name)s] %(message)s",
    )
    fill_master_exposure(
        repo_root=Path(args.repo_root).resolve(),
        input_csv=Path(args.input).resolve(),
        output_csv=Path(args.out).resolve(),
        spd_supplement_xlsx=Path(args.spd_supplement_xlsx).resolve(),
        fda_mapping_csv=str(args.fda_mapping_csv or ""),
        preserve_existing=not bool(args.overwrite_existing),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

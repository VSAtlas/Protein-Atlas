from __future__ import annotations

import argparse
import csv
import logging
import re
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping

import requests  # type: ignore[import-untyped]

from analysis.reporting.ligand_side_effect_cache import (
    build_cache_payload,
    build_ligand_query,
    fetch_ligand_label_exposure_entry,
    fetch_ligand_side_effect_entry,
    load_spd_exposure_rows,
    lookup_ligand_side_effect_entry,
    resolve_mapping_context,
    spd_exposure_for_query,
    write_cache,
)
from config.output_paths import output_root

LOG = logging.getLogger("ligand-side-effects")
_RDK_RE = re.compile(r"\brdk[_ -]?\d+\b", re.IGNORECASE)


def _normalize_text(value: Any) -> str:
    return " ".join(str(value or "").split()).strip()


def _truthy(value: Any) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "y"}


def _has_rdk_identifier(*values: Any) -> bool:
    return any(_RDK_RE.search(str(value or "")) for value in values)


def _read_heatmap_ligands(
    path: Path,
    *,
    fda_only: bool = True,
    include_controls: bool = False,
    rdk_only: bool = True,
) -> list[Dict[str, Any]]:
    rows_by_key: dict[str, Dict[str, Any]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            library = _normalize_text(row.get("library"))
            if fda_only and library.casefold() != "fda":
                continue
            if _truthy(row.get("is_decoy")):
                continue
            if not include_controls and _truthy(row.get("is_control")):
                continue
            ligand_display = _normalize_text(row.get("ligand_display"))
            ligand_base = _normalize_text(row.get("ligand_base"))
            ligand_name = ligand_display or ligand_base
            if rdk_only and not _has_rdk_identifier(ligand_base, ligand_name):
                continue
            if not ligand_name:
                continue
            key = f"{ligand_name.casefold()}|{ligand_base.casefold()}"
            if key in rows_by_key:
                continue
            rows_by_key[key] = {
                "ligand_name": ligand_name,
                "ligand_display": ligand_display or ligand_name,
                "ligand_base": ligand_base,
                "library": library,
            }
    return sorted(rows_by_key.values(), key=lambda item: item["ligand_name"].casefold())


def refresh_cache(
    *,
    repo_root: Path,
    run_id: str,
    heatmap_csv: Path,
    fda_mapping_csv: str,
    limit: int,
    event_limit: int,
    label_limit: int,
    max_side_effects: int,
    sleep_sec: float,
    reuse_existing: bool,
    fda_only: bool,
    include_controls: bool,
    rdk_only: bool,
    fetch_label_exposure: bool,
    spd_supplement_xlsx: str,
) -> Path:
    mapping_csv, fda_index, mapping_rows = resolve_mapping_context(
        repo_root, run_id, fda_mapping_csv
    )
    if mapping_csv is None:
        LOG.warning("FDA mapping CSV unavailable; name-resolution fallback will be limited")
    else:
        LOG.info("FDA mapping CSV loaded path=%s rows=%d", mapping_csv, len(mapping_rows))
    spd_rows = load_spd_exposure_rows(Path(spd_supplement_xlsx).resolve() if spd_supplement_xlsx else None)
    if spd_rows:
        LOG.info("SPD exposure rows loaded path=%s rows=%d", spd_supplement_xlsx, len(spd_rows))

    ligand_rows = _read_heatmap_ligands(
        heatmap_csv,
        fda_only=fda_only,
        include_controls=include_controls,
        rdk_only=rdk_only,
    )
    if limit > 0:
        ligand_rows = ligand_rows[:limit]
    session = requests.Session()
    entries: list[Mapping[str, Any]] = []
    entries_by_query: dict[tuple[Any, ...], Mapping[str, Any]] = {}
    skipped = 0
    reused_query = 0
    hit_count = 0
    for idx, row in enumerate(ligand_rows, start=1):
        cached = (
            lookup_ligand_side_effect_entry(
                repo_root,
                ligand_label=_normalize_text(row.get("ligand_name")),
                ligand_base=_normalize_text(row.get("ligand_base")),
                ligand_display=_normalize_text(row.get("ligand_display")),
            )
            if reuse_existing
            else {}
        )
        exposure_cached = (
            cached.get("openfda_label_status") is not None
            or cached.get("free_cmax_um") not in (None, "")
            or cached.get("cmax_um") not in (None, "")
        )
        query = build_ligand_query(
            repo_root,
            row,
            mapping_csv=mapping_csv,
            fda_index=fda_index,
            mapping_rows=mapping_rows,
        )
        if not query.query_names and not query.brand_names and not query.rxcuis and not query.uniis:
            LOG.info("skip ligand=%s reason=no_query_identifiers", row.get("ligand_name"))
            continue
        if cached and cached.get("side_effects"):
            entry = dict(cached)
            if spd_rows:
                entry.update(spd_exposure_for_query(spd_rows, query=query))
            if fetch_label_exposure and not exposure_cached:
                entry.update(
                    fetch_ligand_label_exposure_entry(
                        query,
                        session=session,
                        label_limit=label_limit,
                        sleep_sec=sleep_sec,
                    )
                )
                # Preserve previously cached side-effect evidence when this pass is
                # only filling label exposure fields.
                for key in (
                    "side_effects",
                    "side_effect_counts",
                    "side_effect_sources",
                    "side_effect_evidence_summary",
                    "safety_buckets",
                    "primary_display_safety",
                    "secondary_safety_buckets",
                    "drug_ae_buckets",
                    "safety_bucket_scores",
                    "openfda_status",
                    "openfda_event_search",
                    "openfda_last_updated",
                ):
                    if key in cached:
                        entry[key] = cached[key]
            entries.append(entry)
            skipped += 1
            continue
        query_signature = (
            tuple(query.rxcuis),
            tuple(query.uniis),
            tuple(query.brand_names),
            tuple(query.query_names),
        )
        prior = entries_by_query.get(query_signature)
        if prior is not None:
            entry = dict(prior)
            entry.update(
                {
                    "ligand_label": query.ligand_label,
                    "ligand_base": query.ligand_base,
                    "ligand_display": query.ligand_display,
                    "query_reused": True,
                    "mapping_used": query.mapping_used,
                }
            )
            reused_query += 1
        else:
            entry = fetch_ligand_side_effect_entry(
                repo_root,
                query,
                session=session,
                event_limit=event_limit,
                label_limit=label_limit,
                max_side_effects=max_side_effects,
                sleep_sec=sleep_sec,
                fetch_label_exposure=fetch_label_exposure,
            )
            if spd_rows:
                # Prefer SPD supplementary exposure fields when present.
                entry.update(spd_exposure_for_query(spd_rows, query=query))
            entries_by_query[query_signature] = entry
        entries.append(entry)
        if entry.get("side_effects"):
            hit_count += 1
        LOG.info(
            "ligand %d/%d name=%s status=%s label_status=%s side_effects=%d free_cmax=%s",
            idx,
            len(ligand_rows),
            query.ligand_display or query.ligand_label,
            entry.get("openfda_status"),
            entry.get("openfda_label_status"),
            len(entry.get("side_effects") or []),
            entry.get("free_cmax_um") or "",
        )
    path = write_cache(repo_root, build_cache_payload(entries))
    LOG.info(
        "wrote ligand side-effect cache path=%s ligands=%d hits=%d reused=%d reused_query=%d",
        path,
        len(entries),
        hit_count,
        skipped,
        reused_query,
    )
    return path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Refresh ligand side-effect and label-exposure annotations from openFDA."
    )
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--run-id", default="pilotstudy")
    parser.add_argument("--heatmap-csv", default="")
    parser.add_argument("--fda-mapping-csv", default="")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--event-limit", type=int, default=50)
    parser.add_argument("--label-limit", type=int, default=5)
    parser.add_argument("--max-side-effects", type=int, default=20)
    parser.add_argument("--sleep-sec", type=float, default=0.25)
    parser.add_argument("--skip-label-exposure", action="store_true")
    parser.add_argument("--include-non-fda", action="store_true")
    parser.add_argument("--include-controls", action="store_true")
    parser.add_argument("--include-non-rdk-fda", action="store_true")
    parser.add_argument("--reuse-existing", action="store_true")
    parser.add_argument(
        "--spd-supplement-xlsx",
        default="",
        help="Optional path to SPD Supplementary Data workbook (MOESM4) to enrich Cmax/PPB coverage.",
    )
    parser.add_argument("--verbose", action="store_true")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s [%(name)s] %(message)s",
    )
    repo_root = Path(args.repo_root).resolve()
    heatmap_csv = (
        Path(args.heatmap_csv).resolve()
        if args.heatmap_csv
        else output_root(repo_root, "data") / args.run_id / "heatmap_input.csv"
    )
    if not heatmap_csv.exists():
        raise FileNotFoundError(f"heatmap CSV not found: {heatmap_csv}")
    refresh_cache(
        repo_root=repo_root,
        run_id=args.run_id,
        heatmap_csv=heatmap_csv,
        fda_mapping_csv=args.fda_mapping_csv,
        limit=max(0, args.limit),
        event_limit=max(1, args.event_limit),
        label_limit=max(1, args.label_limit),
        max_side_effects=max(1, args.max_side_effects),
        sleep_sec=max(0.0, args.sleep_sec),
        reuse_existing=bool(args.reuse_existing),
        fda_only=not bool(args.include_non_fda),
        include_controls=bool(args.include_controls),
        rdk_only=not bool(args.include_non_rdk_fda),
        fetch_label_exposure=not bool(args.skip_label_exposure),
        spd_supplement_xlsx=str(args.spd_supplement_xlsx or ""),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

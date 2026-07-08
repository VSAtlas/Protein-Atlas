from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Iterable

import requests  # type: ignore[import-untyped]

from analysis import pathway_resolver
from analysis.reporting.target_expression_cache import (
    build_cache_payload,
    fetch_target_expression_entry,
    pdb_ids_from_heatmap_csv,
    write_cache,
)
from config.output_paths import output_root

LOG = logging.getLogger("target-expression")


def refresh_cache(
    *,
    repo_root: Path,
    run_id: str,
    heatmap_csv: Path,
    limit: int,
    sleep_sec: float,
) -> Path:
    pdb_ids = pdb_ids_from_heatmap_csv(heatmap_csv)
    if limit > 0:
        pdb_ids = pdb_ids[:limit]
    session = requests.Session()
    pathway_cache = pathway_resolver.Cache(
        cache_dir=repo_root / "pathways" / "cache",
        refresh=False,
        logger=LOG,
    )
    http_client = pathway_resolver.HttpClient()
    entries = []
    hits = 0
    for idx, pdb_id in enumerate(pdb_ids, start=1):
        try:
            uniprots = pathway_resolver.map_pdb_to_uniprots(
                pdb_id, pathway_cache, http_client
            )
        except Exception as exc:
            LOG.warning("skip pdb=%s reason=uniprot_lookup_failed error=%s", pdb_id, exc)
            uniprots = []
        entry = fetch_target_expression_entry(
            uniprots=uniprots,
            session=session,
            sleep_sec=sleep_sec,
        )
        if entry.get("target_tissue_expression_score"):
            hits += 1
        entry["pdb_id"] = pdb_id
        entries.append(entry)
        LOG.info(
            "target %d/%d pdb=%s uniprots=%d status=%s expression=%s",
            idx,
            len(pdb_ids),
            pdb_id,
            len(uniprots),
            entry.get("hpa_query_status"),
            entry.get("target_tissue_expression") or "",
        )
    path = write_cache(repo_root, build_cache_payload(entries))
    LOG.info(
        "wrote target expression cache path=%s run_id=%s targets=%d hits=%d",
        path,
        run_id,
        len(entries),
        hits,
    )
    return path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Refresh target tissue-expression annotations from Human Protein Atlas."
    )
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--run-id", default="pilotstudy")
    parser.add_argument("--heatmap-csv", default="")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--sleep-sec", type=float, default=0.05)
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
        limit=max(0, args.limit),
        sleep_sec=max(0.0, args.sleep_sec),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

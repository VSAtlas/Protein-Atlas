from __future__ import annotations

import argparse
import json
from pathlib import Path

import src.path_router as atlas_path_router

from analysis.external.spl_pk_candidate_review import review_spl_pk_candidates


def _repo_root() -> Path:
    module_file = atlas_path_router.__file__
    if module_file is None:
        return Path.cwd().resolve()
    return Path(module_file).resolve().parents[2]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Create a provenance-first manual review queue for numeric clearance, "
            "absolute bioavailability, observed PK dose, and maximum labeled dose "
            "candidates in the cached DailyMed/openFDA SPL JSON records."
        )
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=(
            _repo_root()
            / "data"
            / "external"
            / "dailymed_spl"
            / "phase1_openfda"
            / "records"
        ),
    )
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument(
        "--drug-id",
        action="append",
        default=None,
        help="Process only this cached drug_id; repeat to select more than one.",
    )
    parser.add_argument(
        "--max-cache-files",
        type=int,
        default=0,
        help="Maximum selected cache payloads to process; zero processes all selected payloads.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.max_cache_files < 0:
        parser.error("--max-cache-files must be zero or greater")
    if not args.cache_dir.is_dir():
        parser.error(f"cache directory not found: {args.cache_dir}")
    manifest = review_spl_pk_candidates(
        cache_dir=args.cache_dir,
        out_dir=args.out_dir,
        drug_ids=args.drug_id,
        max_cache_files=args.max_cache_files,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

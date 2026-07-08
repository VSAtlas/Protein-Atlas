from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Sequence

from analysis.reporting.target_safety_evidence import (
    TargetSpec,
    build_target_safety_cache,
    extract_local_uniprots_from_pdb,
    load_run_targets,
    write_target_safety_cache,
)
from analysis.reporting.target_safety_drift import (
    build_target_safety_drift_summary,
    write_target_safety_drift_summary,
)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Refresh cached target-safety evidence used by heatmap organization."
    )
    parser.add_argument("--repo-root", default=".", help="Repo root containing input_pdbs/ and pathways/cache/.")
    parser.add_argument("--run-id", default="", help="Optional run id to derive targets from data/<run_id>/heatmap_input.csv.")
    parser.add_argument(
        "--pdb-id",
        action="append",
        default=[],
        help="Optional explicit PDB ids to refresh. Can be repeated.",
    )
    parser.add_argument(
        "--target-name",
        action="append",
        default=[],
        help="Optional explicit target names paired positionally with --pdb-id.",
    )
    parser.add_argument(
        "--max-approved-drugs",
        type=int,
        default=12,
        help="Maximum distinct phase-4 known drugs to inspect per target.",
    )
    parser.add_argument(
        "--no-drug-evidence",
        action="store_true",
        help="Disable drug adverse-event evidence and use only direct Open Targets safety liabilities.",
    )
    parser.add_argument(
        "--out",
        default="",
        help="Optional output path. Defaults to pathways/cache/target_safety_aggregated.json under repo root.",
    )
    return parser


def _explicit_targets(repo_root: Path, pdb_ids: Sequence[str], target_names: Sequence[str]) -> list[TargetSpec]:
    out: list[TargetSpec] = []
    for idx, pdb_id in enumerate(pdb_ids):
        name = target_names[idx] if idx < len(target_names) else ""
        out.append(
            TargetSpec(
                pdb_id=str(pdb_id).strip().upper(),
                target_name=str(name).strip(),
                uniprots=tuple(extract_local_uniprots_from_pdb(repo_root, str(pdb_id))),
                gene_symbols=(),
            )
        )
    return [item for item in out if item.pdb_id]


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    repo_root = Path(args.repo_root).resolve()

    targets: list[TargetSpec] = []
    if args.run_id:
        targets.extend(load_run_targets(repo_root, str(args.run_id)))
    if args.pdb_id:
        targets.extend(_explicit_targets(repo_root, list(args.pdb_id), list(args.target_name)))
    deduped: dict[str, TargetSpec] = {}
    for target in targets:
        deduped[target.pdb_id] = target
    if not deduped:
        parser.error("No targets resolved. Provide --run-id or at least one --pdb-id.")

    payload = build_target_safety_cache(
        repo_root,
        list(deduped.values()),
        include_drug_evidence=not bool(args.no_drug_evidence),
        max_approved_drugs=max(1, int(args.max_approved_drugs)),
    )
    previous_payload = {}
    if not args.out:
        default_path = repo_root / "pathways" / "cache" / "target_safety_aggregated.json"
        if default_path.exists():
            try:
                import json

                previous_payload = json.loads(default_path.read_text(encoding="utf-8"))
            except Exception:
                previous_payload = {}
    out_path = write_target_safety_cache(
        repo_root,
        payload,
        out_path=Path(args.out).resolve() if args.out else None,
    )
    drift_summary = build_target_safety_drift_summary(previous_payload, payload)
    drift_path = write_target_safety_drift_summary(repo_root, drift_summary)
    logging.info(
        "[target-safety.refresh.done] targets=%s cached=%s out=%s drift=%s changed=%s",
        len(deduped),
        len(payload.get("entries") or []),
        out_path,
        drift_path,
        drift_summary.get("changed_count"),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

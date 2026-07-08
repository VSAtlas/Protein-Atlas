# -*- coding: utf-8 -*-
import argparse
import sys
from pathlib import Path

from analysis.reporting.fda_name_map import resolve_mapping_csv_path, try_load_fda_index
from analysis.reporting.manifest_utils import load_run_manifest
from analysis.reporting.test_mode_tokens import resolve_test_mode_tokens
from config.output_paths import output_root, run_output_dir

from analysis.reporting.master_schema.collect_rows import _collect_master_rows
from analysis.reporting.master_schema.constants import COMPONENT
from analysis.reporting.master_schema.decoy_helpers import (
    _discover_consensus_files,
    _infer_decoy_prefix_from_tokens,
    _resolve_decoy_prefix,
)
from analysis.reporting.master_schema.fdr_compute import _apply_fdr_to_master_rows
from analysis.reporting.master_schema.io_utils import _configure_logging
from analysis.reporting.master_schema.metadata_loaders import (
    _load_posebusters_map,
    _parse_dud_eval_summary,
)
from analysis.reporting.master_schema.write_exports import _write_master_exports


def main() -> int:
    parser = argparse.ArgumentParser(description="Export master CSV for a docking run.")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--decoy-prefix", default=None)
    parser.add_argument("--fda-mapping-csv", default=None)
    args = parser.parse_args()

    logger = _configure_logging(args.verbose)
    repo_root = Path(args.repo_root).resolve()
    run_id = args.run_id
    decoy_prefix = _resolve_decoy_prefix(repo_root, run_id, args.decoy_prefix)
    tokens = resolve_test_mode_tokens(repo_root, run_id)
    decoy_prefix = _infer_decoy_prefix_from_tokens(decoy_prefix, tokens)
    logger.info(
        "%s action=preflight decoy_prefix=%s tokens=%s",
        COMPONENT,
        decoy_prefix,
        ",".join(tokens),
    )

    mapping_csv = resolve_mapping_csv_path(
        repo_root, run_id, cli_value=args.fda_mapping_csv
    )
    fda_index = try_load_fda_index(mapping_csv)
    if mapping_csv is None:
        logger.warning(
            "%s action=fda_mapping status=missing run_id=%s", COMPONENT, run_id
        )
    elif fda_index is None:
        logger.warning(
            "%s action=fda_mapping status=load_failed path=%s", COMPONENT, mapping_csv
        )

    post_root = output_root(repo_root, "post_docked")
    if not post_root.exists() and (repo_root / "post_docked").exists():
        post_root = repo_root / "post_docked"
    docked_root = output_root(repo_root, "docked")
    if not docked_root.exists() and (repo_root / "docked").exists():
        docked_root = repo_root / "docked"
    data_dir = run_output_dir(repo_root, "data", run_id)
    data_dir.mkdir(parents=True, exist_ok=True)

    output_csv = data_dir / "master_rows.csv"
    if output_csv.exists() and not args.overwrite:
        logger.info("%s action=skip reason=exists path=%s", COMPONENT, output_csv)
        return 0

    # Discovery
    consensus_files = _discover_consensus_files(
        post_root / run_id, tokens, decoy_prefix, logger
    )

    if not consensus_files:
        logger.warning("%s action=exit reason=no_consensus_files", COMPONENT)
        return 0

    # Load metadata sources
    manifest, _manifest_path = load_run_manifest(repo_root, run_id)
    manifest_paths = manifest.get("paths") if isinstance(manifest, dict) else None
    manifest_paths = manifest_paths if isinstance(manifest_paths, dict) else {}
    processed_root = Path(
        str(
            manifest_paths.get("processed_pdb_dir")
            or (output_root(repo_root, "processed_pdbs") / run_id)
        )
    )
    pb_map = _load_posebusters_map(post_root, run_id, decoy_prefix)
    dud_map = _parse_dud_eval_summary(run_id, repo_root, logger, decoy_prefix)

    master_rows, _control_caches = _collect_master_rows(
        consensus_files=consensus_files,
        post_root=post_root,
        repo_root=repo_root,
        run_id=run_id,
        decoy_prefix=decoy_prefix,
        manifest=manifest,
        processed_root=processed_root,
        pb_map=pb_map,
        dud_map=dud_map,
        fda_index=fda_index,
        logger=logger,
    )

    if not master_rows:
        logger.warning("%s action=exit reason=no_rows_collected", COMPONENT)
        return 0

    master_rows, fdr_summary_rows, scorch_fdr_summary_rows = _apply_fdr_to_master_rows(
        master_rows=master_rows,
        post_root=post_root,
        docked_root=docked_root,
        run_id=run_id,
        decoy_prefix=decoy_prefix,
        logger=logger,
    )

    return _write_master_exports(
        master_rows=master_rows,
        output_csv=output_csv,
        data_dir=data_dir,
        fdr_summary_rows=fdr_summary_rows,
        scorch_fdr_summary_rows=scorch_fdr_summary_rows,
        logger=logger,
    )


if __name__ == "__main__":
    sys.exit(main())

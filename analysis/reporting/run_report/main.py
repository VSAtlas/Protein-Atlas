# ruff: noqa: F403, F405
import argparse
from analysis.reporting.run_report.constants import *
from analysis.reporting.run_report.build_report import build_report
from analysis.reporting.run_report.config_resolve import (
    _configure_logging,
    _resolve_filter_invalid,
    _resolve_percent_searched,
    _resolve_report_asset_mode,
    _resolve_report_bool,
    _resolve_report_float,
    _resolve_report_limit,
)
from analysis.reporting.run_report.display import _parse_highlight_queries
from analysis.reporting.run_report.export_io import write_yaml
from analysis.reporting.run_report.heatmap_csv import _write_heatmap_input_csv
from analysis.reporting.run_report.html_report import _write_html_report
from analysis.reporting.run_report.master_csv import (
    _infer_decoy_prefix_from_tokens,
    _normalize_library_slug,
    _resolve_master_csvs,
    _write_split_master_rows,
)
from analysis.reporting.run_report.target_stats import _resolve_decoy_prefix
from analysis.reporting.run_report.pathway_reports import _write_html_reports
from analysis.reporting.run_report.paths import _run_data_dir
from analysis.reporting.run_report.side_effect_caches import _ensure_side_effect_caches
from analysis.reporting.run_report.vina_fallback import (
    _write_heatmap_input_csv_from_vina,
    _write_vina_heatmap_html_report,
)
def main() -> int:
    parser = argparse.ArgumentParser(description="Generate run report YAML")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--overwrite", action="store_true", default=True)
    parser.add_argument("--no-overwrite", action="store_false", dest="overwrite")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--decoy-prefix", default=None)
    parser.add_argument("--top-n", type=int, default=5)
    parser.add_argument("--extended-top-n", type=int, default=5)
    parser.add_argument("--notable-pct", type=float, default=0.01)
    parser.add_argument("--notable-max", type=int, default=25)
    parser.add_argument("--fda-mapping-csv", default=None)
    parser.add_argument("--highlight-ligands", default="imatinib")
    parser.add_argument(
        "--highlight-match",
        choices=["display_contains", "display_exact", "any_contains"],
        default="any_contains",
    )
    parser.add_argument("--highlight-max-per-target", type=int, default=5)
    parser.add_argument("--multi-target-min-targets", type=int, default=2)
    parser.add_argument("--multi-target-max-pct", type=float, default=0.05)
    parser.add_argument("--multi-target-max-hits", type=int, default=50)
    parser.add_argument(
        "--multi-target-sort",
        choices=["worst_pct_then_mean", "mean_pct", "best_rank_sum"],
        default="worst_pct_then_mean",
    )
    parser.add_argument(
        "--emit-html",
        action="store_true",
        default=True,
        dest="emit_html",
    )
    parser.add_argument(
        "--no-emit-html",
        action="store_false",
        dest="emit_html",
    )
    parser.add_argument("--html-path", default=None)
    parser.add_argument(
        "--emit-heatmap-csv",
        action="store_true",
        default=True,
        dest="emit_heatmap_csv",
    )
    parser.add_argument(
        "--no-emit-heatmap-csv",
        action="store_false",
        dest="emit_heatmap_csv",
    )
    parser.add_argument(
        "--emit-pathway-html",
        action="store_true",
        default=False,
        dest="emit_pathway_html",
        help="Also emit one pathway-filtered HTML report per Reactome/pathway group.",
    )
    parser.add_argument(
        "--no-emit-pathway-html",
        action="store_false",
        dest="emit_pathway_html",
        help="Skip pathway-filtered HTML reports and write only the main report.html.",
    )
    parser.add_argument("--heatmap-csv-path", default=None)
    parser.add_argument("--heatmap-top-k", type=int, default=200)
    parser.add_argument(
        "--report-asset-mode",
        choices=["auto", "inline", "relative", "cdn"],
        default="auto",
    )
    parser.add_argument(
        "--report-html-warn-bytes", type=int, default=REPORT_HTML_WARN_BYTES_DEFAULT
    )
    parser.add_argument(
        "--side-effect-cache",
        dest="side_effect_cache",
        choices=["auto", "always", "never"],
        default="auto",
        help=(
            "Refresh ligand openFDA and target/PDB safety caches before HTML "
            "generation when missing/incomplete, always, or never."
        ),
    )
    parser.add_argument(
        "--ligand-side-effect-cache",
        dest="side_effect_cache",
        choices=["auto", "always", "never"],
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--side-effect-cache-min-hits", type=int, default=1)
    parser.add_argument(
        "--ligand-side-effect-cache-min-hits",
        dest="side_effect_cache_min_hits",
        type=int,
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--side-effect-cache-limit", type=int, default=0)
    parser.add_argument(
        "--ligand-side-effect-cache-limit",
        dest="side_effect_cache_limit",
        type=int,
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--side-effect-cache-event-limit", type=int, default=50)
    parser.add_argument(
        "--ligand-side-effect-cache-event-limit",
        dest="side_effect_cache_event_limit",
        type=int,
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--side-effect-cache-label-limit", type=int, default=5)
    parser.add_argument(
        "--ligand-side-effect-cache-label-limit",
        dest="side_effect_cache_label_limit",
        type=int,
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--side-effect-cache-max-terms", type=int, default=20)
    parser.add_argument(
        "--ligand-side-effect-cache-max-side-effects",
        dest="side_effect_cache_max_terms",
        type=int,
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--side-effect-cache-sleep-sec", type=float, default=0.25)
    parser.add_argument(
        "--ligand-side-effect-cache-sleep-sec",
        dest="side_effect_cache_sleep_sec",
        type=float,
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--side-effect-cache-fetch-exposure",
        dest="side_effect_cache_fetch_exposure",
        action="store_true",
        default=False,
        help="Also query openFDA labels for Cmax/protein-binding fields while refreshing.",
    )
    parser.add_argument(
        "--ligand-side-effect-cache-fetch-exposure",
        dest="side_effect_cache_fetch_exposure",
        action="store_true",
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-vina",
        "--vina",
        action="store_true",
        default=False,
        dest="vina_fallback",
        help="Allow heatmap/report fallback from Vina docking summaries when master_rows.csv is unavailable.",
    )
    parser.add_argument("--top-targets-n", type=int, default=3)
    parser.add_argument("--breadth-pct-threshold", type=float, default=0.01)
    parser.add_argument(
        "--coverage-requires-pose-valid",
        action="store_true",
        default=False,
    )
    parser.add_argument("--pct-display-decimals", type=int, default=1)
    parser.add_argument("--tail-median-decimals", type=int, default=0)
    parser.add_argument("--highlights-top-pct", type=float, default=0.01)
    parser.add_argument("--highlights-max-ligands", type=int, default=25)
    parser.add_argument(
        "-combined",
        "--combined",
        action="store_true",
        default=False,
        help="Emit a single combined report instead of per-library reports.",
    )
    args = parser.parse_args()

    logger = _configure_logging(args.verbose)
    repo_root = Path(args.repo_root).resolve()
    run_id = args.run_id
    decoy_prefix = _resolve_decoy_prefix(repo_root, run_id, args.decoy_prefix)
    report_asset_mode = _resolve_report_asset_mode(
        repo_root,
        run_id,
        args.report_asset_mode if "--report-asset-mode" in sys.argv else None,
    )
    tokens = resolve_test_mode_tokens(repo_root, run_id)
    decoy_prefix = _infer_decoy_prefix_from_tokens(decoy_prefix, tokens)
    logger.info(
        "%s action=preflight decoy_prefix=%s tokens=%s",
        COMPONENT,
        decoy_prefix,
        ",".join(tokens),
    )
    percent_searched = _resolve_percent_searched(repo_root)
    notable_pct = args.notable_pct
    multi_target_max_pct = args.multi_target_max_pct
    if percent_searched is not None:
        if "--notable-pct" not in sys.argv:
            notable_pct = percent_searched
        if "--multi-target-max-pct" not in sys.argv:
            multi_target_max_pct = percent_searched
        logger.info(
            "%s action=percent_searched value=%.6f",
            COMPONENT,
            percent_searched,
        )

    filter_invalid = _resolve_filter_invalid(repo_root, run_id)
    if filter_invalid is None:
        filter_invalid = False

    top_n = args.top_n
    extended_top_n = args.extended_top_n
    if "--top-n" not in sys.argv:
        cfg_top_n = _resolve_report_limit(repo_root, run_id, "REPORT_TOP_N")
        if cfg_top_n is not None:
            top_n = cfg_top_n
    if "--extended-top-n" not in sys.argv:
        cfg_extended = _resolve_report_limit(repo_root, run_id, "REPORT_EXTENDED_TOP_N")
        if cfg_extended is not None:
            extended_top_n = cfg_extended

    top_targets_n = args.top_targets_n
    if "--top-targets-n" not in sys.argv:
        cfg_top_targets_n = _resolve_report_limit(repo_root, run_id, "TOP_TARGETS_N")
        if cfg_top_targets_n is not None:
            top_targets_n = cfg_top_targets_n

    breadth_pct_threshold = args.breadth_pct_threshold
    if "--breadth-pct-threshold" not in sys.argv:
        cfg_breadth_pct = _resolve_report_float(
            repo_root, run_id, "BREADTH_PCT_THRESHOLD", allow_percent=True
        )
        if cfg_breadth_pct is not None:
            breadth_pct_threshold = cfg_breadth_pct

    coverage_requires_pose_valid = args.coverage_requires_pose_valid
    if "--coverage-requires-pose-valid" not in sys.argv:
        cfg_coverage_requires = _resolve_report_bool(
            repo_root, run_id, "COVERAGE_REQUIRES_POSE_VALID"
        )
        if cfg_coverage_requires is not None:
            coverage_requires_pose_valid = cfg_coverage_requires

    pct_display_decimals = args.pct_display_decimals
    if "--pct-display-decimals" not in sys.argv:
        cfg_pct_display_decimals = _resolve_report_limit(
            repo_root, run_id, "PCT_DISPLAY_DECIMALS"
        )
        if cfg_pct_display_decimals is not None:
            pct_display_decimals = cfg_pct_display_decimals

    tail_median_decimals = args.tail_median_decimals
    if "--tail-median-decimals" not in sys.argv:
        cfg_tail_median_decimals = _resolve_report_limit(
            repo_root, run_id, "TAIL_MEDIAN_DECIMALS"
        )
        if cfg_tail_median_decimals is not None:
            tail_median_decimals = cfg_tail_median_decimals

    report_html_warn_bytes = args.report_html_warn_bytes
    if "--report-html-warn-bytes" not in sys.argv:
        cfg_report_html_warn_bytes = _resolve_report_limit(
            repo_root, run_id, "REPORT_HTML_WARN_BYTES"
        )
        if cfg_report_html_warn_bytes is not None:
            report_html_warn_bytes = cfg_report_html_warn_bytes

    highlights_top_pct = args.highlights_top_pct
    highlights_max_ligands = args.highlights_max_ligands

    master_csvs = _resolve_master_csvs(repo_root, run_id, tokens, decoy_prefix)
    if not master_csvs:
        if not args.vina_fallback:
            logger.error(
                "%s action=fail error=%s",
                COMPONENT,
                f"Master CSV not found for run_id={run_id} tokens={tokens}",
            )
            return 1
        try:
            top_k = args.heatmap_top_k if args.heatmap_top_k > 0 else None
            data_dir = _run_data_dir(repo_root, run_id)
            if args.heatmap_csv_path:
                fallback_heatmap_path = Path(args.heatmap_csv_path)
            else:
                fallback_heatmap_path = data_dir / "heatmap_input.csv"
            _write_heatmap_input_csv_from_vina(
                repo_root,
                run_id,
                fallback_heatmap_path,
                top_k,
                decoy_prefix,
                filter_invalid,
                args.fda_mapping_csv,
            )
            logger.warning(
                "%s action=vina_fallback status=using_masterless_heatmap run_id=%s",
                COMPONENT,
                run_id,
            )
            if args.emit_html:
                if args.html_path:
                    html_path = Path(args.html_path)
                else:
                    html_path = data_dir / "report.html"
                _ensure_side_effect_caches(
                    repo_root=repo_root,
                    run_id=run_id,
                    heatmap_csv=fallback_heatmap_path,
                    fda_mapping_csv=args.fda_mapping_csv,
                    mode=args.side_effect_cache,
                    min_hits=args.side_effect_cache_min_hits,
                    limit=args.side_effect_cache_limit,
                    event_limit=args.side_effect_cache_event_limit,
                    label_limit=args.side_effect_cache_label_limit,
                    max_side_effects=args.side_effect_cache_max_terms,
                    sleep_sec=args.side_effect_cache_sleep_sec,
                    fetch_label_exposure=args.side_effect_cache_fetch_exposure,
                    logger=logger,
                )
                heatmap_html = render_interactive_heatmap_html(
                    repo_root,
                    run_id,
                    fallback_heatmap_path,
                    top_k=top_k or 100,
                    report_asset_mode=report_asset_mode,
                )
                _write_vina_heatmap_html_report(run_id, html_path, heatmap_html, repo_root)
                logger.info(
                    "%s action=write_html status=ok path=%s library=vina_fallback",
                    COMPONENT,
                    html_path,
                )
            return 0
        except _REPORT_HTML_BUILD_ERRORS as exc:
            logger.error("%s action=fail error=%s", COMPONENT, exc, exc_info=True)
            return 1

    jobs: List[Tuple[Optional[str], List[Path]]] = []
    if args.combined:
        jobs.append((None, master_csvs))
    else:
        split_map = _write_split_master_rows(
            repo_root, run_id, decoy_prefix, master_csvs
        )
        split_keys = sorted(split_map.keys())
        if len(split_keys) > 1:
            for lib_slug in split_keys:
                jobs.append((lib_slug, [split_map[lib_slug]]))
        elif len(split_keys) == 1:
            only = split_keys[0]
            jobs.append((None, [split_map[only]]))
        else:
            jobs.append((None, master_csvs))

    try:
        highlight_queries = _parse_highlight_queries(args.highlight_ligands)
        if args.highlight_max_per_target > 0:
            highlight_queries = highlight_queries[: args.highlight_max_per_target]
        else:
            highlight_queries = []
        top_k = args.heatmap_top_k if args.heatmap_top_k > 0 else None

        for lib_suffix, job_master_csvs in jobs:
            suffix = _normalize_library_slug(lib_suffix) if lib_suffix else ""
            data_dir = _run_data_dir(repo_root, run_id)
            if suffix:
                yaml_path = data_dir / f"report_{suffix}.yaml"
                default_heatmap_path = data_dir / f"heatmap_input_{suffix}.csv"
                default_html_path = data_dir / f"report_{suffix}.html"
            else:
                yaml_path = data_dir / "report.yaml"
                default_heatmap_path = data_dir / "heatmap_input.csv"
                default_html_path = data_dir / "report.html"

            if yaml_path.exists() and not args.overwrite:
                logger.info(
                    "%s action=skip reason=exists path=%s", COMPONENT, yaml_path
                )
                continue

            report = build_report(
                run_id,
                repo_root,
                top_n=top_n,
                extended_top_n=extended_top_n,
                notable_pct=notable_pct,
                notable_max=args.notable_max,
                decoy_prefix=decoy_prefix,
                fda_mapping_csv=args.fda_mapping_csv,
                highlight_ligands=args.highlight_ligands,
                highlight_match=args.highlight_match,
                highlight_max_per_target=args.highlight_max_per_target,
                multi_target_min_targets=args.multi_target_min_targets,
                multi_target_max_pct=multi_target_max_pct,
                multi_target_max_hits=args.multi_target_max_hits,
                multi_target_sort=args.multi_target_sort,
                filter_invalid=filter_invalid,
                top_targets_n=top_targets_n,
                breadth_pct_threshold=breadth_pct_threshold,
                coverage_requires_pose_valid=coverage_requires_pose_valid,
                pct_display_decimals=pct_display_decimals,
                tail_median_decimals=tail_median_decimals,
                highlights_top_pct=highlights_top_pct,
                highlights_max_ligands=highlights_max_ligands,
                test_mode_tokens=tokens,
                master_csvs_override=job_master_csvs,
            )
            write_yaml(report, yaml_path)
            logger.info(
                "%s action=write status=ok path=%s library=%s",
                COMPONENT,
                yaml_path,
                suffix or "combined",
            )

            heatmap_path: Optional[Path] = None
            if args.emit_heatmap_csv:
                if args.heatmap_csv_path and len(jobs) == 1:
                    heatmap_path = Path(args.heatmap_csv_path)
                else:
                    heatmap_path = default_heatmap_path
                try:
                    _write_heatmap_input_csv(
                        repo_root,
                        run_id,
                        heatmap_path,
                        top_k,
                        args.fda_mapping_csv,
                        filter_invalid,
                        master_csv_override=job_master_csvs[0]
                        if job_master_csvs
                        else None,
                    )
                    logger.info(
                        "%s action=write_heatmap_csv status=ok path=%s library=%s",
                        COMPONENT,
                        heatmap_path,
                        suffix or "combined",
                    )
                    try:
                        publication_path, readiness_path = write_publication_exports(
                            heatmap_path,
                            job_master_csvs[0] if job_master_csvs else data_dir / "master_rows.csv",
                            data_dir,
                            suffix=suffix,
                        )
                        logger.info(
                            "%s action=write_publication_exports status=ok publication=%s readiness=%s library=%s",
                            COMPONENT,
                            publication_path,
                            readiness_path,
                            suffix or "combined",
                        )
                    except _REPORT_PUBLICATION_PIPELINE_ERRORS as exc:
                        logger.warning(
                            "%s action=write_publication_exports status=failed heatmap=%s error=%s library=%s",
                            COMPONENT,
                            heatmap_path,
                            exc,
                            suffix or "combined",
                            exc_info=True,
                        )
                except _REPORT_PUBLICATION_PIPELINE_ERRORS as exc:
                    logger.warning(
                        "%s action=write_heatmap_csv status=failed path=%s error=%s library=%s",
                        COMPONENT,
                        heatmap_path,
                        exc,
                        suffix or "combined",
                        exc_info=True,
                    )
                    heatmap_path = None

            if args.emit_html:
                if args.html_path and len(jobs) == 1:
                    html_path = Path(args.html_path)
                else:
                    html_path = default_html_path
                try:
                    heatmap_source_for_html = heatmap_path
                    if heatmap_source_for_html is None:
                        if args.heatmap_csv_path:
                            heatmap_source_for_html = Path(args.heatmap_csv_path)
                        elif job_master_csvs:
                            heatmap_source_for_html = job_master_csvs[0]
                    _ensure_side_effect_caches(
                        repo_root=repo_root,
                        run_id=run_id,
                        heatmap_csv=heatmap_source_for_html,
                        fda_mapping_csv=args.fda_mapping_csv,
                        mode=args.side_effect_cache,
                        min_hits=args.side_effect_cache_min_hits,
                        limit=args.side_effect_cache_limit,
                        event_limit=args.side_effect_cache_event_limit,
                        label_limit=args.side_effect_cache_label_limit,
                        max_side_effects=args.side_effect_cache_max_terms,
                        sleep_sec=args.side_effect_cache_sleep_sec,
                        fetch_label_exposure=args.side_effect_cache_fetch_exposure,
                        logger=logger,
                    )
                    if args.emit_pathway_html:
                        pathway_htmls = _write_html_reports(
                            report,
                            html_path,
                            highlight_queries,
                            repo_root,
                            heatmap_source_override=heatmap_source_for_html,
                            report_asset_mode=report_asset_mode,
                            report_html_warn_bytes=report_html_warn_bytes,
                        )
                    else:
                        _write_html_report(
                            report,
                            html_path,
                            highlight_queries,
                            repo_root,
                            heatmap_source_override=heatmap_source_for_html,
                            report_asset_mode=report_asset_mode,
                            report_html_warn_bytes=report_html_warn_bytes,
                        )
                        pathway_htmls = []
                    logger.info(
                        "%s action=write_html status=ok path=%s library=%s",
                        COMPONENT,
                        html_path,
                        suffix or "combined",
                    )
                    if pathway_htmls:
                        logger.info(
                            "%s action=write_html_pathways status=ok count=%d library=%s",
                            COMPONENT,
                            len(pathway_htmls),
                            suffix or "combined",
                        )
                except _REPORT_HTML_BUILD_ERRORS as exc:
                    logger.warning(
                        "%s action=write_html status=failed path=%s error=%s library=%s",
                        COMPONENT,
                        html_path,
                        exc,
                        suffix or "combined",
                        exc_info=True,
                    )
    except _REPORT_MAIN_ERRORS as e:
        logger.error("%s action=fail error=%s", COMPONENT, e, exc_info=True)
        return 1

    return 0


if __name__ == '__main__':
    sys.exit(main())

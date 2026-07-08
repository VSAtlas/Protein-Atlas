from __future__ import annotations

import argparse
import sys
from typing import Sequence

from config.output_paths import run_output_dir

from cli.qol.status import _cmd_status
from cli.qol import _bindings

def _cmd_report(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="atlas report",
        description="Export master rows and generate the canonical run report HTML/YAML artifacts.",
    )
    parser.add_argument("run_id")
    parser.add_argument("--vina", action="store_true", help="Allow Vina summary fallback.")
    parser.add_argument("--combined", action="store_true")
    parser.add_argument(
        "--skip-master-export",
        action="store_true",
        help="Do not run the master schema export before report generation.",
    )
    parser.add_argument(
        "--overwrite-master",
        action="store_true",
        help="Regenerate outputs/data/<RUN_ID>/master_rows.csv before reporting.",
    )
    parser.add_argument("--open", action="store_true", dest="open_after")
    parser.add_argument("extra", nargs=argparse.REMAINDER)
    args = parser.parse_args(list(argv))

    if not args.skip_master_export:
        export_rc = _run_master_export(
            args.run_id,
            overwrite=args.overwrite_master,
            verbose=False,
        )
        if export_rc != 0:
            return export_rc

    forwarded = ["--run-id", args.run_id, "--repo-root", str(_bindings.repo_root())]
    if args.vina:
        forwarded.append("--vina")
    if args.combined:
        forwarded.append("--combined")
    forwarded.extend(args.extra)
    old_argv = sys.argv[:]
    rc = 1
    try:
        sys.argv = ["atlas report", *forwarded]
        from analysis.reporting import run_report_core

        rc = int(run_report_core.main())
    finally:
        sys.argv = old_argv
    html_path = run_output_dir(_bindings.repo_root(), "data", args.run_id) / "report.html"
    if rc == 0 and html_path.exists():
        print(f"Report HTML: {html_path}")
        if args.open_after:
            print("Open this path in your browser; GUI launching is not done from this CLI.")
    return rc


def _run_master_export(
    run_id: str,
    *,
    overwrite: bool = False,
    verbose: bool = False,
    decoy_prefix: str | None = None,
    fda_mapping_csv: str | None = None,
) -> int:
    forwarded = ["--run-id", str(run_id), "--repo-root", str(_bindings.repo_root())]
    if overwrite:
        forwarded.append("--overwrite")
    if verbose:
        forwarded.append("--verbose")
    if decoy_prefix:
        forwarded.extend(["--decoy-prefix", decoy_prefix])
    if fda_mapping_csv:
        forwarded.extend(["--fda-mapping-csv", fda_mapping_csv])
    old_argv = sys.argv[:]
    rc = 1
    try:
        sys.argv = ["atlas analysis export-master", *forwarded]
        from analysis.reporting import master_schema_export

        rc = int(master_schema_export.main())
    finally:
        sys.argv = old_argv
    out_path = run_output_dir(_bindings.repo_root(), "data", str(run_id)) / "master_rows.csv"
    if out_path.exists():
        print(f"Master rows CSV: {out_path}")
    return rc


def _cmd_analysis(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="atlas analysis",
        description="Discoverable wrappers for maintained analysis/reporting workflows.",
    )
    sub = parser.add_subparsers(dest="analysis_cmd", required=True)
    report = sub.add_parser("report", help="Export master rows and generate the run report.")
    report.add_argument("run_id")
    report.add_argument("--vina", action="store_true")
    report.add_argument("--combined", action="store_true")
    report.add_argument("--overwrite-master", action="store_true")
    report.add_argument(
        "--status-html",
        action="store_true",
        help="Also write outputs/data/<RUN_ID>/status.html.",
    )
    export = sub.add_parser(
        "export-master",
        help="Export outputs/data/<RUN_ID>/master_rows.csv; normally included in `atlas analysis report`.",
    )
    export.add_argument("run_id")
    export.add_argument("--overwrite", action="store_true")
    export.add_argument("--verbose", action="store_true")
    export.add_argument("--decoy-prefix")
    export.add_argument("--fda-mapping-csv")
    from cli.analysis_qol import add_analysis_subcommands, run_analysis_subcommand

    add_analysis_subcommands(sub)
    args = parser.parse_args(list(argv))

    if args.analysis_cmd == "report":
        report_args = [args.run_id]
        if args.vina:
            report_args.append("--vina")
        if args.combined:
            report_args.append("--combined")
        if args.overwrite_master:
            report_args.append("--overwrite-master")
        rc = _cmd_report(report_args)
        if rc == 0 and args.status_html:
            rc = _cmd_status([args.run_id, "--html"])
        return rc

    if args.analysis_cmd == "export-master":
        return _run_master_export(
            args.run_id,
            overwrite=args.overwrite,
            verbose=args.verbose,
            decoy_prefix=args.decoy_prefix,
            fda_mapping_csv=args.fda_mapping_csv,
        )
    routed_rc = run_analysis_subcommand(args, _bindings.repo_root())
    if routed_rc is not None:
        return routed_rc
    return 1

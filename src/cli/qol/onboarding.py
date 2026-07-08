from __future__ import annotations

import argparse
import csv
import html
import importlib
import json
import shutil
import sys
from pathlib import Path
from typing import Any, Sequence

from config.output_paths import output_root
from config.tool_resolver import resolve_tool

from cli.qol.constants import DEMO_RUN_ID, DEMO_ROWS
from cli.qol.utils import _display_path
from cli.qol import _bindings
from cli.qol.status import _cmd_status
from cli.qol.pipeline import _cmd_run

def _cmd_help(argv: Sequence[str]) -> int:
    _ = argv
    print(
        "\n".join(
            [
                "Atlas user commands:",
                "  atlas init                 create/update config.txt and auto-detect local tools",
                "  atlas demo                 write a tiny report demo under outputs/data/atlas_demo/",
                "  atlas demo --status-html   also write the demo status dashboard HTML",
                "  atlas smoke public         run the publication-safe public smoke",
                "  atlas setup-report         summarize install readiness and missing optional tools",
                "  atlas new-run              validate inputs and build the first-run command",
                "  atlas first-run            select/install a target, prepare ligands, and plan/run",
                "  atlas status [RUN_ID]      summarize run progress and key output paths",
                "  atlas status RUN_ID --errors --explain --watch 30",
                "  atlas runs                 list recent run IDs and report status",
                "  atlas report RUN_ID        export master rows and generate outputs/data/<RUN_ID>/report.html",
                "  atlas analysis report RUN_ID",
                "  atlas analysis dud-eval RUN_ID",
                "  atlas analysis throughput integrity RUN_ID",
                "  atlas analysis interactions export RUN_ID",
                "  atlas ml train --run-id RUN_ID",
                "  atlas ml audit --dataset DATASET --label LABEL",
                "  atlas ml doctor --run-id RUN_ID --strict",
                "  atlas ml hpo --dataset DATASET --label LABEL",
                "  atlas ml chemprop --dataset DATASET --label LABEL",
                "  atlas ml reinvent --generated-smiles generated.csv",
                "  atlas debug RUN_ID --deep  summarize structured run errors and operator actions",
                "  atlas artifacts measure --run-id RUN_ID",
                "  atlas throughput bench --profile smoke --mode local",
                "  atlas throughput bench --profile smoke --mode slurm-sim --array 0-2%3 --cpus-per-task 4",
                "  atlas reproduce bundle RUN_ID",
                "  atlas slurm submit --run-id RUN_ID --array 0-31%4 --cpus-per-task 8 --dry-run",
                "  atlas slurm submit --bench2-canary --run-id RUN_ID --dry-run",
                "  atlas slurm progress RUN_ID",
                "  atlas slurm finalize RUN_ID --reconcile-only",
                "  atlas screenshot RUN_ID --pdb 1ABC --top 20",
                "  atlas targets guide        prompt through gene target search/install",
                "  atlas targets genes        write a reusable gene/category CSV",
                "  atlas targets from-genes   select PDB structures from gene symbols",
                "  atlas targets search       full-text RCSB target search",
                "  atlas targets install      select and download PDB targets into input_pdbs/",
                "  atlas run-panel kinases --ligands chembl --fast",
                "  atlas ligands sources      list built-in downloadable ligand libraries",
                "  atlas ligands install chembl download/cache SDFs and prepare PDBQT ligands",
                "  atlas doctor --pdb 1ABC --ligands chembl",
                "  atlas chembl --pdb 1ABC   install ChEMBL if needed, then dock ChEMBL",
                "  atlas dev verify --full --fix --smoke",
                "  atlas run [OPTIONS]        run the legacy docking pipeline explicitly",
                "",
                "Existing pipeline flags still work directly, for example:",
                "  atlas --pdb 1BN1 --fast",
                "  atlas --verify-tools",
                "  atlas --doctor",
            ]
        )
    )
    return 0


def _cmd_init(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="atlas init",
        description="Create config.txt and auto-detect local external tool paths.",
    )
    parser.add_argument("--config", default="config.txt")
    parser.add_argument("--example", default="config.example.txt")
    parser.add_argument("--set", action="append", default=[], metavar="KEY=VALUE")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--no-detect", action="store_true")
    parser.add_argument("--verify-tools", action="store_true")
    parser.add_argument(
        "--write-example-inputs",
        action="store_true",
        help="Write lightweight docs/examples onboarding notes without touching runtime inputs.",
    )
    args = parser.parse_args(list(argv))

    from tools.installers import configure_local_tools

    root = _bindings.repo_root()
    config_path = (root / args.config).resolve() if not Path(args.config).is_absolute() else Path(args.config)
    example_path = (root / args.example).resolve() if not Path(args.example).is_absolute() else Path(args.example)
    if not config_path.exists():
        if not example_path.exists():
            parser.error(f"example config not found: {example_path}")
        if args.dry_run:
            print(f"Would create {config_path} from {example_path}")
        else:
            shutil.copyfile(example_path, config_path)
            print(f"Created {config_path} from {example_path}")

    forwarded = ["--config", str(config_path), "--example", str(example_path)]
    for item in args.set:
        forwarded.extend(["--set", item])
    if args.dry_run:
        forwarded.append("--dry-run")
    if args.no_detect:
        forwarded.append("--no-detect")
    rc = configure_local_tools.main(forwarded)
    print("Next: run `atlas doctor` for a quick environment check.")
    if args.verify_tools:
        print("Next: run `atlas --verify-tools` after licensed tools are installed.")
    if args.write_example_inputs:
        _write_example_input_docs(root, dry_run=args.dry_run)
    return int(rc)


def _write_example_input_docs(root: Path, *, dry_run: bool = False) -> None:
    examples_dir = root / "docs" / "examples"
    examples_readme = examples_dir / "new_user_inputs.md"
    text = (
        "# New User Inputs\n\n"
        "Atlas reads protein structures from the configured input root, which "
        "defaults to `input_pdbs/`. Put receptor structure files there as "
        "`<PDB_ID>.pdb`, for example `1ABC.pdb`.\n\n"
        "Do not put generated docking outputs in source control. Run "
        "`atlas new-run --pdb <PDB_ID> --small-library --fast --dry-run` to "
        "validate that Atlas can see the target before launching docking.\n\n"
        "Prepared ligand libraries live under the configured ligand library root. "
        "Start with `atlas demo` to inspect the report UI without requiring a real "
        "library or external docking tools.\n"
    )
    if dry_run:
        print(f"Would write example guide: {examples_readme}")
        return
    examples_dir.mkdir(parents=True, exist_ok=True)
    if examples_readme.exists():
        print(f"Example guide already exists: {examples_readme}")
        return
    examples_readme.write_text(text, encoding="utf-8")
    print(f"Wrote example guide: {examples_readme}")


def _cmd_demo(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="atlas demo",
        description="Create a tiny report demo without external docking tools.",
    )
    parser.add_argument("--run-id", default=DEMO_RUN_ID)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--top-k", type=int, default=25)
    parser.add_argument(
        "--report-asset-mode",
        choices=["auto", "inline", "relative", "cdn"],
        default="inline",
    )
    parser.add_argument(
        "--status-html",
        action="store_true",
        help="Also write the static status dashboard for the demo run.",
    )
    args = parser.parse_args(list(argv))

    root = _bindings.repo_root()
    data_dir = output_root(root, "data") / str(args.run_id)
    csv_path = data_dir / "heatmap_input.csv"
    html_path = data_dir / "report.html"
    if html_path.exists() and not args.force:
        print(f"Demo already exists: {html_path}")
        print("Use `atlas demo --force` to regenerate it.")
        if args.status_html:
            return _cmd_status([args.run_id, "--html"])
        return 0

    data_dir.mkdir(parents=True, exist_ok=True)
    fieldnames = list(DEMO_ROWS[0].keys())
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(DEMO_ROWS)

    from analysis.reporting import heatmap_html as heatmap_module

    original_uniprots = heatmap_module._fetch_target_uniprots_by_pdb
    original_pathways = heatmap_module._fetch_pathway_memberships_by_pdb
    heatmap_module._fetch_target_uniprots_by_pdb = _offline_demo_uniprots
    heatmap_module._fetch_pathway_memberships_by_pdb = _offline_demo_pathways
    try:
        heatmap_html = heatmap_module.render_interactive_heatmap_html(
            root,
            args.run_id,
            csv_path,
            top_k=args.top_k,
            report_asset_mode=args.report_asset_mode,
        )
    finally:
        heatmap_module._fetch_target_uniprots_by_pdb = original_uniprots
        heatmap_module._fetch_pathway_memberships_by_pdb = original_pathways
    title = f"Atlas demo report: {args.run_id}"
    body = (
        "<!doctype html>\n"
        "<html>\n"
        "<head>\n"
        '  <meta charset="utf-8">\n'
        '  <meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f"  <title>{html.escape(title)}</title>\n"
        "</head>\n"
        "<body>\n"
        "<section class=\"section\" id=\"heatmap\">\n"
        f"{heatmap_html}\n"
        "</section>\n"
        "</body>\n"
        "</html>\n"
    )
    html_path.write_text(body, encoding="utf-8")
    print(f"Demo heatmap CSV: {csv_path}")
    print(f"Demo report HTML: {html_path}")
    if args.status_html:
        return _cmd_status([args.run_id, "--html"])
    print("Next: open the HTML report in a browser or run `atlas status atlas_demo`.")
    return 0


def _cmd_smoke(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="atlas smoke",
        description="Run maintained smoke workflows.",
    )
    sub = parser.add_subparsers(dest="smoke_cmd", required=True)
    sub.add_parser("public", help="Run the publication-safe public smoke.")
    sub.add_parser("internal", help="Run the legacy-compatible --test -fast smoke.")
    args = parser.parse_args(list(argv))

    if args.smoke_cmd == "internal":
        return _cmd_run(["--test", "-fast"])
    root = _bindings.repo_root()
    script = root / "tools" / "public_smoke_check.sh"
    return _bindings.subprocess().run(["bash", str(script)], cwd=root, check=False).returncode


def _offline_demo_uniprots(workspace_root: Path, pdb_ids: Any) -> dict[str, list[str]]:
    _ = workspace_root
    return {str(pdb_id).upper(): [] for pdb_id in pdb_ids}


def _offline_demo_pathways(
    workspace_root: Path,
    pdb_ids: Any,
    *,
    uniprots_by_pdb: dict[str, list[str]] | None = None,
) -> dict[str, list[str]]:
    _ = workspace_root
    _ = uniprots_by_pdb
    return {str(pdb_id).upper(): [] for pdb_id in pdb_ids}


def _cmd_setup_report(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="atlas setup-report",
        description="Summarize install readiness and missing optional/BYOL tools.",
    )
    parser.add_argument("--json", action="store_true", dest="as_json")
    args = parser.parse_args(list(argv))

    report = _bindings.build_setup_report()
    if args.as_json:
        print(json.dumps(report, indent=2, sort_keys=True, default=str))
        return 0 if report["ready"]["demo"] else 1
    _print_setup_report(report)
    return 0 if report["ready"]["demo"] else 1


def _build_setup_report() -> dict[str, Any]:
    root = _bindings.repo_root()
    cfg = _bindings.load_effective_config()
    imports = _check_imports(
        {
            "yaml": True,
            "numpy": True,
            "pandas": True,
            "rdkit": True,
            "Bio": True,
            "meeko": True,
            "openbabel": False,
        }
    )
    tools = {
        "vina": resolve_tool(cfg, "VINA_EXE", "vina"),
        "openbabel": resolve_tool(cfg, "OPENBABEL_PATH", "obabel"),
        "p2rank": resolve_tool(cfg, "P2RANK_PATH", "prank"),
        "scorch": resolve_tool(cfg, "SCORCH", "scorch.py"),
    }
    for tool_name, tool_row in tools.items():
        tool_row["fix"] = (
            "none" if tool_row.get("resolved_path") else _tool_fix(tool_name)
        )
    demo_ready = _all_required_ok(imports, required_only=True)
    report_ready = all(
        imports[name]["ok"] for name in ("yaml", "numpy", "pandas") if name in imports
    )
    vina_ready = bool(tools["vina"].get("resolved_path")) and (
        bool(tools["openbabel"].get("resolved_path")) or imports["openbabel"]["ok"]
    )
    prep_ready = bool(imports["meeko"]["ok"])
    optional_missing = [
        label
        for label, key in (
            ("P2Rank pocket planner", "p2rank"),
            ("SCORCH rescoring", "scorch"),
        )
        if not tools[key].get("resolved_path")
    ]
    blocked = []
    if not demo_ready:
        blocked.append("demo: missing required Python imports")
    if not vina_ready:
        blocked.append("Vina docking: missing Vina or Open Babel")
    if not prep_ready:
        blocked.append("full receptor prep: missing Meeko")
    return {
        "repo_root": str(root),
        "config_path": str(root / "config.txt"),
        "python_executable": sys.executable,
        "imports": imports,
        "tools": tools,
        "ready": {
            "demo": demo_ready,
            "report": report_ready,
            "vina_docking": vina_ready,
            "full_receptor_prep": prep_ready,
        },
        "optional_missing": optional_missing,
        "blocked": blocked,
        "next": _setup_next_steps(demo_ready, vina_ready, prep_ready),
    }


def _check_imports(modules: dict[str, bool]) -> dict[str, dict[str, Any]]:
    rows: dict[str, dict[str, Any]] = {}
    for module_name, required in modules.items():
        try:
            importlib.import_module(module_name)
        except Exception as exc:
            rows[module_name] = {
                "ok": False,
                "required": required,
                "detail": str(exc),
                "fix": _import_fix(module_name),
            }
            continue
        rows[module_name] = {
            "ok": True,
            "required": required,
            "detail": "ok",
            "fix": "none",
        }
    return rows


def _all_required_ok(rows: dict[str, dict[str, Any]], *, required_only: bool) -> bool:
    for row in rows.values():
        if required_only and not row.get("required"):
            continue
        if not row.get("ok"):
            return False
    return True


def _import_fix(module_name: str) -> str:
    conda_names = {
        "yaml": "pyyaml",
        "Bio": "biopython",
        "rdkit": "rdkit",
        "openbabel": "openbabel",
        "meeko": "meeko",
    }
    package = conda_names.get(module_name, module_name)
    return f"rerun: bash tools/installers/install_atlas_publication_stack.sh (or install conda-forge::{package})"


def _tool_fix(tool_name: str) -> str:
    fixes = {
        "vina": "rerun one-command installer, or set VINA_EXE=/path/to/vina",
        "openbabel": "rerun one-command installer, or set OPENBABEL_PATH=/path/to/obabel",
        "p2rank": "optional/BYOL: register P2Rank with install_p2rank.sh --p2rank-root",
        "scorch": "optional: register SCORCH with install_scorch.sh --scorch-source",
    }
    return fixes.get(tool_name, "set the matching config.txt path")


def _setup_next_steps(demo_ready: bool, vina_ready: bool, prep_ready: bool) -> list[str]:
    if not demo_ready:
        return ["rerun the one-command installer or inspect `atlas doctor`"]
    steps = ["run `atlas demo` then `atlas status atlas_demo`"]
    if vina_ready and prep_ready:
        steps.append("put a PDB in input_pdbs/ and run `atlas new-run --pdb <ID> --small-library --fast --dry-run`")
    elif vina_ready:
        steps.append(
            "Vina is present; configure BYOL receptor prep tools before full docking"
        )
    else:
        steps.append("configure Vina/Open Babel for real docking, or continue with report demos")
    return steps


def _print_setup_report(report: dict[str, Any]) -> None:
    root = Path(str(report.get("repo_root") or _bindings.repo_root())).expanduser()
    print("Atlas setup report")
    print(f"python: {_display_path(report['python_executable'], root=root)}")
    print(f"config: {_display_path(report['config_path'], root=root)}")
    print("readiness:")
    for key, value in report["ready"].items():
        print(f"  {key}: {'ready' if value else 'not ready'}")
    print("detected tools:")
    tools: dict[str, Any] = report["tools"]
    for name in sorted(tools):
        path = _display_path(tools[name].get("resolved_path") or "", root=root)
        source = str(tools[name].get("source") or "missing")
        print(f"  {name}: {path or 'missing'} ({source})")
        if not path:
            print(f"    fix: {tools[name].get('fix')}")
    optional_missing = report.get("optional_missing") or []
    if optional_missing:
        print("optional/BYOL tools not configured:")
        for item in optional_missing:
            print(f"  {item}")
        print("Atlas demo/report workflows still run without optional/BYOL tools.")
    blocked = report.get("blocked") or []
    if blocked:
        print("blocked:")
        for item in blocked:
            print(f"  {item}")
    for step in report.get("next") or []:
        print(f"next: {step}")

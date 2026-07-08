from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.spd_panel import stage_spd_full_panel
from src.cli.target_install import write_target_genes_csv
from src.cli import target_install


def _coerce_target_rows(
    rows: list[dict[str, object]],
) -> list[dict[str, str]]:
    seen: set[str] = set()
    cleaned: list[dict[str, str]] = []
    for row in rows:
        gene = str(row.get("target_id", "")).strip().upper()
        if not gene or gene in seen:
            continue
        seen.add(gene)
        cleaned.append({"target_id": gene})
    return cleaned


def _collect_spd_targets_for_install(args: argparse.Namespace) -> list[dict[str, str]]:
    from analysis.ml.spd_panel import build_spd_target_assay_manifest

    sources: list[Path] = []
    if args.spd_assays is not None and args.spd_assays.exists():
        sources.append(args.spd_assays)
    elif args.spd is not None and args.spd.exists():
        sources.append(args.spd)

    rows: list[dict[str, str]] = []
    for source in sources:
        if source is None or not source.exists():
            continue
        rows.extend(_coerce_target_rows(build_spd_target_assay_manifest(source).to_dict("records")))
    if rows:
        return rows

    pair = pd.read_csv(args.pair_table, low_memory=False)
    for target_id in pair.get("target_id", pd.Series(dtype="object")).astype(str):
        target_id = target_id.strip().upper()
        if target_id and target_id not in {r["target_id"] for r in rows}:
            rows.append({"target_id": target_id})
    return rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Stage full SPD panel coverage against Atlas pairs.")
    parser.add_argument("--spd", required=True, type=Path)
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--mapping", type=Path, default=None)
    parser.add_argument(
        "--spd-assays",
        type=Path,
        default=None,
        help="Optional SPD assay catalog for target→assay provenance mapping.",
    )
    parser.add_argument(
        "--mode",
        choices=("strict", "relaxed", "both"),
        default="strict",
        help="strict: strict holo-only production targets; relaxed: apo+relaxed ligand filter;",
    )
    parser.add_argument(
        "--install-targets",
        action="store_true",
        help="Run target selection and download (or dry-run by default).",
    )
    parser.add_argument(
        "--max-candidates",
        type=int,
        default=200,
        help="PDB query candidate cap per SPD target.",
    )
    parser.add_argument(
        "--max-return",
        type=int,
        default=1,
        help="PDB selections per target; set to 1 for one PDB per target.",
    )
    parser.add_argument("--margin-strong", type=float, default=10.0)
    parser.add_argument("--margin-weak", type=float, default=100.0)
    parser.add_argument(
        "--install-dry-run",
        action="store_true",
        default=False,
        help="Select candidates without downloading (safe default).",
    )
    args = parser.parse_args(argv)
    # local import to keep this CLI lightweight when used only for staging
    from src.config.runtime_config import load_inputs

    cfg = load_inputs()
    out_root = args.out_dir
    modes: list[str] = ["strict", "relaxed"] if args.mode == "both" else [args.mode]
    mode_results: dict[str, Any] = {}
    strict_target_count: int | None = None
    relaxed_target_count: int | None = None
    for mode in modes:
        mode_out = out_root / ("spd_full_panel_" + mode)
        include_assays = True if mode == "strict" else False
        manifest = stage_spd_full_panel(
            args.spd,
            args.pair_table,
            mode_out,
            mapping_path=args.mapping,
            margin_strong=args.margin_strong,
            margin_weak=args.margin_weak,
            include_assay_manifest=include_assays,
        )
        mode_results[mode] = manifest

        if args.spd_assays is not None and args.spd_assays.exists():
            manifest["target_assay_map"] = str(mode_out / "spd_target_assay_manifest.csv")

        # Build explicit target→assay manifest first if assay source is available.
        target_rows = _collect_spd_targets_for_install(args)
        # fallback to full-panel-only mode keeps the command operable even without assay metadata
        if not target_rows:
            pair = pd.read_csv(args.pair_table, low_memory=False)
            for t in sorted(set(pair.get("target_id", pd.Series(dtype="object")).astype(str).str.upper())):
                t = str(t).strip()
                if t and t != "NAN":
                    target_rows.append({"target_id": t})

        if target_rows:
            genes = [row["target_id"] for row in target_rows if row.get("target_id")]
            genes_path = mode_out / "spd_target_genes.csv"
            write_target_genes_csv(genes_path, genes, category="SPD")
            target_install_result = target_install.install_targets(
                cfg,
                genes_file=genes_path,
                selected_out=mode_out / "spd_target_selected.csv",
                manifest_out=mode_out / "spd_target_install_manifest.json",
                missing_out=mode_out / "spd_target_install_missing.csv",
                out_dir=mode_out / "structures",
                dry_run=not args.install_targets or bool(args.install_dry_run),
                max_return=args.max_return,
                max_candidates=args.max_candidates,
                allow_apo=mode == "relaxed",
                ligand_filter_mode="relaxed" if mode == "relaxed" else "strict",
            )
            selected_count = len(target_install_result.selected_ids)
            if mode == "strict":
                strict_target_count = selected_count
            if mode == "relaxed":
                relaxed_target_count = selected_count
            manifest["target_install"] = {
                "genes": len(genes),
                "selected_count": selected_count,
                "selected_out": str(target_install_result.selected_out),
                "manifest_path": str(target_install_result.manifest_path),
                "installed": len(target_install_result.installed_ids),
                "mode": mode,
                "allow_apo": mode == "relaxed",
                "ligand_filter_mode": "relaxed" if mode == "relaxed" else "strict",
            }
        else:
            manifest["target_install"] = {"status": "skipped", "reason": "no targets to stage"}

    summary = {
        "status": "written",
        "modes": args.mode,
        "results": mode_results,
        "strict_target_count": strict_target_count,
        "relaxed_target_count": relaxed_target_count,
        "commands": {
            "install_targets_enabled": bool(args.install_targets),
            "install_dry_run": bool(args.install_dry_run),
            "max_candidates": args.max_candidates,
            "max_return": args.max_return,
        },
    }
    (out_root / "stage_spd_panel_manifest.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

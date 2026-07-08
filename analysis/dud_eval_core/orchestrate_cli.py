from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd  # type: ignore[import-untyped]

import analysis.dud_eval_core.types as dud_eval_types
from analysis.dud_eval_core.discovery import (
    _is_probable_run_id_dirname,
    _is_under,
    _normalize_pdb_ids,
    _resolve_reranked_scorch_path,
    _resolve_scan_roots,
    derive_target_name,
    make_target_key,
)
from analysis.dud_eval_core.engine import (
    emit_consensus_summary,
    emit_reranked_scorch_summary,
    evaluate_target,
    evaluate_target_consensus,
    evaluate_target_post_docked_reranked_scorch,
)
from analysis.dud_eval_core.labels import _make_placeholder_row
from analysis.dud_eval_core.log import dbg, set_log_level
from analysis.dud_eval_core.orchestrate_discovery import (
    CONTROL_CENTERS_RE,
    CONTROL_REDOCK_RE,
    _CONTROL_PATTERNS_LOGGED,
    _build_control_records,
    _candidate_protein_logs,
    _compute_analysis_root,
    _format_run_label,
    _load_default_cfg,
    _manifest_library_for_target,
    _manifest_protein_entries,
    _manifest_proteins_by_pdb,
    _resolve_run_label,
    select_default_run_id,
)
from analysis.dud_eval_core.orchestrate_metrics import (
    DUD_E_PROTEIN_CLASSES,
    _class_metric_base_columns,
    _class_metric_stat_columns,
    _class_slug,
    _emit_class_aggregate_plots,
    _target_keys_from_summary_df,
    _write_placeholder_metrics,
)
from analysis.dud_eval_core.reporting import (
    _write_pretty_summary,
    _write_pretty_table_noformat,
)
from analysis.dud_eval_core.types import (
    CONSENSUS_CSV_BASENAME,
    CSV_BASENAMES,
    TargetEvaluation,
    TargetSpec,
)
from cli.run_manifest_runtime import _load_manifest, get_manifest_paths
from src.path_router.path_router import expand_variants, make_paths

def main():
    ap = argparse.ArgumentParser(
        description="Atlas VS benchmark evaluator (filename-labeled actives/decoys)."
    )
    ap.add_argument(
        "--docked-root",
        type=str,
        default="docked",
        help="Root folder (or a single target folder) to scan for docking_score_long.csv.",
    )
    ap.add_argument(
        "--out-dir",
        type=str,
        default="analysis/dud_eval",
        help="Output root directory (run_id is appended automatically when provided).",
    )
    ap.add_argument(
        "--post-docked-root",
        type=str,
        default="post_docked",
        help="Root folder for post-docked outputs (consensus reranked with SCORCH).",
    )
    ap.add_argument(
        "--post-docked-basename",
        type=str,
        default="consensus_reranked_scorch.csv",
        help="Filename for post-docked reranked SCORCH CSV (default: consensus_reranked_scorch.csv).",
    )
    try:
        ap.add_argument(
            "--eval-post-docked-scorch",
            action=argparse.BooleanOptionalAction,
            default=True,
            help="Evaluate post-docked consensus reranked with SCORCH (default: on).",
        )
    except Exception:
        ap.add_argument(
            "--eval-post-docked-scorch",
            action="store_true",
            default=True,
            help="Evaluate post-docked consensus reranked with SCORCH (default: on).",
        )
    ap.add_argument(
        "--lig-col",
        type=str,
        default=None,
        help="Column containing ligand filename/path (auto-detected if omitted).",
    )
    ap.add_argument(
        "--score-col",
        type=str,
        default=None,
        help="Score column (lower is better). Auto-detected if omitted; for reranked SCORCH, defaults to final_score for full-library evaluation (override for subset analyses, e.g., --score-col scorch_composite).",
    )
    ap.add_argument(
        "--run-id",
        type=str,
        default=None,
        help="Filter docking_score_long.csv rows to a specific run identifier.",
    )
    ap.add_argument(
        "--pdb-id",
        action="append",
        default=None,
        help="Restrict evaluation to specific PDB IDs (repeatable or comma-separated).",
    )
    ap.add_argument(
        "-valid",
        "--valid",
        dest="valid_only",
        action="store_true",
        default=False,
        help="Only evaluate poses marked valid in the docking_score_long.csv file.",
    )
    ap.add_argument(
        "--valid-col",
        type=str,
        default=None,
        help="Validity column name to use with --valid (auto-detected if omitted).",
    )
    ap.add_argument("--bedroc-alpha", type=float, default=20.0)
    ap.add_argument("--logauc-lambda", type=float, default=1e-3)
    ap.add_argument(
        "--class-min-targets",
        type=int,
        default=1,
        help="Minimum class-member targets required before class metrics are computed (default: 1).",
    )
    try:
        ap.add_argument(
            "--class-placeholder-images",
            action=argparse.BooleanOptionalAction,
            default=False,
            help="When class metrics are unavailable, emit placeholder PNGs (default: off).",
        )
    except Exception:
        ap.add_argument(
            "--class-placeholder-images",
            action="store_true",
            default=False,
            help="When class metrics are unavailable, emit placeholder PNGs.",
        )
    ap.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=("DEBUG", "INFO", "WARN", "ERROR"),
        help="Logging verbosity (default: INFO).",
    )
    ap.add_argument(
        "--target-name-from-pdb",
        action="store_true",
        default=True,
        help="If set, add a target_name column derived from PDB headers.",
    )
    ap.add_argument(
        "--target-name-prefer",
        type=str,
        default="auto",
        choices=("auto", "compnd", "uniprot"),
        help="Preference order when selecting target_name (default: auto).",
    )
    ap.add_argument(
        "--pdb-root",
        type=str,
        default=None,
        help="Optional override root for processed PDB folders (processed_pdbs/<target>/).",
    )
    ap.add_argument(
        "--report-library",
        action="store_true",
        default=True,
        help="If set, attempt to infer library_name from prepped_ligands.",
    )
    ap.add_argument(
        "--prepped-root",
        type=str,
        default=None,
        help="Optional override root for prepped ligands (prepped_ligands/<target>/).",
    )
    ap.add_argument(
        "--emit-control-report",
        action="store_true",
        default=True,
        help="If set, parse protein.log control entries and emit control_redock TSVs.",
    )
    ap.add_argument(
        "--log-root",
        type=str,
        default=None,
        help="Optional override root containing <PDB>/protein.log (default: docked/).",
    )
    try:
        ap.add_argument(
            "--pretty-summary",
            action=argparse.BooleanOptionalAction,
            default=True,
            help="Also write a human-readable aligned text summary (default: on).",
        )
    except Exception:
        ap.add_argument(
            "--pretty-summary",
            dest="pretty_summary",
            action="store_true",
            help="Write a human-readable aligned text summary.",
        )
        ap.add_argument(
            "--no-pretty-summary", dest="pretty_summary", action="store_false"
        )
        ap.set_defaults(pretty_summary=True)

    args = ap.parse_args()

    if getattr(args, "post_docked_basename", None):
        dud_eval_types.RERANKED_SCORCH_BASENAME = args.post_docked_basename

    set_log_level(args.log_level)
    if args.class_min_targets < 1:
        dbg(
            "WARN",
            "class.aggregate",
            f"class_min_targets={args.class_min_targets} invalid; using 1",
        )
        args.class_min_targets = 1
    dbg(
        "DEBUG",
        "args",
        f"log_level={args.log_level} target_name_from_pdb={'ON' if args.target_name_from_pdb else 'OFF'} report_library={'ON' if args.report_library else 'OFF'} control_report={'ON' if args.emit_control_report else 'OFF'} run_id={args.run_id or 'none'}",
    )
    pdb_id_filter = _normalize_pdb_ids(args.pdb_id)
    active_run_id: Optional[str] = args.run_id
    csv_basenames: Tuple[str, ...] = CSV_BASENAMES

    cfg = _load_default_cfg()
    docked_root = Path(args.docked_root)
    out_root = Path(args.out_dir)
    post_docked_root = Path(args.post_docked_root)
    pdb_root_override = Path(args.pdb_root) if args.pdb_root else None
    prepped_root_override = Path(args.prepped_root) if args.prepped_root else None
    log_root_override = Path(args.log_root) if getattr(args, "log_root", None) else None

    analysis_root = _compute_analysis_root(args, cfg)
    dbg(
        "DEBUG",
        "paths",
        f"docked_root={docked_root} out_root={out_root} analysis_root={analysis_root} pdb_root={pdb_root_override or 'none'} prepped_root={prepped_root_override or 'none'}",
    )
    scan_roots = _resolve_scan_roots(docked_root, args.run_id)
    dbg("INFO", "paths", f"scan_roots={[str(p) for p in scan_roots]}")
    manifest_docked_root_primary = scan_roots[0] if scan_roots else docked_root
    legacy_root = docked_root

    if args.run_id:
        os.environ["ATLAS_RUN_ID"] = str(args.run_id)
        if cfg is not None:
            cfg["RUN_ID"] = str(args.run_id)

    def _infer_variant_ph_from_csv_path(
        pdb_id: str, csv_path: Path
    ) -> tuple[Optional[str], Optional[str]]:
        """
        Infer (variant, pH_tag) from the chosen docking_score_long.csv path.

        We treat 'docked/<PDB>/' from path_router as the root and look at
        the extra directory segments between that root and the CSV's parent.
        Any 'APO'/'HOLO' segment becomes the variant, any segment starting
        with 'pH'/'ph' becomes the pH tag.
        """
        try:
            paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
            base = paths.docked_pdb_root()  # docked/<PDB>/
        except Exception:
            base = docked_root / pdb_id

        try:
            rel = csv_path.parent.relative_to(base)
            parts = [p for p in rel.parts if p]
        except Exception:
            parts = []

        variant: Optional[str] = None
        ph_tag: Optional[str] = None
        for part in parts:
            up = part.upper()
            if variant is None and up in {"APO", "HOLO"}:
                variant = up
                continue
            low = part.lower()
            if ph_tag is None and low.startswith("ph"):
                ph_tag = part
        return variant, ph_tag

    def _resolve_docking_csv(pdb_id: str, fallback_dir: Path) -> Optional[Path]:
        candidates_tried: List[str] = []
        picked: Optional[Path] = None

        def _record(path: Path) -> bool:
            nonlocal picked
            exists = path.exists()
            state = "hit" if exists else "miss"
            candidates_tried.append(f"{path} ({state})")
            if exists and picked is None:
                picked = path
                return True
            return False

        def _record_basenames(base_path: Path) -> None:
            """
            Try all CSV_BASENAMES at base_path, e.g. base_path / name,
            in order. The first existing file wins via _record().
            """
            nonlocal picked
            for basename in csv_basenames:
                if picked is not None:
                    break
                candidate = base_path / basename
                _record(candidate)

        def _scan_tree(root: Path) -> None:
            """
            Look for CSV_BASENAMES starting at root.

            Search order:
            1) root/<basename> for each basename in CSV_BASENAMES.
            2) root/<variant>/<basename> (variant is any subdir name).
            3) root/<variant>/<pH>/<basename> (pH is any sub-subdir name).
            We stop as soon as _record() finds a hit.
            """
            nonlocal picked

            # 1) root/<basename>
            if picked is None:
                _record_basenames(root)

            # 2) root/<variant>/<basename> and 3) root/<variant>/<pH>/<basename>
            if picked is not None:
                return

            if not root.is_dir():
                return

            for variant_dir in sorted(root.iterdir()):
                if picked is not None:
                    break
                if not variant_dir.is_dir():
                    continue

                # 2) docked/<PDB>/<VARIANT>/<basename>
                _record_basenames(variant_dir)

                if picked is not None:
                    break

                # 3) docked/<PDB>/<VARIANT>/<pH>/<basename>
                for ph_dir in sorted(variant_dir.iterdir()):
                    if picked is not None:
                        break
                    if not ph_dir.is_dir():
                        continue
                    _record_basenames(ph_dir)

        if not cfg:
            _record_basenames(fallback_dir)
            if picked is None:
                _scan_tree(fallback_dir)
        else:
            try:
                mode = str(cfg.get("APO_HOLO_MODE", "")).strip()
                variants = expand_variants(mode)  # returns [None] | ["APO","HOLO"]
                paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
            except Exception as exc:
                dbg("WARN", "resolve", f"pdb={pdb_id} router_err={exc}")
                _record_basenames(fallback_dir)
                if picked is None:
                    _scan_tree(fallback_dir)
            else:
                docked_root_cfg = paths.docked_pdb_root()
                summary_csv = paths.docking_score_summary_csv()
                long_csv = paths.docking_score_long_csv()

                if docked_root_cfg.exists():
                    # 1) Try path_router's long_csv location, but with both basenames
                    for basename in csv_basenames:
                        if picked is not None:
                            break
                        candidate = long_csv.with_name(basename)
                        _record(candidate)

                    # 2) If still nothing, try variant roots (and their pH subdirs)
                    if picked is None:
                        for variant in variants:
                            variant_root = paths.docked_variant_root(variant)
                            _record_basenames(variant_root)
                            if picked is not None:
                                break
                            if variant_root.is_dir():
                                for ph_dir in sorted(variant_root.iterdir()):
                                    if picked is not None:
                                        break
                                    if ph_dir.is_dir():
                                        _record_basenames(ph_dir)

                    # 3) Final fallback: CLI fallback_dir (includes nested scan)
                    if picked is None:
                        _scan_tree(fallback_dir)

                    # Touch summary path to exercise router (no fallback to summary file for eval)
                    try:
                        summary_csv.exists()
                    except Exception:
                        pass
                else:
                    dbg(
                        "WARN",
                        "resolve",
                        f"pdb={pdb_id} docked_root_missing root={docked_root_cfg}",
                    )
                    _record_basenames(fallback_dir)

        if picked is not None:
            dbg(
                "DEBUG",
                "resolve",
                f"pdb={pdb_id} tried={len(candidates_tried)} candidates={candidates_tried}",
            )
            dbg("INFO", "resolve", f"pdb={pdb_id} picked={picked}")

            return picked

        dbg(
            "WARN",
            "resolve",
            f"pdb={pdb_id} no_csv_found tried={candidates_tried or ['<none>']} search_root={fallback_dir}",
        )
        return None

    def _resolve_docking_csv_for_manifest(
        pdb_id: str,
        variant_hint: Optional[str],
        ph_hint: Optional[str],
        docked_root: Path,
        cfg: Dict,
        fallback_root: Optional[Path] = None,
    ) -> Optional[Path]:
        """
        Variant/pH-aware wrapper around _resolve_docking_csv.

        Strategy:
          1) If cfg and variant_hint/ph_hint are provided, try the exact
             docked/<PDB>/<VARIANT>/<PH>/<CSV_BASENAME> locations first.
          2) If that fails, fall back to _resolve_docking_csv(pdb_id, docked_root / pdb_id).
          3) If still missing and fallback_root provided, try legacy fallback_root/pdb_id.
        """
        try:
            paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
            base_root = paths.docked_pdb_root()
        except Exception:
            paths = None
            base_root = docked_root / pdb_id

        if variant_hint or ph_hint:
            var = (variant_hint or "").strip()
            ph = (ph_hint or "").strip()
            candidate_roots: List[Path] = []
            if var:
                try:
                    variant_root = (
                        paths.docked_variant_root(var) if paths else base_root / var
                    )
                except Exception:
                    variant_root = base_root / var
                if ph:
                    candidate_roots.append(variant_root / ph)
                else:
                    candidate_roots.append(variant_root)
            elif ph:
                candidate_roots.append(base_root / ph)

            for root in candidate_roots:
                for basename in csv_basenames:
                    candidate = root / basename
                    if candidate.exists():
                        dbg(
                            "INFO",
                            "resolve",
                            f"pdb={pdb_id} picked={candidate} via=manifest_hint",
                        )
                        return candidate

        picked = _resolve_docking_csv(pdb_id, docked_root / pdb_id)
        if picked is None and fallback_root and fallback_root != docked_root:
            picked = _resolve_docking_csv(pdb_id, fallback_root / pdb_id)
        return picked

    def _discover_docking_csvs(
        pdb_id: str,
        docked_root: Path,
        cfg: Optional[dict],
        csv_basenames: Tuple[str, ...],
    ) -> List[TargetSpec]:
        """
        Return all docking CSVs for a pdb_id across variants/pH subfolders.
        """
        try:
            paths = (
                make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
                if cfg
                else None
            )
        except Exception:
            paths = None

        # Primary discovery root must come from the explicit scan tree passed
        # into this function (derived from --docked-root / run scan roots).
        scan_base_root = docked_root / pdb_id
        if docked_root.name.upper() == pdb_id.upper():
            scan_base_root = docked_root

        # Secondary fallback root from path-router/config. Keep it as optional
        # backup only, because configs may point at a different filesystem than
        # the CLI-provided docked root (e.g., /work2 config vs /scratch outputs).
        router_base_root: Optional[Path] = None
        if paths:
            try:
                router_base_root = paths.docked_pdb_root()
            except Exception:
                router_base_root = None

        candidate_roots: List[Path] = [scan_base_root]
        if router_base_root is not None:
            try:
                same_root = router_base_root.resolve() == scan_base_root.resolve()
            except Exception:
                same_root = router_base_root == scan_base_root
            if not same_root:
                candidate_roots.append(router_base_root)

        seen_paths: Set[str] = set()
        ranked: List[Tuple[TargetSpec, int]] = []
        basename_rank = {name: idx for idx, name in enumerate(csv_basenames)}

        def _add(
            candidate: Path, variant: Optional[str], ph_tag: Optional[str], source: str
        ) -> None:
            if not candidate.exists():
                return
            try:
                key = str(candidate.resolve())
            except Exception:
                key = str(candidate)
            if key in seen_paths:
                return
            seen_paths.add(key)
            variant_norm = variant.upper() if variant else None
            rank = basename_rank.get(candidate.name, len(csv_basenames))
            ranked.append(
                (
                    TargetSpec(
                        target_key=make_target_key(pdb_id, variant_norm, ph_tag),
                        pdb_id=pdb_id,
                        variant=variant_norm,
                        ph_tag=ph_tag,
                        csv_path=candidate,
                        source=source,
                    ),
                    rank,
                )
            )

        for root_idx, base_root in enumerate(candidate_roots):
            source_tag = "scan" if root_idx == 0 else "router_fallback"

            # Legacy root: docked/<PDB>/<basename>
            for basename in csv_basenames:
                _add(base_root / basename, None, None, source_tag)

            if base_root.is_dir():
                for ph_dir in sorted(base_root.iterdir()):
                    if not ph_dir.is_dir():
                        continue
                    ph_label = ph_dir.name
                    for basename in csv_basenames:
                        _add(ph_dir / basename, None, ph_label, source_tag)

            for variant in ("APO", "HOLO"):
                variant_dir = base_root / variant
                for basename in csv_basenames:
                    _add(variant_dir / basename, variant, None, source_tag)
                if variant_dir.is_dir():
                    for ph_dir in sorted(variant_dir.iterdir()):
                        if not ph_dir.is_dir():
                            continue
                        ph_label = ph_dir.name
                        for basename in csv_basenames:
                            _add(ph_dir / basename, variant, ph_label, source_tag)

        ranked.sort(
            key=lambda pair: (
                pair[0].pdb_id,
                pair[0].variant or "",
                pair[0].ph_tag or "",
                pair[1],
                str(pair[0].csv_path),
            )
        )
        specs = [spec for spec, _rank in ranked]
        pairs = [
            f"{spec.variant or 'legacy'}|{spec.ph_tag or 'base'}" for spec in specs
        ]
        variant_label = ";".join(pairs) if pairs else "<none>"
        dbg(
            "INFO",
            "discover",
            f"pdb={pdb_id} targets={len(specs)} variants={variant_label}",
        )
        return specs

    manifest_data: dict | None = None
    manifest_entries: List[Dict[str, str]] = []
    manifest_proteins: Dict[str, Dict[str, str]] = {}
    manifest_driven = False
    if cfg and active_run_id:
        if (
            len(scan_roots) == 1
            and scan_roots[0].is_dir()
            and _is_probable_run_id_dirname(scan_roots[0].name)
        ):
            pass
        try:
            _manifest_dir, manifest_path = get_manifest_paths(cfg, str(active_run_id))
            manifest_data = _load_manifest(manifest_path) or {}
            manifest_entries = _manifest_protein_entries(manifest_data or {})
            manifest_proteins = _manifest_proteins_by_pdb(manifest_data or {})
            manifest_driven = bool(active_run_id and manifest_entries)
            test_mode_enable = ""
            cmd_block = (
                manifest_data.get("command") if isinstance(manifest_data, dict) else {}
            )
            if isinstance(cmd_block, dict):
                test_mode_enable = (
                    str(cmd_block.get("TEST_MODE_ENABLE", "")).strip().lower()
                )
            if test_mode_enable == "dud":
                csv_basenames = ("docking_score_long.csv", "dud_docking_score_long.csv")
                dbg(
                    "INFO",
                    "manifest",
                    f"run_id={active_run_id} test_mode=dud csv_priority={';'.join(csv_basenames)}",
                )
            dbg(
                "INFO",
                "manifest",
                f"run_id={active_run_id} manifest={manifest_path} proteins={len(manifest_entries)}",
            )
        except Exception as exc:
            dbg("WARN", "manifest", f"run_id={active_run_id} load_err={exc}")
            manifest_data = {}
            manifest_entries = []
            manifest_proteins = {}
            manifest_driven = False
    else:
        dbg(
            "DEBUG",
            "manifest",
            "pre-discover: no cfg or active_run_id; manifest lookup disabled",
        )

    targets: List[TargetSpec] = []

    if manifest_driven:
        for entry in manifest_entries:
            pdb_id = str(entry.get("pdb_id", "")).strip()
            if not pdb_id:
                continue
            variant_hint = str(entry.get("variant", "") or "").strip()
            ph_hint = str(entry.get("ph", "") or "").strip()
            variant_upper = variant_hint.upper() if variant_hint else None
            csv_path = _resolve_docking_csv_for_manifest(
                pdb_id=pdb_id,
                variant_hint=variant_hint or None,
                ph_hint=ph_hint or None,
                docked_root=manifest_docked_root_primary,
                cfg=cfg,
                fallback_root=legacy_root,
            )
            target_key = make_target_key(pdb_id, variant_upper, ph_hint or None)
            targets.append(
                TargetSpec(
                    target_key=target_key,
                    pdb_id=pdb_id,
                    variant=variant_upper,
                    ph_tag=ph_hint or None,
                    csv_path=csv_path,
                    source="manifest",
                )
            )
            if csv_path is None:
                dbg(
                    "WARN",
                    "discover",
                    f"pdb={pdb_id} variant={variant_upper or 'base'} ph={ph_hint or 'base'} reason=no_docked_csv_for_manifest_entry",
                )
    elif docked_root.is_dir():
        # Iterate over candidate scan roots (run-scoped and legacy).
        for root in scan_roots:
            if not root.is_dir():
                # Treat as single-target fallback
                pdb_id = root.name
                targets.extend(_discover_docking_csvs(pdb_id, root, cfg, csv_basenames))
                continue

            try:
                root_is_docked_root = root.resolve() == docked_root.resolve()
            except Exception:
                root_is_docked_root = root == docked_root
            any_subdirs = False
            any_added = False
            for sub in sorted(root.iterdir()):
                if not sub.is_dir():
                    continue
                any_subdirs = True
                pdb_id = sub.name

                # When scanning the top-level docked root, avoid treating
                # run-id-like folders as PDB targets.
                if root_is_docked_root and _is_probable_run_id_dirname(pdb_id):
                    dbg(
                        "DEBUG",
                        "discover.skip",
                        f"pdb={pdb_id} reason=run_like_subdir_under_docked_root",
                    )
                    continue

                pdb_path = Path("input_pdbs") / f"{pdb_id}.pdb"
                if not pdb_path.exists():
                    dbg(
                        "WARN",
                        "discover.fallback",
                        f"pdb={pdb_id} reason=no_input_pdb path={pdb_path} action=scan_outputs_anyway",
                    )

                discovered = _discover_docking_csvs(pdb_id, root, cfg, csv_basenames)
                targets.extend(discovered)
                if discovered:
                    any_added = True

            if not any_subdirs:
                pdb_id = root.name
                targets.extend(_discover_docking_csvs(pdb_id, root, cfg, csv_basenames))
            elif not any_added:
                dbg(
                    "WARN",
                    "discover",
                    f"root={root} reason=no_targets_in_subdirs",
                )
    else:
        pdb_id = docked_root.name
        targets.extend(_discover_docking_csvs(pdb_id, docked_root, cfg, csv_basenames))

    unique_targets: List[TargetSpec] = []
    seen_keys: Set[str] = set()
    for spec in targets:
        if spec.target_key in seen_keys:
            dbg(
                "WARN",
                "discover",
                f"target_key={spec.target_key} csv={spec.csv_path} source={spec.source} action=skip_duplicate",
            )
            continue
        seen_keys.add(spec.target_key)
        unique_targets.append(spec)
    targets = unique_targets

    if pdb_id_filter:
        filter_set = set(pdb_id_filter)
        before = len(targets)
        targets = [t for t in targets if t.pdb_id.upper() in filter_set]
        dbg(
            "INFO",
            "dud-eval",
            f"filter pdb_id={sorted(filter_set)} kept={len(targets)}/{before}",
        )

    if not targets:
        if pdb_id_filter:
            dbg(
                "ERROR",
                "discover",
                f"no targets after pdb_id filter pdb_id={pdb_id_filter}",
            )
        else:
            dbg(
                "ERROR",
                "discover",
                f"no docking_score_long.csv / dud_docking_score_long.csv under {docked_root}",
            )
        raise SystemExit(2)

    dbg("INFO", "discover", f"targets={len(targets)} root={docked_root}")

    rows: List[pd.Series] = []
    consensus_rows: List[pd.Series] = []
    reranked_rows: List[pd.Series] = []
    target_eval_results: Dict[str, TargetEvaluation] = {}
    consensus_eval_results: Dict[str, TargetEvaluation] = {}
    reranked_eval_results: Dict[str, TargetEvaluation] = {}
    target_specs_by_key: Dict[str, TargetSpec] = {spec.target_key: spec for spec in targets}
    csv_paths_by_target: Dict[str, Path] = {}
    for spec in targets:
        if spec.csv_path:
            csv_paths_by_target[spec.target_key] = spec.csv_path
    targets_with_csv = [t for t in targets if t.csv_path]

    if not active_run_id:
        try:
            active_run_id = select_default_run_id(
                targets_with_csv, csv_paths_by_target, args.lig_col, args.score_col
            )
        except Exception as _exc:
            active_run_id = None
            dbg("WARN", "run", f"auto_select_failed err={_exc}")
    dbg(
        "INFO",
        "run",
        f"active={active_run_id or '(none)'} source={'CLI' if args.run_id else 'auto'}",
    )
    run_dir = Path(args.docked_root) / str(active_run_id) if active_run_id else None
    if active_run_id:
        analysis_root = _compute_analysis_root(args, cfg, run_id_override=active_run_id)
        dbg(
            "INFO",
            "paths",
            f"analysis_root_updated run_id={active_run_id} path={analysis_root}",
        )
    (analysis_root / "consensus").mkdir(parents=True, exist_ok=True)
    (analysis_root / "post_docked").mkdir(parents=True, exist_ok=True)
    new_layout_root_exists = bool(run_dir and run_dir.exists())
    dbg(
        "INFO",
        "run.layout",
        f"run_id={active_run_id or '(none)'} new_layout_root_exists={new_layout_root_exists}",
    )
    run_label_for_files = _format_run_label(active_run_id)

    if manifest_data is None and cfg and active_run_id:
        try:
            _manifest_dir, manifest_path = get_manifest_paths(cfg, str(active_run_id))
            manifest_data = _load_manifest(manifest_path) or {}
            manifest_entries = _manifest_protein_entries(manifest_data or {})
            manifest_proteins = _manifest_proteins_by_pdb(manifest_data or {})
            manifest_driven = manifest_driven or bool(
                active_run_id and manifest_entries
            )
            test_mode_enable = ""
            cmd_block = (
                manifest_data.get("command") if isinstance(manifest_data, dict) else {}
            )
            if isinstance(cmd_block, dict):
                test_mode_enable = (
                    str(cmd_block.get("TEST_MODE_ENABLE", "")).strip().lower()
                )
            if test_mode_enable == "dud":
                csv_basenames = ("docking_score_long.csv", "dud_docking_score_long.csv")
                dbg(
                    "INFO",
                    "manifest",
                    f"run_id={active_run_id} test_mode=dud csv_priority={';'.join(csv_basenames)}",
                )
            dbg(
                "INFO",
                "manifest",
                f"run_id={active_run_id} manifest={manifest_path} proteins={len(manifest_entries)}",
            )
        except Exception as exc:
            dbg("WARN", "manifest", f"run_id={active_run_id} load_err={exc}")
            manifest_data = {}
            manifest_entries = []
            manifest_proteins = {}
            manifest_driven = bool(active_run_id and manifest_proteins)
    elif manifest_data is None:
        dbg("DEBUG", "manifest", "no cfg or active_run_id; manifest lookup disabled")

    if not manifest_entries:
        manifest_entries = _manifest_protein_entries(manifest_data or {})
    if not manifest_proteins:
        manifest_proteins = _manifest_proteins_by_pdb(manifest_data or {})
    manifest_driven = manifest_driven or bool(active_run_id and manifest_entries)

    for spec in targets:
        evaluated: Optional[TargetEvaluation] = None
        status_reason = "missing_csv"
        out_dir = analysis_root / spec.target_key
        consensus_out_dir = analysis_root / "consensus" / spec.target_key
        reranked_out_dir = (
            analysis_root / "post_docked" / "consensus_reranked_scorch" / spec.target_key
        )
        if spec.csv_path:
            evaluated = evaluate_target(
                pdb_id=spec.pdb_id,
                csv_path=spec.csv_path,
                out_dir=out_dir,
                lig_col_cli=args.lig_col,
                score_col_cli=args.score_col,
                bedroc_alpha=args.bedroc_alpha,
                logauc_lambda=args.logauc_lambda,
                run_id=active_run_id,
                valid_only=args.valid_only,
                valid_col_cli=args.valid_col,
                new_layout_active=bool(
                    run_dir and spec.csv_path and _is_under(spec.csv_path, run_dir)
                ),
                layout_label="new_run"
                if (run_dir and spec.csv_path and _is_under(spec.csv_path, run_dir))
                else "legacy",
            )
            if evaluated is not None:
                status_reason = evaluated.status_reason or "ok"
            else:
                status_reason = "read_error"
        consensus_candidate: Optional[Path] = None
        if spec.csv_path:
            candidate = spec.csv_path.parent / CONSENSUS_CSV_BASENAME
            if candidate.exists():
                consensus_candidate = candidate
        if consensus_candidate is not None:
            cons_eval = evaluate_target_consensus(
                pdb_id=spec.pdb_id,
                csv_path=consensus_candidate,
                out_dir=consensus_out_dir,
                lig_col_cli=args.lig_col,
                score_col_cli=args.score_col,
                bedroc_alpha=args.bedroc_alpha,
                logauc_lambda=args.logauc_lambda,
                run_id=active_run_id,
                variant_hint=spec.variant,
                ph_tag=spec.ph_tag,
                new_layout_active=bool(
                    run_dir
                    and consensus_candidate
                    and _is_under(consensus_candidate, run_dir)
                ),
                layout_label="new_run"
                if (
                    run_dir
                    and consensus_candidate
                    and _is_under(consensus_candidate, run_dir)
                )
                else "legacy",
            )
            if cons_eval is not None:
                consensus_eval_results[spec.target_key] = cons_eval
            if cons_eval is not None and cons_eval.metrics is not None:
                cons_metrics = cons_eval.metrics.copy()
                for col in ("variant", "pH"):
                    if col in cons_metrics.index:
                        cons_metrics.drop(index=col, inplace=True)
                cons_metrics["pdb_id"] = spec.pdb_id
                cons_metrics["variant"] = spec.variant or ""
                cons_metrics["pH"] = spec.ph_tag or ""
                if "status_reason" not in cons_metrics:
                    cons_metrics["status_reason"] = cons_eval.status_reason or "ok"
                consensus_rows.append(cons_metrics)
            elif cons_eval is not None:
                placeholder_cons = _make_placeholder_row(
                    spec.pdb_id,
                    spec.variant,
                    spec.ph_tag,
                    active_run_id,
                    args.bedroc_alpha,
                    cons_eval.status_reason or "empty_after_filter",
                )
                consensus_rows.append(placeholder_cons)
                _write_placeholder_metrics(consensus_out_dir, placeholder_cons)
            else:
                dbg(
                    "WARN",
                    "consensus.eval",
                    f"pdb={spec.pdb_id} path={consensus_candidate} reason=evaluation_failed",
                )
                placeholder_cons = _make_placeholder_row(
                    spec.pdb_id,
                    spec.variant,
                    spec.ph_tag,
                    active_run_id,
                    args.bedroc_alpha,
                    "read_error",
                )
                consensus_rows.append(placeholder_cons)
                _write_placeholder_metrics(consensus_out_dir, placeholder_cons)
        else:
            placeholder_cons = _make_placeholder_row(
                spec.pdb_id,
                spec.variant,
                spec.ph_tag,
                active_run_id,
                args.bedroc_alpha,
                "missing_consensus_csv",
            )
            consensus_rows.append(placeholder_cons)
            _write_placeholder_metrics(consensus_out_dir, placeholder_cons)
        reranked_candidate: Optional[Path] = None
        if getattr(args, "eval_post_docked_scorch", True):
            reranked_candidate = _resolve_reranked_scorch_path(
                spec, docked_root, post_docked_root, active_run_id, run_dir
            )
        if reranked_candidate is not None:
            post_run_root = (
                post_docked_root / str(active_run_id) if active_run_id else None
            )
            new_layout_reranked = bool(
                post_run_root and _is_under(reranked_candidate, post_run_root)
            )
            reranked_eval = evaluate_target_post_docked_reranked_scorch(
                pdb_id=spec.pdb_id,
                csv_path=reranked_candidate,
                out_dir=reranked_out_dir,
                lig_col_cli=args.lig_col,
                score_col_cli=args.score_col,
                bedroc_alpha=args.bedroc_alpha,
                logauc_lambda=args.logauc_lambda,
                run_id=active_run_id,
                variant_hint=spec.variant,
                ph_tag=spec.ph_tag,
                new_layout_active=new_layout_reranked,
                layout_label="post_docked_new"
                if new_layout_reranked
                else "post_docked_legacy",
            )
            if reranked_eval is not None:
                reranked_eval_results[spec.target_key] = reranked_eval
            if reranked_eval is not None and reranked_eval.metrics is not None:
                rer_metrics = reranked_eval.metrics.copy()
                for col in ("variant", "pH"):
                    if col in rer_metrics.index:
                        rer_metrics.drop(index=col, inplace=True)
                rer_metrics["pdb_id"] = spec.pdb_id
                rer_metrics["variant"] = spec.variant or ""
                rer_metrics["pH"] = spec.ph_tag or ""
                rer_metrics["run_id"] = active_run_id or ""
                if "status_reason" not in rer_metrics:
                    rer_metrics["status_reason"] = reranked_eval.status_reason or "ok"
                reranked_rows.append(rer_metrics)
            elif reranked_eval is not None:
                placeholder_reranked = _make_placeholder_row(
                    spec.pdb_id,
                    spec.variant,
                    spec.ph_tag,
                    active_run_id,
                    args.bedroc_alpha,
                    reranked_eval.status_reason or "empty_after_filter",
                )
                reranked_rows.append(placeholder_reranked)
                _write_placeholder_metrics(reranked_out_dir, placeholder_reranked)
            else:
                placeholder_reranked = _make_placeholder_row(
                    spec.pdb_id,
                    spec.variant,
                    spec.ph_tag,
                    active_run_id,
                    args.bedroc_alpha,
                    "read_error",
                )
                reranked_rows.append(placeholder_reranked)
                _write_placeholder_metrics(reranked_out_dir, placeholder_reranked)
        elif getattr(args, "eval_post_docked_scorch", True):
            placeholder_reranked = _make_placeholder_row(
                spec.pdb_id,
                spec.variant,
                spec.ph_tag,
                active_run_id,
                args.bedroc_alpha,
                "missing_reranked_scorch_csv",
            )
            reranked_rows.append(placeholder_reranked)
            _write_placeholder_metrics(reranked_out_dir, placeholder_reranked)
        if evaluated is not None:
            target_eval_results[spec.target_key] = evaluated
            if evaluated.metrics is not None:
                metrics = evaluated.metrics.copy()
                for col in ("variant", "pH"):
                    if col in metrics.index:
                        metrics.drop(index=col, inplace=True)
                metrics["pdb_id"] = spec.pdb_id
                metrics["variant"] = spec.variant or ""
                metrics["pH"] = spec.ph_tag or ""
                if "status_reason" not in metrics:
                    metrics["status_reason"] = status_reason
                rows.append(metrics)
            else:
                placeholder = _make_placeholder_row(
                    spec.pdb_id,
                    spec.variant,
                    spec.ph_tag,
                    active_run_id,
                    args.bedroc_alpha,
                    status_reason,
                )
                rows.append(placeholder)
                _write_placeholder_metrics(out_dir, placeholder)
        else:
            placeholder = _make_placeholder_row(
                spec.pdb_id,
                spec.variant,
                spec.ph_tag,
                active_run_id,
                args.bedroc_alpha,
                status_reason,
            )
            rows.append(placeholder)
            _write_placeholder_metrics(out_dir, placeholder)

    if not rows:
        dbg("ERROR", "metrics", "no targets produced evaluable rows")
        raise SystemExit(3)

    df_all = pd.DataFrame(rows)
    df_consensus_all = (
        pd.DataFrame(consensus_rows) if consensus_rows else pd.DataFrame()
    )
    df_reranked_all = pd.DataFrame(reranked_rows) if reranked_rows else pd.DataFrame()
    for col in ("variant", "pH"):
        if col not in df_all.columns:
            df_all[col] = ""
    has_variant = (
        "variant" in df_all.columns
        and df_all["variant"].astype(str).str.strip().ne("").any()
    )
    has_ph = (
        "pH" in df_all.columns and df_all["pH"].astype(str).str.strip().ne("").any()
    )
    if not has_variant and "variant" in df_all.columns:
        df_all = df_all.drop(columns=["variant"])
    if not has_ph and "pH" in df_all.columns:
        df_all = df_all.drop(columns=["pH"])
    sort_cols = [c for c in ("pdb_id", "variant", "pH") if c in df_all.columns]
    df_all = df_all.sort_values(sort_cols or ["pdb_id"])

    if not df_consensus_all.empty:
        for col in ("variant", "pH"):
            if col not in df_consensus_all.columns:
                df_consensus_all[col] = ""
        cons_has_variant = (
            "variant" in df_consensus_all.columns
            and df_consensus_all["variant"].astype(str).str.strip().ne("").any()
        )
        cons_has_ph = (
            "pH" in df_consensus_all.columns
            and df_consensus_all["pH"].astype(str).str.strip().ne("").any()
        )
        if not cons_has_variant and "variant" in df_consensus_all.columns:
            df_consensus_all = df_consensus_all.drop(columns=["variant"])
        if not cons_has_ph and "pH" in df_consensus_all.columns:
            df_consensus_all = df_consensus_all.drop(columns=["pH"])
        sort_cols_cons = [
            c for c in ("pdb_id", "variant", "pH") if c in df_consensus_all.columns
        ]
        df_consensus_all = df_consensus_all.sort_values(sort_cols_cons or ["pdb_id"])

    if manifest_driven:
        expected_keys = [
            make_target_key(
                str(entry.get("pdb_id", "")).strip(),
                (entry.get("variant") or "").strip() or None,
                (entry.get("ph") or "").strip() or None,
            )
            for entry in manifest_entries
            if str(entry.get("pdb_id", "")).strip()
        ]
        observed_keys: List[str] = []
        for _, row in df_all.iterrows():
            variant_val = ""
            ph_val = ""
            if "variant" in df_all.columns:
                raw_var = row.get("variant", "")
                variant_val = "" if pd.isna(raw_var) else str(raw_var)
            if "pH" in df_all.columns:
                raw_ph = row.get("pH", "")
                ph_val = "" if pd.isna(raw_ph) else str(raw_ph)
            observed_keys.append(
                make_target_key(str(row.get("pdb_id", "")), variant_val, ph_val)
            )
        missing = sorted(set(expected_keys) - set(observed_keys))
        extra = sorted(set(observed_keys) - set(expected_keys))
        dbg(
            "INFO",
            "manifest.coverage",
            f"expected={len(expected_keys)} observed={len(observed_keys)} missing={len(missing)} extra={len(extra)}",
        )
        if missing:
            suffix = "..." if len(missing) > 5 else ""
            dbg("WARN", "manifest.coverage", f"missing_keys={missing[:5]}{suffix}")

    annotate_targets = list(df_all["pdb_id"].astype(str).unique())
    if not df_consensus_all.empty:
        for pdb_id in df_consensus_all["pdb_id"].astype(str).unique():
            if pdb_id not in annotate_targets:
                annotate_targets.append(pdb_id)

    names: Dict[str, str] = {p: "" for p in annotate_targets}
    if args.target_name_from_pdb:
        for pdb_id in annotate_targets:
            names[pdb_id] = derive_target_name(
                pdb_id,
                prefer=args.target_name_prefer,
                pdb_root_override=pdb_root_override,
                cfg=cfg,
            )
        if not df_all.empty:
            df_all["target_name"] = df_all["pdb_id"].map(names).fillna("")
    else:
        dbg("DEBUG", "target", "target_name_from_pdb=OFF")

    if args.report_library:
        if manifest_data:
            df_all["library_name"] = df_all.apply(
                lambda r: _manifest_library_for_target(
                    manifest=manifest_data,
                    pdb_id=str(r.get("pdb_id", "")),
                    variant_label=(
                        ""
                        if "variant" not in df_all.columns
                        else (
                            ""
                            if pd.isna(r.get("variant", ""))
                            else str(r.get("variant", ""))
                        )
                    ),
                    ph_tag=(
                        ""
                        if "pH" not in df_all.columns
                        else ("" if pd.isna(r.get("pH", "")) else str(r.get("pH", "")))
                    ),
                ),
                axis=1,
            )
            df_all["library_name"] = df_all["library_name"].fillna("")
        else:
            df_all["library_name"] = ""
    else:
        dbg("DEBUG", "library", "report_library=OFF")
        if "library_name" not in df_all.columns:
            df_all["library_name"] = ""

    if not df_all.empty and "target_name" not in df_all.columns:
        df_all["target_name"] = df_all["pdb_id"].map(names).fillna("")
    if "library_name" not in df_all.columns:
        df_all["library_name"] = ""

    if not df_consensus_all.empty:
        if args.report_library:
            if manifest_data:
                df_consensus_all["library_name"] = df_consensus_all.apply(
                    lambda r: _manifest_library_for_target(
                        manifest=manifest_data,
                        pdb_id=str(r.get("pdb_id", "")),
                        variant_label=(
                            ""
                            if "variant" not in df_consensus_all.columns
                            else (
                                ""
                                if pd.isna(r.get("variant", ""))
                                else str(r.get("variant", ""))
                            )
                        ),
                        ph_tag=(
                            ""
                            if "pH" not in df_consensus_all.columns
                            else (
                                "" if pd.isna(r.get("pH", "")) else str(r.get("pH", ""))
                            )
                        ),
                    ),
                    axis=1,
                )
                df_consensus_all["library_name"] = df_consensus_all[
                    "library_name"
                ].fillna("")
            else:
                df_consensus_all["library_name"] = df_consensus_all.get(
                    "library_name", pd.Series("", index=df_consensus_all.index)
                )
        else:
            if "library_name" not in df_consensus_all.columns:
                df_consensus_all["library_name"] = ""
        if "target_name" not in df_consensus_all.columns:
            df_consensus_all["target_name"] = (
                df_consensus_all["pdb_id"].map(names).fillna("")
            )

    if not df_reranked_all.empty:
        for col in ("variant", "pH"):
            if col not in df_reranked_all.columns:
                df_reranked_all[col] = ""
        reranked_has_variant = (
            "variant" in df_reranked_all.columns
            and df_reranked_all["variant"].astype(str).str.strip().ne("").any()
        )
        reranked_has_ph = (
            "pH" in df_reranked_all.columns
            and df_reranked_all["pH"].astype(str).str.strip().ne("").any()
        )
        if not reranked_has_variant and "variant" in df_reranked_all.columns:
            df_reranked_all = df_reranked_all.drop(columns=["variant"])
        if not reranked_has_ph and "pH" in df_reranked_all.columns:
            df_reranked_all = df_reranked_all.drop(columns=["pH"])
        sort_cols_reranked = [
            c for c in ("pdb_id", "variant", "pH") if c in df_reranked_all.columns
        ]
        df_reranked_all = df_reranked_all.sort_values(sort_cols_reranked or ["pdb_id"])

        if args.report_library:
            if manifest_data:
                df_reranked_all["library_name"] = df_reranked_all.apply(
                    lambda r: _manifest_library_for_target(
                        manifest=manifest_data,
                        pdb_id=str(r.get("pdb_id", "")),
                        variant_label=(
                            ""
                            if "variant" not in df_reranked_all.columns
                            else (
                                ""
                                if pd.isna(r.get("variant", ""))
                                else str(r.get("variant", ""))
                            )
                        ),
                        ph_tag=(
                            ""
                            if "pH" not in df_reranked_all.columns
                            else (
                                "" if pd.isna(r.get("pH", "")) else str(r.get("pH", ""))
                            )
                        ),
                    ),
                    axis=1,
                )
                df_reranked_all["library_name"] = df_reranked_all[
                    "library_name"
                ].fillna("")
            else:
                df_reranked_all["library_name"] = df_reranked_all.get(
                    "library_name", pd.Series("", index=df_reranked_all.index)
                )
        else:
            if "library_name" not in df_reranked_all.columns:
                df_reranked_all["library_name"] = ""
        if "target_name" not in df_reranked_all.columns:
            df_reranked_all["target_name"] = (
                df_reranked_all["pdb_id"].map(names).fillna("")
            )

    control_targets = [t for t in targets if t.target_key in target_eval_results]

    df = df_all.copy()
    df_consensus = df_consensus_all.copy()

    analysis_root.mkdir(parents=True, exist_ok=True)

    # --- Exclude rows where the inferred library is FDA or the PDB-specific library ---
    # Normalize library names and PDB IDs for robust comparison
    lib_norm = df["library_name"].astype(str).str.strip()
    pdb_norm = df["pdb_id"].astype(str).str.strip()
    # FDA library (accept both "fda_library" and "fda")
    mask_fda = lib_norm.str.lower().isin({"fda_library", "fda"})
    # Per-PDB library: library name equals the PDB code (case-insensitive)
    mask_pdb = lib_norm.str.upper() == pdb_norm.str.upper()
    exclude_mask = mask_fda | mask_pdb
    # Emit a visible exclusions report with reasons
    if exclude_mask.any():
        excluded = df.loc[exclude_mask, ["pdb_id", "library_name"]].copy()
        # Reason per row
        excluded["reason"] = np.where(
            mask_fda.loc[exclude_mask],
            "library name = fda library",
            "library name = pdb",
        )
        excl_path = analysis_root / "excluded.tsv"
        excluded.to_csv(excl_path, sep="\t", index=False)
        dbg("INFO", "exclude", f"excluded={len(excluded)} out={excl_path}")
        # Also print each line for quick visibility
        for _, r in excluded.iterrows():
            print(
                f"[exclude] pdb={r['pdb_id']} library={r['library_name']} reason={r['reason']}"
            )
    else:
        dbg("DEBUG", "exclude", "excluded=0")

    # Keep only the non-excluded rows for downstream metrics/summary
    df = df.loc[~exclude_mask].copy()
    if not df_consensus.empty:
        cons_lib_norm = df_consensus["library_name"].astype(str).str.strip()
        cons_pdb_norm = df_consensus["pdb_id"].astype(str).str.strip()
        cons_mask_fda = cons_lib_norm.str.lower().isin({"fda_library", "fda"})
        cons_mask_pdb = cons_lib_norm.str.upper() == cons_pdb_norm.str.upper()
        cons_exclude_mask = cons_mask_fda | cons_mask_pdb
        dropped_consensus = int(cons_exclude_mask.sum())
        if dropped_consensus:
            dbg(
                "INFO",
                "consensus.exclude",
                f"excluded={dropped_consensus} reason=library_filter",
            )
        df_consensus = df_consensus.loc[~cons_exclude_mask].copy()
        sort_cols_cons = [
            c for c in ("pdb_id", "variant", "pH") if c in df_consensus.columns
        ]
        df_consensus = df_consensus.sort_values(sort_cols_cons or ["pdb_id"])
    if not df_reranked_all.empty:
        rer_lib_norm = df_reranked_all["library_name"].astype(str).str.strip()
        rer_pdb_norm = df_reranked_all["pdb_id"].astype(str).str.strip()
        rer_mask_fda = rer_lib_norm.str.lower().isin({"fda_library", "fda"})
        rer_mask_pdb = rer_lib_norm.str.upper() == rer_pdb_norm.str.upper()
        rer_exclude_mask = rer_mask_fda | rer_mask_pdb
        if rer_exclude_mask.any():
            dbg(
                "INFO",
                "reranked.exclude",
                f"excluded={int(rer_exclude_mask.sum())} reason=library_filter",
            )
        df_reranked_all = df_reranked_all.loc[~rer_exclude_mask].copy()
        sort_cols_rer = [
            c for c in ("pdb_id", "variant", "pH") if c in df_reranked_all.columns
        ]
        df_reranked_all = df_reranked_all.sort_values(sort_cols_rer or ["pdb_id"])

    docking_target_keys = _target_keys_from_summary_df(df)
    consensus_target_keys = _target_keys_from_summary_df(df_consensus)
    reranked_target_keys = _target_keys_from_summary_df(df_reranked_all)

    docking_class_results, class_generation_rows = _emit_class_aggregate_plots(
        analysis_root=analysis_root,
        run_id=active_run_id,
        target_specs_by_key=target_specs_by_key,
        eval_results=target_eval_results,
        included_target_keys=docking_target_keys,
        mode_slug="docking",
        score_high_is_better=False,
        hist_xlabel="Best docking score (lower = better)",
        title_suffix="",
        mode_label="docking",
        bedroc_alpha=args.bedroc_alpha,
        logauc_lambda=args.logauc_lambda,
        class_min_targets=args.class_min_targets,
        write_placeholder_images=args.class_placeholder_images,
    )
    _, consensus_class_rows = _emit_class_aggregate_plots(
        analysis_root=analysis_root,
        run_id=active_run_id,
        target_specs_by_key=target_specs_by_key,
        eval_results=consensus_eval_results,
        included_target_keys=consensus_target_keys,
        mode_slug="consensus",
        score_high_is_better=True,
        hist_xlabel="Consensus score (higher = better; negated for plot)",
        title_suffix=" (consensus)",
        mode_label="consensus",
        bedroc_alpha=args.bedroc_alpha,
        logauc_lambda=args.logauc_lambda,
        class_min_targets=args.class_min_targets,
        write_placeholder_images=args.class_placeholder_images,
    )
    _, reranked_class_rows = _emit_class_aggregate_plots(
        analysis_root=analysis_root,
        run_id=active_run_id,
        target_specs_by_key=target_specs_by_key,
        eval_results=reranked_eval_results,
        included_target_keys=reranked_target_keys,
        mode_slug="post_docked_scorch",
        score_high_is_better=True,
        hist_xlabel="Reranked SCORCH score (higher = better)",
        title_suffix=" (reranked_scorch)",
        mode_label="post_docked_scorch",
        bedroc_alpha=args.bedroc_alpha,
        logauc_lambda=args.logauc_lambda,
        class_min_targets=args.class_min_targets,
        write_placeholder_images=args.class_placeholder_images,
    )
    class_generation_rows.extend(consensus_class_rows)
    class_generation_rows.extend(reranked_class_rows)
    class_report_columns = [
        "run_id",
        "mode",
        "class_name",
        "class_slug",
        "status_reason",
        "aggregation_method",
        "targets_total_in_class",
        "targets_present_in_run",
        "targets_evaluable",
        "pooled_rows",
        "pooled_actives",
        *_class_metric_base_columns(args.bedroc_alpha),
        *_class_metric_stat_columns(args.bedroc_alpha),
    ]
    class_report_df = pd.DataFrame(class_generation_rows, columns=class_report_columns)
    class_report_path = analysis_root / "class_aggregates" / "class_generation_report.tsv"
    class_report_path.parent.mkdir(parents=True, exist_ok=True)
    class_report_df.to_csv(class_report_path, sep="\t", index=False)
    if getattr(args, "pretty_summary", True):
        _write_pretty_table_noformat(
            class_report_df,
            class_report_path.with_name(class_report_path.stem + "_pretty.txt"),
            title="Class aggregate generation report",
        )

    control_long_records: List[Dict] = []
    control_summary_records: List[Dict] = []
    scanned_ok: List[str] = []
    missing_logs: List[str] = []
    parsed_counts: Dict[str, Dict[str, int]] = {}
    if getattr(args, "emit_control_report", False):
        global _CONTROL_PATTERNS_LOGGED
        if not _CONTROL_PATTERNS_LOGGED:
            dbg(
                "DEBUG",
                "control",
                f'searching_pattern.control_centers="{CONTROL_CENTERS_RE.pattern}"',
            )
            dbg(
                "DEBUG",
                "control",
                f'searching_pattern.control_redock="{CONTROL_REDOCK_RE.pattern}"',
            )
            _CONTROL_PATTERNS_LOGGED = True
        for spec in control_targets:
            eval_meta = target_eval_results.get(spec.target_key)
            run_label = _resolve_run_label(
                spec.pdb_id, eval_meta, args.run_id, active_run_id
            )
            target_label = names.get(spec.pdb_id, "")
            library_label = ""
            if args.report_library and manifest_data:
                library_label = _manifest_library_for_target(
                    manifest=manifest_data,
                    pdb_id=spec.pdb_id,
                    variant_label=spec.variant,
                    ph_tag=spec.ph_tag,
                )
            if not library_label and not df.empty and "library_name" in df.columns:
                mask = df["pdb_id"].astype(str) == str(spec.pdb_id)
                if "variant" in df.columns:
                    mask &= df["variant"].astype(str) == (spec.variant or "")
                if "pH" in df.columns:
                    mask &= df["pH"].astype(str) == (spec.ph_tag or "")
                if mask.any():
                    library_label = str(df.loc[mask, "library_name"].iloc[0])

            csv_path = csv_paths_by_target.get(spec.target_key)
            if csv_path is None:
                continue
            parsed_counts[spec.target_key] = {
                "centers": 0,
                "redock_lines": 0,
                "report_rows": 0,
            }
            selected_log, candidates = _candidate_protein_logs(
                spec.pdb_id,
                csv_path,
                docked_root,
                log_root_override,
                cfg,
            )
            dbg(
                "DEBUG",
                "control",
                f"target={spec.target_key} log_candidates={[str(c) for c in candidates]}",
            )
            if selected_log:
                dbg(
                    "INFO",
                    "control",
                    f"target={spec.target_key} protein_log={selected_log}",
                )
                scanned_ok.append(spec.target_key)
            else:
                missing_logs.append(spec.target_key)
                missing_hint = (
                    str(candidates[0])
                    if candidates
                    else str(Path("docked") / spec.pdb_id / "protein.log")
                )
                dbg(
                    "WARN",
                    "control",
                    f"target={spec.target_key} protein_log_missing={missing_hint} action=skip_control_parse",
                )

            long_rows, summary_row, meta = _build_control_records(
                spec.pdb_id,
                run_label,
                target_label,
                library_label,
                selected_log,
                spec.variant,
                spec.ph_tag,
            )
            if spec.target_key not in parsed_counts:
                parsed_counts[spec.target_key] = {
                    "centers": 0,
                    "redock_lines": 0,
                    "report_rows": 0,
                }
            parsed_counts[spec.target_key]["centers"] = meta.get("centers_found", 0)
            parsed_counts[spec.target_key]["redock_lines"] = meta.get("redock_lines", 0)
            parsed_counts[spec.target_key]["report_rows"] = meta.get("report_rows", 0)
            dbg(
                "INFO",
                "control",
                (
                    f"target={spec.target_key} centers_found={parsed_counts[spec.target_key]['centers']} "
                    f"redock_lines={parsed_counts[spec.target_key]['redock_lines']} "
                    f"report_rows={parsed_counts[spec.target_key]['report_rows']}"
                ),
            )
            if long_rows:
                control_long_records.extend(long_rows)
            control_summary_records.append(summary_row)
        for spec in control_targets:
            if spec.target_key not in parsed_counts:
                parsed_counts[spec.target_key] = {
                    "centers": 0,
                    "redock_lines": 0,
                    "report_rows": 0,
                }
                if (
                    spec.target_key not in missing_logs
                    and spec.target_key not in scanned_ok
                ):
                    missing_logs.append(spec.target_key)

    summary_name = "summary.tsv"
    if args.run_id:
        run_suffix = str(args.run_id).replace(" ", "")
        summary_name = f"summary{run_suffix}.tsv"
    summary_path = analysis_root / summary_name

    variant_col_present = "variant" in df.columns
    variant_upper = (
        df["variant"].astype(str).str.upper() if variant_col_present else None
    )
    apo_mask = (variant_upper == "APO") if variant_col_present else None
    holo_mask = (variant_upper == "HOLO") if variant_col_present else None
    has_sections = bool(
        variant_col_present
        and (
            (apo_mask is not None and apo_mask.any())
            or (holo_mask is not None and holo_mask.any())
        )
    )

    preferred_cols = (
        "variant",
        "pH",
        "run_id",
        "target_name",
        "library_name",
        "pdb_id",
    )
    excluded_cols = {
        "run_id",
        "target_name",
        "library_name",
        "pdb_id",
        "N",
        "n_actives",
        "actives_fraction",
        "variant",
        "pH",
        "status_reason",
    }

    with open(summary_path, "w", newline="") as fh:
        if has_sections:
            _df = df.copy()
            cols = [c for c in preferred_cols if c in _df.columns] + [
                c for c in _df.columns if c not in preferred_cols
            ]
            apo_rows = (
                df.loc[apo_mask].sort_values("pdb_id")
                if apo_mask is not None
                else pd.DataFrame()
            )
            holo_rows = (
                df.loc[holo_mask].sort_values("pdb_id")
                if holo_mask is not None
                else pd.DataFrame()
            )
            n_apo = len(apo_rows)
            n_holo = len(holo_rows)

            sections = []
            leftover_rows = pd.DataFrame()
            if n_apo:
                apo_out = apo_rows
                a_cols = [c for c in preferred_cols if c in apo_out.columns] + [
                    c for c in apo_out.columns if c not in preferred_cols
                ]
                fh.write("Apo\n")
                apo_out[a_cols].to_csv(fh, sep="\t", index=False, float_format="%.3f")
                sections.append(("Apo", apo_out[a_cols]))

            if n_holo:
                if n_apo:
                    fh.write("\n")
                holo_out = holo_rows
                h_cols = [c for c in preferred_cols if c in holo_out.columns] + [
                    c for c in holo_out.columns if c not in preferred_cols
                ]
                fh.write("Holo\n")
                holo_out[h_cols].to_csv(fh, sep="\t", index=False, float_format="%.3f")
                sections.append(("Holo", holo_out[h_cols]))

            leftover_mask = (
                ~(apo_mask | holo_mask)
                if (apo_mask is not None and holo_mask is not None)
                else pd.Series(False, index=df.index)
            )
            leftover_rows = (
                df.loc[leftover_mask].sort_values("pdb_id")
                if not df.empty
                else pd.DataFrame()
            )
            if not leftover_rows.empty:
                if n_apo or n_holo:
                    fh.write("\n")
                lo_base = leftover_rows
                l_cols = [c for c in preferred_cols if c in lo_base.columns] + [
                    c for c in lo_base.columns if c not in preferred_cols
                ]
                lo_base[l_cols].to_csv(fh, sep="\t", index=False)
                sections.append(("Unlabeled", lo_base[l_cols]))
            dbg(
                "INFO",
                "summary",
                f"sections=Apo:{n_apo} Holo:{n_holo} other={len(leftover_rows)}",
            )

            # Pretty summary file (same run-id naming as TSV)
            if args.pretty_summary and sections:
                pretty_name = summary_path.with_name(summary_path.stem + "_pretty.txt")
                _write_pretty_summary(sections, pretty_name)
                print(f"[eval] pretty_summary_out={pretty_name}")
        else:
            _df = df.copy()
            cols = [c for c in preferred_cols if c in _df.columns] + [
                c for c in _df.columns if c not in preferred_cols
            ]
            _df[cols].to_csv(fh, sep="\t", index=False)
            # Pretty summary alongside TSV
            if args.pretty_summary:  # if you added the flag; otherwise remove the 'if'
                pretty_name = summary_path.with_name(summary_path.stem + "_pretty.txt")
                _write_pretty_summary([(None, _df[cols])], pretty_name)
                print(f"[eval] pretty_summary_out={pretty_name}")

    metric_cols = [c for c in df.columns if c not in excluded_cols]
    class_macro_metric_cols = metric_cols + [
        c for c in _class_metric_stat_columns(args.bedroc_alpha) if c not in metric_cols
    ]
    macro_rows: List[Dict[str, Any]] = []
    macro = df[metric_cols].mean(numeric_only=True).to_dict()
    macro_status = "ok" if not df.empty else "empty_summary"
    macro_rows.append(
        {
            "pdb_id": "macro_avg",
            **{k: macro.get(k) for k in class_macro_metric_cols},
            "status_reason": macro_status,
        }
    )

    docking_class_status: Dict[str, str] = {}
    for row in class_generation_rows:
        if str(row.get("mode", "")) != "docking":
            continue
        class_name = str(row.get("class_name", "")).strip()
        if class_name:
            docking_class_status[class_name] = str(
                row.get("status_reason", "no_targets_in_run")
            )

    for class_name in DUD_E_PROTEIN_CLASSES:
        class_slug = _class_slug(class_name)
        class_eval = docking_class_results.get(class_name)
        metric_payload: Dict[str, Any] = {}
        status_reason = docking_class_status.get(class_name, "no_targets_in_run")
        if class_eval is not None and class_eval.metrics is not None:
            source = (
                class_eval.metrics.to_dict()
                if hasattr(class_eval.metrics, "to_dict")
                else dict(class_eval.metrics)
            )
            metric_payload = {k: source.get(k) for k in class_macro_metric_cols}
            status_reason = str(source.get("status_reason") or class_eval.status_reason)
        else:
            metric_payload = {k: float("nan") for k in class_macro_metric_cols}
            if class_eval is not None and class_eval.status_reason:
                status_reason = class_eval.status_reason

        macro_rows.append(
            {
                "pdb_id": f"class_{class_slug}",
                **metric_payload,
                "status_reason": status_reason,
            }
        )

    macro_df = pd.DataFrame(macro_rows)
    macro_columns = ["pdb_id", *class_macro_metric_cols, "status_reason"]
    for col in macro_columns:
        if col not in macro_df.columns:
            macro_df[col] = "" if col in {"pdb_id", "status_reason"} else float("nan")
    macro_df = macro_df.loc[:, macro_columns]
    macro_path = analysis_root / "summary_macro.tsv"
    macro_df.to_csv(macro_path, sep="\t", index=False)
    emit_consensus_summary(
        df_consensus,
        analysis_root,
        run_label_for_files,
        getattr(args, "pretty_summary", True),
    )
    emit_reranked_scorch_summary(
        df_reranked_all,
        analysis_root,
        run_label_for_files,
        getattr(args, "pretty_summary", True),
    )

    if getattr(args, "emit_control_report", False):
        long_columns = [
            "variant",
            "pH",
            "pdb_id",
            "run_id",
            "target_name",
            "library_name",
            "controls_found",
            "control_id",
            "rmsd_to_crystal_A",
            "redock_best_energy_kcal_mol",
            "chosen",
            "max_spread_A",
            "policy",
        ]
        report_df = pd.DataFrame(control_long_records, columns=long_columns)
        control_report_path = analysis_root / "control_redock_report.tsv"
        report_df.to_csv(control_report_path, sep="\t", index=False)

        control_summary_df = pd.DataFrame(control_summary_records, columns=long_columns)
        for col in ("variant", "pH"):
            if col in control_summary_df.columns:
                s = control_summary_df[col]
                if s.isna().all() or (s.astype(str).str.strip() == "").all():
                    control_summary_df.drop(columns=[col], inplace=True)
        control_summary_path = analysis_root / "control_redock_summary.tsv"
        control_summary_df.to_csv(control_summary_path, sep="\t", index=False)

        dbg(
            "INFO", "control", f"report_rows={len(report_df)} out={control_report_path}"
        )
        dbg(
            "INFO",
            "control",
            f"summary_rows={len(control_summary_df)} out={control_summary_path}",
        )
        # Pretty control summary (no schema/value changes)
        if getattr(args, "pretty_summary", True):
            control_pretty_path = control_summary_path.with_name(
                control_summary_path.stem + "_pretty.txt"
            )
            _write_pretty_table_noformat(
                control_summary_df, control_pretty_path, title="Control Redock Summary"
            )
            print(f"[control] pretty_summary_out={control_pretty_path}")

    total_targets = len(control_targets)
    total_scanned = len(scanned_ok)
    total_missing = len(missing_logs)
    scanned_label = ",".join(scanned_ok) if scanned_ok else "(none)"
    missing_label = ",".join(missing_logs) if missing_logs else "(none)"
    total_redock = sum(stats.get("redock_lines", 0) for stats in parsed_counts.values())
    total_centers = sum(stats.get("centers", 0) for stats in parsed_counts.values())
    total_reports = sum(stats.get("report_rows", 0) for stats in parsed_counts.values())

    dbg("INFO", "summary", f"n_targets_total={total_targets}")
    dbg(
        "INFO",
        "summary",
        f"n_logs_scanned={total_scanned} n_logs_missing={total_missing}",
    )
    dbg("INFO", "summary", f"scanned_ok={scanned_label}")
    dbg("INFO", "summary", f"missing_logs={missing_label}")
    dbg(
        "INFO",
        "summary",
        (
            "totals: "
            f"redock_lines={total_redock} "
            f"control_centers_found={total_centers} "
            f"report_rows={total_reports}"
        ),
    )

    scan_summary_path = analysis_root / "control_redock_scan_summary.txt"
    per_target_rows: List[Dict[str, int | str]] = []
    for spec in control_targets:
        stats = parsed_counts.get(
            spec.target_key, {"centers": 0, "redock_lines": 0, "report_rows": 0}
        )
        per_target_rows.append(
            {
                "target_key": spec.target_key,
                "pdb_id": spec.pdb_id,
                "variant": spec.variant or "",
                "pH": spec.ph_tag or "",
                "centers_found": stats.get("centers", 0),
                "redock_lines": stats.get("redock_lines", 0),
                "report_rows": stats.get("report_rows", 0),
            }
        )
    with open(scan_summary_path, "w", encoding="utf-8") as fh:
        fh.write("Control Log Scan Summary\n")
        fh.write(f"n_targets_total: {total_targets}\n")
        fh.write(f"n_logs_scanned: {total_scanned}\n")
        fh.write(f"n_logs_missing: {total_missing}\n")
        fh.write(f"scanned_ok: {scanned_label}\n")
        fh.write(f"missing_logs: {missing_label}\n")
        fh.write(
            (
                "totals: "
                f"redock_lines={total_redock} "
                f"control_centers_found={total_centers} "
                f"report_rows={total_reports}\n"
            )
        )
        fh.write("\nPer-target counts\n")
        if per_target_rows:
            counts_df = pd.DataFrame(
                per_target_rows,
                columns=[
                    "target_key",
                    "pdb_id",
                    "variant",
                    "pH",
                    "centers_found",
                    "redock_lines",
                    "report_rows",
                ],
            )
            fh.write(counts_df.to_string(index=False))
            fh.write("\n")
        else:
            fh.write("(none)\n")

    dbg("INFO", "summary", f"scan_summary_pretty={scan_summary_path}")

    print(f"[dbg.summary] header={list(_df[cols].columns)}")
    dbg("INFO", "summary", f"out={summary_path} columns={metric_cols}")
    dbg("DEBUG", "summary", f"targets_written={len(df)} macro_path={macro_path}")

# ruff: noqa: F403, F405
from analysis.reporting.run_report.constants import *
from analysis.reporting.fda_name_map import (
    mapping_file_provenance_path,
    mapping_file_sha256,
)
from analysis.reporting.run_report.config_resolve import _resolve_target_name_cfg
from analysis.reporting.run_report.display import (
    _matches_query,
    _parse_highlight_queries,
    _resolve_ligand_display,
    _summarize_display_resolution,
)
from analysis.reporting.run_report.export_io import _format_num
from analysis.reporting.run_report.master_csv import _resolve_master_csvs, _safe_relpath
from analysis.reporting.run_report.target_stats import (
    _compute_fdr_stats,
    _format_pct_display,
    _ligand_identity_for_library,
    _load_scorch_stats,
    _pick_repeated_value,
    _resolve_manifest_pocket,
    _resolve_target_name,
    _row_z_selected_rank_sort_key,
    _target_display_name,
    _total_library_n,
    _vector_from_row,
)
def build_report(
    run_id: str,
    repo_root: Path,
    top_n: int = 5,
    extended_top_n: int = 15,
    notable_pct: float = 0.01,
    notable_max: int = 25,
    decoy_prefix: str = DECOY_PREFIX_DEFAULT,
    fda_mapping_csv: Optional[str] = None,
    highlight_ligands: Optional[str] = "imatinib",
    highlight_match: str = "any_contains",
    highlight_max_per_target: int = 5,
    multi_target_min_targets: int = 2,
    multi_target_max_pct: float = 0.05,
    multi_target_max_hits: int = 50,
    multi_target_sort: str = "worst_pct_then_mean",
    filter_invalid: bool = False,
    top_targets_n: int = 3,
    breadth_pct_threshold: float = 0.01,
    coverage_requires_pose_valid: bool = False,
    pct_display_decimals: int = 1,
    tail_median_decimals: int = 0,
    highlights_top_pct: float = 0.01,
    highlights_max_ligands: int = 25,
    test_mode_tokens: Optional[List[str]] = None,
    master_csvs_override: Optional[List[Path]] = None,
) -> Dict[str, Any]:
    # Avoid cross-run leakage when multiple repos/tests reuse the same process.
    _TARGET_NAME_CACHE.clear()

    tokens = test_mode_tokens or resolve_test_mode_tokens(repo_root, run_id)
    if master_csvs_override:
        master_csvs = [Path(p) for p in master_csvs_override]
    else:
        master_csvs = _resolve_master_csvs(repo_root, run_id, tokens, decoy_prefix)
    if not master_csvs:
        raise FileNotFoundError(
            "Report input check failed: master_rows CSV is missing. "
            f"run_id={run_id} tokens={tokens} checked_dir={repo_root / 'data' / run_id}. "
            "Generate docking/report inputs first, then re-run report generation. "
            "Use `atlas --doctor` to confirm run paths."
        )
    primary_master_csv = master_csvs[0]

    canonical_manifest_rel = Path("manifests") / run_id / "run_manifest.yaml"
    manifest_data, _manifest_path = load_run_manifest(repo_root, run_id)

    top_n = max(0, int(top_n))
    extended_top_n = max(0, int(extended_top_n))
    notable_max = max(0, int(notable_max))
    highlight_max_per_target = max(0, int(highlight_max_per_target))
    multi_target_min_targets = max(0, int(multi_target_min_targets))
    multi_target_max_hits = max(0, int(multi_target_max_hits))
    multi_target_max_pct = (
        float(multi_target_max_pct) if multi_target_max_pct is not None else 0.0
    )
    top_targets_n = max(1, int(top_targets_n))
    breadth_pct_threshold = (
        float(breadth_pct_threshold) if breadth_pct_threshold is not None else 0.01
    )
    if not math.isfinite(breadth_pct_threshold) or breadth_pct_threshold < 0:
        breadth_pct_threshold = 0.01
    pct_display_decimals = max(0, int(pct_display_decimals))
    tail_median_decimals = max(0, int(tail_median_decimals))
    highlights_max_ligands = max(0, int(highlights_max_ligands))
    highlights_top_pct = (
        float(highlights_top_pct) if highlights_top_pct is not None else 0.01
    )
    if not math.isfinite(highlights_top_pct) or highlights_top_pct < 0:
        highlights_top_pct = 0.01

    logger = logging.getLogger("run-report")
    rows: List[Dict[str, Any]] = []
    has_ligand_display = False
    for master_csv in master_csvs:
        with master_csv.open("r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            if "ligand_display" in (reader.fieldnames or []):
                has_ligand_display = True
            for r in reader:
                rows.append(r)

    highlight_queries = _parse_highlight_queries(highlight_ligands)
    if highlight_max_per_target:
        highlight_queries = highlight_queries[:highlight_max_per_target]
    else:
        highlight_queries = []

    fda_index = None
    mapping_csv = resolve_mapping_csv_path(repo_root, run_id, cli_value=fda_mapping_csv)
    if mapping_csv is None:
        logger.warning(
            "%s action=fda_mapping status=missing run_id=%s", COMPONENT, run_id
        )
    else:
        fda_index = try_load_fda_index(mapping_csv)
        if fda_index is None:
            logger.warning(
                "%s action=fda_mapping status=load_failed path=%s",
                COMPONENT,
                mapping_csv,
            )
        else:
            logger.info(
                "%s action=fda_mapping status=loaded path=%s",
                COMPONENT,
                mapping_csv,
            )

    _summarize_display_resolution(
        rows,
        has_ligand_display,
        fda_index,
        mapping_csv,
        logger,
    )
    target_name_cfg = _resolve_target_name_cfg(repo_root, run_id)

    # Group by target combo
    # Key: "<pdb_id>|<variant>|<ph_label>"
    target_groups: Dict[str, List[Dict[str, Any]]] = {}
    ligand_top5_hits: Dict[str, List[Dict[str, Any]]] = {}
    ligand_notable_hits: Dict[str, List[Dict[str, Any]]] = {}
    ligand_display_map: Dict[str, str] = {}
    multi_target_entries: Dict[str, List[Dict[str, Any]]] = {}
    compact_ligand_entries: Dict[str, List[Dict[str, Any]]] = {}

    unique_ligands = set()

    for r in rows:
        pdb = r.get("pdb_id", "").strip()
        variant = r.get("variant", "").strip()
        ph = r.get("ph_label", "").strip()
        key = build_target_id(pdb, variant, ph)
        target_groups.setdefault(key, []).append(r)

        lig_base = clean_report_text(r.get("ligand_base"))
        if lig_base:
            unique_ligands.add(lig_base)

    summary = {
        "n_target_combos": len(target_groups),
        "n_rows": len(rows),
        "n_unique_ligands": len(unique_ligands),
    }

    targets_out = {}

    sorted_keys = sorted(target_groups.keys())

    for key in sorted_keys:
        group_rows = target_groups[key]
        tid = parse_target_id(key)
        pdb, variant, ph = tid.pdb_id, tid.variant, tid.ph_label
        target_name = _resolve_target_name(pdb, target_name_cfg)
        target_display = _target_display_name(pdb, target_name)

        # QC Stats
        n_rows = len(group_rows)
        n_decoys = sum(1 for r in group_rows if as_report_row_bool(r.get("is_decoy")))
        n_controls = sum(1 for r in group_rows if as_report_row_bool(r.get("is_control")))

        # Extract metrics from first row (assuming they are repeated/consistent per target)
        first = group_rows[0]
        ef1 = as_float(first.get("ef1"))
        roc_auc = as_float(first.get("roc_auc"))
        roc_auc_adj = as_float(first.get("roc_auc_adj")) or roc_auc
        if ef1 is None:
            ef1 = 0.0
        if roc_auc is None:
            roc_auc = 0.0
        if roc_auc_adj is None:
            roc_auc_adj = roc_auc
        dud_reason = first.get("dud_eval_status_reason") or None

        # Pocket info
        # Prefer row data if present (it was joined in master export)
        # master export columns: pocket_method, center_x, etc.
        p_method = (first.get("pocket_method") or "").strip() or None
        center = _vector_from_row(first, "center")
        box = _vector_from_row(first, "box")

        manifest_pocket = _resolve_manifest_pocket(
            manifest_data, pdb.upper(), variant.upper(), ph
        )
        if not p_method:
            p_method = manifest_pocket.get("method")
        if center is None:
            center = manifest_pocket.get("center")
        if box is None:
            box = manifest_pocket.get("box")

        # Rank non-decoys with valid z_selected
        rankable: List[tuple[Dict[str, Any], float]] = []
        for r in group_rows:
            if as_report_row_bool(r.get("is_decoy")):
                continue
            if filter_invalid and not as_report_row_bool(r.get("pose_valid_any")):
                continue
            t_val = as_float(r.get("z_selected"))
            if t_val is None or not math.isfinite(t_val):
                continue
            rankable.append((r, t_val))

        sorted_rankable = sorted(rankable, key=_row_z_selected_rank_sort_key)
        ranked_entries = []
        ranked_rows: List[Dict[str, Any]] = []
        total_library_n = _total_library_n(group_rows)
        for idx, (r, t_sel) in enumerate(sorted_rankable, start=1):
            lig_base = clean_report_text(r.get("ligand_base"))
            lig_display = _resolve_ligand_display(r, has_ligand_display, fda_index)
            t_src = r.get("z_selected_source") or None
            pct = idx / total_library_n

            entry = {
                "ligand_base": lig_base,
                "ligand_display": lig_display,
                "z_selected": t_sel,
                "z_selected_source": t_src,
                "pose_valid_any": as_report_row_bool(r.get("pose_valid_any")),
                "pose_invalid_reason_top": clean_report_text(
                    r.get("pose_invalid_reason_top")
                ),
                "rank": idx,
                "pct_rank": float(f"{pct:.6f}"),
                "is_control": as_report_row_bool(r.get("is_control")),
                "is_decoy": as_report_row_bool(r.get("is_decoy")),
            }
            ranked_entries.append(entry)
            ranked_rows.append({"row": r, "entry": entry})
            if lig_base and lig_base not in ligand_display_map:
                ligand_display_map[lig_base] = lig_display or lig_base

        compact_best_by_ligand: Dict[str, Dict[str, Any]] = {}
        for r in group_rows:
            if as_report_row_bool(r.get("is_decoy")):
                continue
            t_val = as_float(r.get("z_selected"))
            if t_val is None or not math.isfinite(t_val):
                continue
            if coverage_requires_pose_valid and not as_report_row_bool(r.get("pose_valid_any")):
                continue
            lig_base = _ligand_identity_for_library(r)
            if not lig_base:
                continue
            lig_display = _resolve_ligand_display(r, has_ligand_display, fda_index)
            lig_file = clean_report_text(r.get("ligand_file") or r.get("ligand") or "")
            candidate = {
                "ligand_base": lig_base,
                "ligand_display": lig_display or lig_base,
                "ligand_file": lig_file,
                "z_selected": t_val,
                "z_selected_source": r.get("z_selected_source") or None,
            }
            current = compact_best_by_ligand.get(lig_base)
            if current is None:
                compact_best_by_ligand[lig_base] = candidate
                continue
            current_t = as_float(current.get("z_selected"))
            if current_t is None or t_val > current_t:
                compact_best_by_ligand[lig_base] = candidate
                continue
            if (
                current_t is not None
                and t_val == current_t
                and lig_file < str(current.get("ligand_file") or "")
            ):
                compact_best_by_ligand[lig_base] = candidate

        def compact_rank_key(item: Dict[str, Any]) -> tuple:
            t_val = as_float(item.get("z_selected"))
            return (
                -(t_val if t_val is not None else float("-inf")),
                clean_report_text(item.get("ligand_base")),
                clean_report_text(item.get("ligand_file")),
            )

        compact_rankable = sorted(compact_best_by_ligand.values(), key=compact_rank_key)
        for idx, item in enumerate(compact_rankable, start=1):
            lig_base = clean_report_text(item.get("ligand_base"))
            if not lig_base:
                continue
            lig_display = clean_report_text(item.get("ligand_display")) or lig_base
            if lig_base not in ligand_display_map:
                ligand_display_map[lig_base] = lig_display
            pct = idx / total_library_n
            compact_ligand_entries.setdefault(lig_base, []).append(
                {
                    "target_id": key,
                    "pdb_id": pdb,
                    "target_name": target_name,
                    "target_display": target_display,
                    "rank": idx,
                    "pct_rank": float(f"{pct:.6f}"),
                    "z_selected": float(as_float(item.get("z_selected")) or 0.0),
                    "z_selected_source": item.get("z_selected_source"),
                    "ligand_display": lig_display,
                }
            )

        top_list = ranked_entries[:top_n] if top_n > 0 else []
        extended_list = ranked_entries[:extended_top_n] if extended_top_n > 0 else []

        notable_candidates = []
        for entry in ranked_entries:
            pct_rank = entry.get("pct_rank")
            is_control = bool(entry.get("is_control"))
            if isinstance(pct_rank, (int, float)):
                if pct_rank <= notable_pct or is_control:
                    notable_candidates.append(entry)
        if notable_max <= 0:
            notable_list = []
        elif len(notable_candidates) <= notable_max:
            notable_list = notable_candidates
        else:
            controls = [e for e in notable_candidates if e["is_control"]]
            non_controls = [e for e in notable_candidates if not e["is_control"]]
            controls_sorted = sorted(controls, key=lambda e: int(e.get("rank") or 0))
            non_controls_sorted = sorted(
                non_controls, key=lambda e: int(e.get("rank") or 0)
            )
            notable_list = controls_sorted[:notable_max]
            remaining = notable_max - len(notable_list)
            if remaining > 0:
                notable_list.extend(non_controls_sorted[:remaining])

        highlights = []
        if highlight_queries:
            for query in highlight_queries:
                candidates = []
                for item in ranked_rows:
                    row = item["row"]
                    entry = item["entry"]
                    ligand_raw = clean_report_text(
                        row.get("ligand") or row.get("ligand_file") or ""
                    )
                    if _matches_query(
                        query,
                        entry["ligand_display"],
                        entry["ligand_base"],
                        ligand_raw,
                        highlight_match,
                        fda_index,
                    ):
                        candidates.append(item)

                if not candidates:
                    highlights.append(
                        {
                            "query": query,
                            "found": False,
                            "pose_valid_any": False,
                            "pose_invalid_reason_top": "",
                        }
                    )
                    continue

                def _highlight_sort_key(item: Dict[str, Any]) -> tuple:
                    entry = item["entry"]
                    t_val = entry["z_selected"]
                    t_sort = -(t_val if t_val is not None else float("-inf"))
                    lig_base = entry["ligand_base"]
                    return (entry["rank"], t_sort, lig_base)

                best = min(candidates, key=_highlight_sort_key)
                row = best["row"]
                entry = best["entry"]
                highlight_entry = {
                    "query": query,
                    "found": True,
                    "ligand_display": entry["ligand_display"],
                    "ligand_base": entry["ligand_base"],
                    "z_selected": entry["z_selected"],
                    "z_selected_source": entry["z_selected_source"],
                    "pose_valid_any": entry["pose_valid_any"],
                    "pose_invalid_reason_top": entry["pose_invalid_reason_top"],
                    "rank": entry["rank"],
                    "pct_rank": entry["pct_rank"],
                    "is_control": entry["is_control"],
                    "is_decoy": entry["is_decoy"],
                }
                library = clean_report_text(row.get("library"))
                if library:
                    highlight_entry["library"] = library
                highlights.append(highlight_entry)

        for item in ranked_rows:
            row = item["row"]
            entry = item["entry"]
            lig_base = entry["ligand_base"]
            if not lig_base:
                continue
            detail = {
                "target_id": key,
                "rank": entry["rank"],
                "pct_rank": entry["pct_rank"],
                "z_selected": entry["z_selected"],
                "z_selected_source": entry["z_selected_source"],
                "ligand_display": entry["ligand_display"],
                "is_control": entry["is_control"],
            }
            library = clean_report_text(row.get("library"))
            if library:
                detail["library"] = library
            multi_target_entries.setdefault(lig_base, []).append(detail)

        for entry in top_list:
            lig_base = clean_report_text(entry.get("ligand_base"))
            if lig_base:
                ligand_top5_hits.setdefault(lig_base, []).append(
                    {
                        "target": key,
                        "z_selected": entry["z_selected"],
                        "z_selected_source": entry["z_selected_source"],
                    }
                )

        for entry in notable_list:
            lig_base = clean_report_text(entry.get("ligand_base"))
            if lig_base:
                ligand_notable_hits.setdefault(lig_base, []).append(
                    {
                        "target": key,
                        "z_selected": entry["z_selected"],
                        "z_selected_source": entry["z_selected_source"],
                        "rank": entry["rank"],
                        "pct_rank": entry["pct_rank"],
                    }
                )

        fdr_stats = _compute_fdr_stats(group_rows)
        fdr_value: Optional[Any] = None
        if fdr_stats:
            best_q = as_float(fdr_stats.get("best_q"))
            if best_q is None:
                best_q = 1.0
            reliable = bool(fdr_stats.get("reliable"))
            fdr_value = float(best_q) if reliable else f"SMALL {best_q:.6g}"

        targets_out[key] = {
            "target_name": target_name,
            "qc": {
                "ef1": ef1,
                "roc_auc": roc_auc,
                "roc_auc_adj": roc_auc_adj,
                "dud_eval_status_reason": dud_reason,
                "n_rows": n_rows,
                "n_decoys": n_decoys,
                "n_controls": n_controls,
                "fdr": fdr_value,
                "decoy_stats": {
                    "consensus": {
                        "mu": _pick_repeated_value(group_rows, "consensus_mu_decoy"),
                        "sigma": _pick_repeated_value(
                            group_rows, "consensus_sigma_decoy"
                        ),
                        "n": _pick_repeated_value(
                            group_rows, "consensus_n_decoys", is_int=True
                        ),
                    },
                    "blend": {
                        "mu": _pick_repeated_value(group_rows, "blend_mu_decoy"),
                        "sigma": _pick_repeated_value(group_rows, "blend_sigma_decoy"),
                        "n": _pick_repeated_value(
                            group_rows, "blend_n_decoys", is_int=True
                        ),
                    },
                    "scorch": _load_scorch_stats(
                        repo_root,
                        run_id,
                        pdb,
                        variant,
                        ph,
                        decoy_prefix=decoy_prefix,
                    ),
                },
                "fdr_stats": fdr_stats if fdr_stats else None,
            },
            "pocket": {"method": p_method, "center": center, "box": box},
            "top5_ligands": top_list,
            "top_ligands_extended": extended_list,
            "notable_ligands": notable_list,
            "highlights": highlights,
        }

    # Ligands section
    ligands_out = {}
    all_bases = sorted(set(ligand_top5_hits.keys()) | set(ligand_notable_hits.keys()))
    for base in all_bases:
        entry = {"ligand_display": ligand_display_map.get(base, base)}

        if base in ligand_top5_hits:
            hits = ligand_top5_hits[base]
            hits.sort(key=lambda x: x["target"])
            entry["targets_in_top5"] = hits

        if base in ligand_notable_hits:
            hits = ligand_notable_hits[base]
            hits.sort(key=lambda x: x["target"])
            entry["targets_notable"] = hits

        ligands_out[base] = entry

    multi_target_hits = []
    if multi_target_max_hits > 0 and multi_target_entries:
        for lig_base, entries in sorted(multi_target_entries.items()):
            qualifying = [
                e
                for e in entries
                if isinstance(e.get("pct_rank"), (int, float))
                and e["pct_rank"] <= multi_target_max_pct
            ]
            if multi_target_min_targets and len(qualifying) < multi_target_min_targets:
                continue
            if not qualifying:
                continue

            best_entry = min(
                qualifying,
                key=lambda e: (e["rank"], e.get("ligand_display") or ""),
            )
            display = best_entry.get("ligand_display") or lig_base

            pct_values = [float(e["pct_rank"]) for e in qualifying]
            rank_values = [int(e["rank"]) for e in qualifying]
            t_values = [float(e["z_selected"]) for e in qualifying]
            worst_pct = max(pct_values)
            mean_pct = sum(pct_values) / len(pct_values)
            best_rank = min(rank_values)
            best_z_selected = max(t_values)
            rank_sum = sum(rank_values)

            targets_sorted = sorted(
                qualifying, key=lambda e: (e["rank"], e["target_id"])
            )
            targets_out_list = []
            for entry in targets_sorted:
                target_entry = {
                    "target_id": entry["target_id"],
                    "rank": entry["rank"],
                    "pct_rank": entry["pct_rank"],
                    "z_selected": entry["z_selected"],
                    "z_selected_source": entry["z_selected_source"],
                }
                library = clean_report_text(entry.get("library"))
                if library:
                    target_entry["library"] = library
                targets_out_list.append(target_entry)

            multi_target_hits.append(
                {
                    "ligand_base": lig_base,
                    "ligand_display": display,
                    "targets_qualified": len(qualifying),
                    "worst_pct": float(f"{worst_pct:.6f}"),
                    "mean_pct": float(f"{mean_pct:.6f}"),
                    "best_rank": best_rank,
                    "best_z_selected": best_z_selected,
                    "targets": targets_out_list,
                    "_rank_sum": rank_sum,
                }
            )

        def _multi_sort_key(entry: Dict[str, Any]) -> tuple:
            if multi_target_sort == "mean_pct":
                return (
                    entry["mean_pct"],
                    entry["worst_pct"],
                    entry["best_rank"],
                    entry["ligand_base"],
                )
            if multi_target_sort == "best_rank_sum":
                return (
                    entry["_rank_sum"],
                    entry["worst_pct"],
                    entry["mean_pct"],
                    entry["ligand_base"],
                )
            return (
                entry["worst_pct"],
                entry["mean_pct"],
                entry["best_rank"],
                entry["ligand_base"],
            )

        multi_target_hits = sorted(multi_target_hits, key=_multi_sort_key)[
            :multi_target_max_hits
        ]
        for entry in multi_target_hits:
            entry.pop("_rank_sum", None)

    ligand_highlights = []
    total_targets = len(target_groups)
    for lig_base, entries in sorted(compact_ligand_entries.items()):
        if not entries:
            continue
        entries_sorted = sorted(
            entries,
            key=lambda item: (
                float(as_float(item.get("pct_rank")) or 1.0),
                -(float(as_float(item.get("z_selected")) or float("-inf"))),
                clean_report_text(item.get("target_id")),
            ),
        )
        best = entries_sorted[0]
        best_score = float(as_float(best.get("z_selected")) or 0.0)
        best_pct_rank = float(as_float(best.get("pct_rank")) or 1.0)
        best_target = clean_report_text(best.get("target_display")) or clean_report_text(
            best.get("pdb_id")
        )

        top_entries = entries_sorted[:top_targets_n]
        top_lines = []
        for idx, entry in enumerate(top_entries, start=1):
            target_display = clean_report_text(entry.get("pdb_id"))
            z_text = _format_num(entry.get("z_selected"), ".6g")
            p_text = _format_pct_display(
                entry.get("pct_rank"),
                max(2, pct_display_decimals),
                min_nonzero_pct=0.01,
            )
            top_lines.append(f"{idx}) {target_display}  z={z_text}  pct={p_text}")
        top_targets_text = "\n".join(top_lines) if top_lines else "—"

        tail_entries = entries_sorted[top_targets_n:]
        if tail_entries:
            tail_pcts = [float(e.get("pct_rank") or 0.0) for e in tail_entries]
            tail_median = median_floats(tail_pcts)
            if tail_median is None:
                tail_summary = f"+{len(tail_entries)} more"
            else:
                tail_pct = _format_pct_display(tail_median, tail_median_decimals)
                tail_summary = f"+{len(tail_entries)} more (median p={tail_pct})"
        else:
            tail_summary = "—"

        breadth_count = sum(
            1
            for entry in entries_sorted
            if float(as_float(entry.get("pct_rank")) or 1.0) <= breadth_pct_threshold
        )
        coverage_k = len(entries_sorted)
        ligand_name = ligand_display_map.get(lig_base, lig_base)
        ligand_highlights.append(
            {
                "ligand": ligand_name,
                "ligand_display": ligand_name,
                "ligand_base": lig_base,
                "best_target": best_target,
                "best_score": best_score,
                "best_percentile": _format_pct_display(
                    best_pct_rank,
                    max(2, pct_display_decimals),
                    min_nonzero_pct=0.01,
                ),
                "best_pct_rank": float(f"{best_pct_rank:.6f}"),
                "top_targets": top_targets_text,
                "coverage": f"{coverage_k} / {total_targets} targets",
                "coverage_k": coverage_k,
                "coverage_m": total_targets,
                "breadth": f"{breadth_count}",
                "breadth_count": breadth_count,
                "tail_summary": tail_summary,
            }
        )

    ligand_highlights.sort(
        key=lambda item: (
            float(as_float(item.get("best_pct_rank")) or 1.0),
            -(float(as_float(item.get("best_score")) or 0.0)),
            clean_report_text(item.get("ligand_base")),
        )
    )
    highlights_top_count = len(ligand_highlights)
    if highlights_top_count and highlights_top_pct > 0:
        pct_fraction = min(float(highlights_top_pct), 1.0)
        highlights_top_count = max(
            1, int(math.ceil(highlights_top_count * pct_fraction))
        )
    filtered_highlights = ligand_highlights[:highlights_top_count]
    if highlights_max_ligands > 0:
        ligand_highlights = filtered_highlights[:highlights_max_ligands]
    else:
        ligand_highlights = filtered_highlights

    sources: Dict[str, Any] = {
        "master_rows_csv": _safe_relpath(primary_master_csv, repo_root),
        "manifest_yaml": str(canonical_manifest_rel),
    }
    if mapping_csv is not None:
        sources["fda_mapping_csv"] = mapping_file_provenance_path(
            mapping_csv, repo_root
        )
        sources["fda_mapping_sha256"] = mapping_file_sha256(mapping_csv)
    if len(master_csvs) > 1:
        sources["master_rows_csvs"] = [
            _safe_relpath(path, repo_root) for path in master_csvs
        ]

    report = {
        "run_id": run_id,
        "generated_at": datetime.datetime.utcnow().isoformat() + "Z",
        "sources": sources,
        "summary": summary,
        "targets": targets_out,
        "ligands": ligands_out,
        "ligand_highlights": ligand_highlights,
        "ligand_highlight_config": {
            "top_targets_n": top_targets_n,
            "breadth_pct_threshold": float(f"{breadth_pct_threshold:.6f}"),
            "coverage_requires_pose_valid": bool(coverage_requires_pose_valid),
            "pct_display_decimals": pct_display_decimals,
            "tail_median_decimals": tail_median_decimals,
            "highlights_top_pct": float(f"{highlights_top_pct:.6f}"),
            "highlights_max_ligands": highlights_max_ligands,
        },
        "multi_target_hits": multi_target_hits,
    }

    return report

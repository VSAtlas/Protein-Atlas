# ruff: noqa: F403, F405
from analysis.reporting.run_report.constants import *
from analysis.reporting.run_report.config_resolve import _resolve_target_name_cfg
from analysis.reporting.run_report.display import (
    _resolve_ligand_display,
    _summarize_display_resolution,
)
from analysis.reporting.run_report.paths import _run_data_dir
from analysis.reporting.run_report.target_stats import (
    _resolve_target_name,
    _row_z_selected_rank_sort_key,
    _total_library_n,
)
def _write_heatmap_input_csv(
    repo_root: Path,
    run_id: str,
    out_path: Path,
    top_k: Optional[int],
    fda_mapping_csv: Optional[str],
    filter_invalid: bool,
    master_csv_override: Optional[Path] = None,
) -> None:
    logger = logging.getLogger("run-report")
    master_csv = master_csv_override or (
        _run_data_dir(repo_root, run_id) / "master_rows.csv"
    )
    if not master_csv.exists():
        raise FileNotFoundError(
            "Heatmap/report input check failed: master_rows CSV missing at "
            f"{master_csv}. Generate the run outputs (or pass an override CSV), "
            "then run `atlas --doctor` to verify paths."
        )

    existing_target_names: Dict[str, str] = {}
    if out_path.exists():
        try:
            with out_path.open("r", encoding="utf-8", newline="") as existing_handle:
                existing_reader = csv.DictReader(existing_handle)
                for existing_row in existing_reader:
                    target_id = clean_report_text(existing_row.get("target_id"))
                    target_name = clean_report_text(existing_row.get("target_name"))
                    pdb = clean_report_text(existing_row.get("pdb_id")).upper()
                    variant = clean_report_text(existing_row.get("variant"))
                    ph = clean_report_text(existing_row.get("ph_label"))
                    if target_name:
                        if target_id:
                            existing_target_names[target_id] = target_name
                        if pdb:
                            existing_target_names.setdefault(
                                build_target_id(pdb, variant, ph), target_name
                            )
                            existing_target_names.setdefault(pdb, target_name)
        except _REPORT_IO_ERRORS as exc:
            logger.warning(
                "%s action=target_name_backfill status=failed path=%s error=%s",
                COMPONENT,
                out_path,
                exc,
            )

    rows = []
    with master_csv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        has_ligand_display = "ligand_display" in (reader.fieldnames or [])
        for r in reader:
            rows.append(r)
    target_name_cfg = _resolve_target_name_cfg(repo_root, run_id)

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

    rows_by_target: Dict[str, List[Dict[str, Any]]] = {}
    row_rank: Dict[int, Dict[str, Any]] = {}
    for row in rows:
        pdb = clean_report_text(row.get("pdb_id"))
        variant = clean_report_text(row.get("variant"))
        ph = clean_report_text(row.get("ph_label"))
        target_id = build_target_id(pdb, variant, ph)
        rows_by_target.setdefault(target_id, []).append(row)

    for target_id, group_rows in rows_by_target.items():
        rankable: List[tuple[Dict[str, Any], float]] = []
        for row in group_rows:
            if as_report_row_bool(row.get("is_decoy")):
                continue
            if filter_invalid and not as_report_row_bool(row.get("pose_valid_any")):
                continue
            t_val = as_float(row.get("z_selected"))
            if t_val is None or not math.isfinite(t_val):
                continue
            rankable.append((row, t_val))

        sorted_rankable = sorted(rankable, key=_row_z_selected_rank_sort_key)
        total_library_n = _total_library_n(group_rows)
        for idx, (row, _t_val) in enumerate(sorted_rankable, start=1):
            pct = idx / total_library_n
            row_rank[id(row)] = {
                "rank": idx,
                "pct_rank": float(f"{pct:.6f}"),
                "target_id": target_id,
            }

    if top_k is not None and top_k > 0:
        filtered_rows = []
        for target_id, group_rows in rows_by_target.items():
            for row in group_rows:
                if filter_invalid and not as_report_row_bool(row.get("pose_valid_any")):
                    continue
                rank_info = row_rank.get(id(row))
                is_control = as_report_row_bool(row.get("is_control"))
                if is_control:
                    filtered_rows.append(row)
                elif rank_info and rank_info["rank"] <= top_k:
                    filtered_rows.append(row)
        rows = filtered_rows

    if filter_invalid and (top_k is None or top_k <= 0):
        rows = [row for row in rows if as_report_row_bool(row.get("pose_valid_any"))]

    def sort_key(row: Dict[str, Any]) -> tuple:
        rank_info = row_rank.get(id(row)) or {}
        rank_val = rank_info.get("rank")
        sort_rank = rank_val if isinstance(rank_val, int) else 1_000_000
        target_id = rank_info.get("target_id")
        if not target_id:
            pdb = clean_report_text(row.get("pdb_id"))
            variant = clean_report_text(row.get("variant"))
            ph = clean_report_text(row.get("ph_label"))
            target_id = build_target_id(pdb, variant, ph)
        lig_base = clean_report_text(row.get("ligand_base"))
        lig_file = clean_report_text(row.get("ligand_file") or row.get("ligand") or "")
        return (target_id, sort_rank, lig_base, lig_file)

    rows.sort(key=sort_key)

    fieldnames = [
        "target_id",
        "target_name",
        "pdb_id",
        "variant",
        "ph_label",
        "ligand_base",
        "ligand_display",
        "library",
        "z_selected",
        "z_selected_source",
        "fdr_score_field",
        "fdr_null_source",
        "fdr_n_decoys",
        "fdr_unique_decoy_scores",
        "fdr_reliable",
        "fdr_p_empirical",
        "fdr_q_target",
        "fdr_hit_q05",
        "fdr_hit_q10",
        "scorch_fdr_scope",
        "scorch_fdr_n_decoys",
        "scorch_fdr_n_tested",
        "scorch_fdr_q_bh",
        "scorch_fdr_decoy_competition_q",
        "scorch_fdr_decoy_competition_q_plus1",
        "pose_valid_any",
        "pose_invalid_reason_top",
        "rank",
        "pct_rank",
        "is_decoy",
        "is_control",
    ]

    _summarize_display_resolution(
        rows,
        has_ligand_display,
        fda_index,
        mapping_csv,
        logger,
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            rank_info = row_rank.get(id(row)) or {}
            pdb = clean_report_text(row.get("pdb_id"))
            variant = clean_report_text(row.get("variant"))
            ph = clean_report_text(row.get("ph_label"))
            target_id = build_target_id(pdb, variant, ph)
            display = _resolve_ligand_display(row, has_ligand_display, fda_index)
            resolved_target_name = (
                clean_report_text(row.get("target_name"))
                or _resolve_target_name(pdb, target_name_cfg)
                or existing_target_names.get(target_id)
                or existing_target_names.get(pdb.upper())
                or ""
            )
            writer.writerow(
                {
                    "target_id": target_id,
                    "target_name": resolved_target_name,
                    "pdb_id": pdb,
                    "variant": variant,
                    "ph_label": ph,
                    "ligand_base": clean_report_text(row.get("ligand_base")),
                    "ligand_display": display,
                    "library": clean_report_text(row.get("library")),
                    "z_selected": clean_report_text(row.get("z_selected")),
                    "z_selected_source": clean_report_text(row.get("z_selected_source")),
                    "fdr_score_field": clean_report_text(row.get("fdr_score_field")),
                    "fdr_null_source": clean_report_text(row.get("fdr_null_source")),
                    "fdr_n_decoys": clean_report_text(row.get("fdr_n_decoys")),
                    "fdr_unique_decoy_scores": clean_report_text(
                        row.get("fdr_unique_decoy_scores")
                    ),
                    "fdr_reliable": clean_report_text(row.get("fdr_reliable")),
                    "fdr_p_empirical": clean_report_text(row.get("fdr_p_empirical")),
                    "fdr_q_target": clean_report_text(row.get("fdr_q_target")),
                    "fdr_hit_q05": clean_report_text(row.get("fdr_hit_q05")),
                    "fdr_hit_q10": clean_report_text(row.get("fdr_hit_q10")),
                    "scorch_fdr_scope": clean_report_text(row.get("scorch_fdr_scope")),
                    "scorch_fdr_n_decoys": clean_report_text(row.get("scorch_fdr_n_decoys")),
                    "scorch_fdr_n_tested": clean_report_text(row.get("scorch_fdr_n_tested")),
                    "scorch_fdr_q_bh": clean_report_text(row.get("scorch_fdr_q_bh")),
                    "scorch_fdr_decoy_competition_q": clean_report_text(
                        row.get("scorch_fdr_decoy_competition_q")
                    ),
                    "scorch_fdr_decoy_competition_q_plus1": clean_report_text(
                        row.get("scorch_fdr_decoy_competition_q_plus1")
                    ),
                    "pose_valid_any": clean_report_text(row.get("pose_valid_any")),
                    "pose_invalid_reason_top": clean_report_text(
                        row.get("pose_invalid_reason_top")
                    ),
                    "rank": rank_info.get("rank", ""),
                    "pct_rank": rank_info.get("pct_rank", ""),
                    "is_decoy": clean_report_text(row.get("is_decoy")),
                    "is_control": clean_report_text(row.get("is_control")),
                }
            )


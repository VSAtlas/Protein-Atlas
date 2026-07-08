# ruff: noqa: F403, F405
from analysis.reporting.run_report.constants import *
from analysis.reporting.run_report.config_resolve import _resolve_target_name_cfg
from analysis.reporting.run_report.master_csv import (
    _is_decoy_filename,
)
from analysis.reporting.run_report.heatmap_csv import _write_heatmap_input_csv
from analysis.reporting.run_report.html_sections import _render_statistical_analysis_section
from analysis.reporting.run_report.paths import (
    _run_docked_dir,
)
from analysis.reporting.run_report.target_stats import _resolve_target_name
def _iter_vina_heatmap_targets(repo_root: Path, run_id: str) -> List[Tuple[str, str, str]]:
    targets: List[Tuple[str, str, str]] = []
    seen: Set[Tuple[str, str, str]] = set()
    manifest_data, _manifest_path = load_run_manifest(repo_root, run_id)
    proteins = (manifest_data or {}).get("proteins") if isinstance(manifest_data, dict) else None
    if isinstance(proteins, dict):
        for key, entry in proteins.items():
            pdb_id = ""
            variant = ""
            ph_label = ""
            key_text = str(key or "").strip()
            if "|" in key_text:
                parts = key_text.split("|", 2)
                if len(parts) >= 1:
                    pdb_id = clean_report_text(parts[0]).upper()
                if len(parts) >= 2:
                    variant = clean_report_text(parts[1])
                if len(parts) >= 3:
                    ph_label = clean_report_text(parts[2])
            if isinstance(entry, dict):
                if not pdb_id:
                    pdb_id = clean_report_text(entry.get("pdb_id")).upper()
                if not variant:
                    variant = clean_report_text(entry.get("variant"))
                if not ph_label:
                    ph_label = clean_report_text(entry.get("ph") or entry.get("ph_label"))
            if not pdb_id:
                continue
            if not variant:
                variant = "HOLO"
            item = (pdb_id, variant, ph_label)
            if item in seen:
                continue
            seen.add(item)
            targets.append(item)

    if targets:
        return sorted(targets)

    run_root = _run_docked_dir(repo_root, run_id)
    if not run_root.exists():
        return []
    for pdb_dir in sorted(run_root.iterdir()):
        if not pdb_dir.is_dir() or pdb_dir.name == "logs":
            continue
        pdb_id = clean_report_text(pdb_dir.name).upper()
        if not pdb_id:
            continue
        for variant_dir in sorted(pdb_dir.iterdir()):
            if not variant_dir.is_dir():
                continue
            variant = clean_report_text(variant_dir.name)
            if not variant:
                continue
            ph_dirs = [d for d in variant_dir.iterdir() if d.is_dir()]
            if not ph_dirs:
                item = (pdb_id, variant, "")
                if item not in seen:
                    seen.add(item)
                    targets.append(item)
                continue
            for ph_dir in sorted(ph_dirs):
                item = (pdb_id, variant, clean_report_text(ph_dir.name))
                if item in seen:
                    continue
                seen.add(item)
                targets.append(item)
    return sorted(targets)

def _strip_ligand_suffixes(raw_name: str) -> str:
    value = os.path.basename(str(raw_name or "").strip())
    if not value:
        return ""
    lowered = value.lower()
    suffixes = (
        ".pdbqt.gz",
        ".sdf.gz",
        ".mol2.gz",
        ".pdbqt",
        ".sdf",
        ".mol2",
        ".gz",
    )
    changed = True
    while changed and value:
        changed = False
        for suffix in suffixes:
            if lowered.endswith(suffix):
                value = value[: -len(suffix)]
                lowered = value.lower()
                changed = True
                break
    return clean_report_text(value) or clean_report_text(Path(raw_name).stem)

def _extract_vina_score_from_pdbqt(path: Path) -> Optional[float]:
    try:
        gz = path.suffix.lower() == ".gz"
        ctx = (
            gzip.open(path, "rt", encoding="utf-8", errors="ignore")
            if gz
            else path.open("r", encoding="utf-8", errors="ignore")
        )
        with ctx as handle:
            for idx, line in enumerate(handle):
                if idx >= 80:
                    break
                match = _VINA_RESULT_RE.search(line)
                if not match:
                    continue
                value = as_float(match.group(1))
                if value is not None and math.isfinite(value):
                    return float(value)
    except (OSError, UnicodeDecodeError):
        return None
    return None

def _strip_stage_suffix(name: str) -> str:
    text = clean_report_text(name)
    if not text:
        return ""
    stripped = re.sub(r"_stage\d+$", "", text, flags=re.IGNORECASE)
    stripped = re.sub(r"_(?:dud|fda_dud)$", "", stripped, flags=re.IGNORECASE)
    return clean_report_text(stripped) or text

def _collect_stage_scores(stage_dir: Path) -> Dict[str, float]:
    scores: Dict[str, float] = {}
    if not stage_dir.exists() or not stage_dir.is_dir():
        return scores
    pose_paths = sorted(stage_dir.glob("*.pdbqt")) + sorted(
        stage_dir.glob("*.pdbqt.gz")
    )
    for pose_path in pose_paths:
        stem = _strip_ligand_suffixes(pose_path.stem)
        ligand_base = _strip_stage_suffix(stem)
        if not ligand_base:
            continue
        score = _extract_vina_score_from_pdbqt(pose_path)
        if score is None:
            continue
        prev = scores.get(ligand_base)
        if prev is None or score < prev:
            scores[ligand_base] = score
    return scores

def _write_master_rows_csv_from_vina(
    repo_root: Path,
    run_id: str,
    out_path: Path,
    decoy_prefix: str,
    stage_required: str = "stage3",
) -> Dict[str, int]:
    logger = logging.getLogger("run-report")
    run_root = _run_docked_dir(repo_root, run_id)
    if not run_root.exists():
        raise FileNotFoundError(f"Docked run root not found: {run_root}")

    target_name_cfg = _resolve_target_name_cfg(repo_root, run_id)
    all_targets = _iter_vina_heatmap_targets(repo_root, run_id)

    stage_dirs = ("stage3", "stage2", "stage1")
    required_stage = clean_report_text(stage_required) or "stage3"
    stage_ready = 0
    rows: List[Dict[str, Any]] = []
    targets_with_scores = 0

    for pdb_id, variant, ph_label in all_targets:
        combo_dir = run_root / pdb_id / variant / ph_label
        if not combo_dir.exists():
            continue
        if not (combo_dir / required_stage).exists():
            continue
        stage_ready += 1

        per_target_scores: Dict[str, float] = {}
        for stage_name in stage_dirs:
            source_scores = _collect_stage_scores(combo_dir / stage_name)
            for ligand_base, score in source_scores.items():
                if ligand_base in per_target_scores:
                    continue
                per_target_scores[ligand_base] = score

        if not per_target_scores:
            continue
        targets_with_scores += 1

        target_id = build_target_id(pdb_id, variant, ph_label)
        target_name = _resolve_target_name(pdb_id, target_name_cfg)
        for ligand_base, score in sorted(per_target_scores.items()):
            if not ligand_base:
                continue
            is_decoy = _is_decoy_filename(ligand_base, decoy_prefix) or ligand_base.lower().startswith(
                "decoys_"
            )
            rows.append(
                {
                    "run_id": run_id,
                    "target_id": target_id,
                    "target_name": target_name,
                    "pdb_id": pdb_id,
                    "variant": variant,
                    "ph_label": ph_label,
                    "ligand_base": ligand_base,
                    "ligand_file": f"{ligand_base}.pdbqt",
                    "ligand_display": ligand_base,
                    "library": "decoy" if is_decoy else "fda",
                    "z_selected": -float(score),
                    "pose_valid_any": "1",
                    "pose_invalid_reason_top": "",
                    "is_decoy": "1" if is_decoy else "0",
                    "is_control": "0",
                }
            )

    if not rows:
        raise ValueError(
            (
                "No usable Vina heatmap rows found "
                f"(targets_seen={len(all_targets)} stage_ready={stage_ready} "
                f"targets_with_scores={targets_with_scores})"
            )
        )

    rows.sort(
        key=lambda row: (
            clean_report_text(row.get("target_id")),
            clean_report_text(row.get("ligand_base")),
        )
    )
    fieldnames = [
        "run_id",
        "target_id",
        "target_name",
        "pdb_id",
        "variant",
        "ph_label",
        "ligand_base",
        "ligand_file",
        "ligand_display",
        "library",
        "z_selected",
        "pose_valid_any",
        "pose_invalid_reason_top",
        "is_decoy",
        "is_control",
    ]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    logger.info(
        (
            "%s action=write_vina_master_rows status=ok path=%s "
            "targets_total=%d stage_ready=%d targets_with_scores=%d rows=%d"
        ),
        COMPONENT,
        out_path,
        len(all_targets),
        stage_ready,
        targets_with_scores,
        len(rows),
    )
    return {
        "targets_total": len(all_targets),
        "stage_ready": stage_ready,
        "targets_with_scores": targets_with_scores,
        "rows": len(rows),
    }

def _write_vina_fallback_fdr(
    master_csv: Path,
    summary_csv: Path,
    *,
    score_col: str = "z_selected",
) -> None:
    """Compute target-level decoy FDR for Vina fallback report exports."""

    if not master_csv.exists():
        return
    with master_csv.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = [dict(row) for row in reader]
        fieldnames = list(reader.fieldnames or [])
    if not rows or score_col not in fieldnames:
        return

    metric_fields = [
        "fdr_score_field",
        "fdr_null_source",
        "fdr_n_decoys",
        "fdr_unique_decoy_scores",
        "fdr_reliable",
        "fdr_p_empirical",
        "fdr_q_target",
        "fdr_hit_q05",
        "fdr_hit_q10",
        "fdr_decoy_competition_q",
        "fdr_decoy_competition_q_plus1",
    ]
    for field in metric_fields:
        if field not in fieldnames:
            fieldnames.append(field)
        for row in rows:
            row.setdefault(field, "")

    grouped: Dict[Tuple[str, str, str], List[int]] = {}
    for idx, row in enumerate(rows):
        grouped.setdefault(
            (
                clean_report_text(row.get("pdb_id")),
                clean_report_text(row.get("variant")),
                clean_report_text(row.get("ph_label")),
            ),
            [],
        ).append(idx)

    summary_rows: List[Dict[str, Any]] = []
    for key, indices in sorted(grouped.items()):
        decoys: List[float] = []
        tested: List[Tuple[int, float]] = []
        for idx in indices:
            score = as_float(rows[idx].get(score_col))
            if score is None or not math.isfinite(score):
                continue
            if as_report_row_bool(rows[idx].get("is_decoy")):
                decoys.append(score)
            else:
                tested.append((idx, score))

        result = compute_score_fdr(tested, decoys, higher_is_better=True)
        n_decoys = int(result.summary.get("n_decoys", 0))
        n_tested = int(result.summary.get("n_tested", 0))
        unique_decoys = len({round(score, 8) for score in decoys})
        reliable = n_decoys >= MIN_DECOYS_FOR_FDR
        for idx, metrics in result.row_metrics.items():
            q_bh = float(metrics.get("q_bh", 1.0))
            row = rows[idx]
            row["fdr_score_field"] = score_col
            row["fdr_null_source"] = "vina_fallback_decoy"
            row["fdr_n_decoys"] = str(n_decoys)
            row["fdr_unique_decoy_scores"] = str(unique_decoys)
            row["fdr_reliable"] = "1" if reliable else "0"
            row["fdr_p_empirical"] = f"{float(metrics.get('p_empirical', 1.0)):.6g}"
            row["fdr_q_target"] = f"{q_bh:.6g}"
            row["fdr_hit_q05"] = "1" if q_bh <= 0.05 else "0"
            row["fdr_hit_q10"] = "1" if q_bh <= 0.10 else "0"
            row["fdr_decoy_competition_q"] = (
                f"{float(metrics.get('decoy_competition_q', 1.0)):.6g}"
            )
            row["fdr_decoy_competition_q_plus1"] = (
                f"{float(metrics.get('decoy_competition_q_plus1', 1.0)):.6g}"
            )

        min_possible = float(result.summary.get("min_possible_bh_q", 1.0))
        summary_rows.append(
            {
                "pdb_id": key[0],
                "variant": key[1],
                "ph_label": key[2],
                "fdr_score_field": score_col,
                "fdr_null_source": "vina_fallback_decoy",
                "fdr_n_decoys": n_decoys,
                "fdr_unique_decoy_scores": unique_decoys,
                "fdr_reliable": 1 if reliable else 0,
                "n_tested": n_tested,
                "n_hits_q05": sum(rows[idx].get("fdr_hit_q05") == "1" for idx, _ in tested),
                "n_hits_q10": sum(rows[idx].get("fdr_hit_q10") == "1" for idx, _ in tested),
                "best_q": f"{float(result.summary.get('best_bh_q', 1.0)):.6g}",
                "min_possible_bh_q": f"{min_possible:.6g}",
                "bh_resolvable_q05": 1 if min_possible <= 0.05 else 0,
                "bh_resolvable_q10": 1 if min_possible <= 0.10 else 0,
                "n_decoy_competition_hits_q05": int(
                    result.summary.get("n_decoy_competition_hits_q05", 0)
                ),
                "n_decoy_competition_hits_q10": int(
                    result.summary.get("n_decoy_competition_hits_q10", 0)
                ),
                "best_decoy_competition_q": (
                    f"{float(result.summary.get('best_decoy_competition_q', 1.0)):.6g}"
                ),
            }
        )

    with master_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    summary_fieldnames = [
        "pdb_id",
        "variant",
        "ph_label",
        "fdr_score_field",
        "fdr_null_source",
        "fdr_n_decoys",
        "fdr_unique_decoy_scores",
        "fdr_reliable",
        "n_tested",
        "n_hits_q05",
        "n_hits_q10",
        "best_q",
        "min_possible_bh_q",
        "bh_resolvable_q05",
        "bh_resolvable_q10",
        "n_decoy_competition_hits_q05",
        "n_decoy_competition_hits_q10",
        "best_decoy_competition_q",
    ]
    summary_csv.parent.mkdir(parents=True, exist_ok=True)
    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

def _write_heatmap_input_csv_from_vina(
    repo_root: Path,
    run_id: str,
    out_path: Path,
    top_k: Optional[int],
    decoy_prefix: str,
    filter_invalid: bool,
    fda_mapping_csv: Optional[str],
) -> Dict[str, int]:
    master_out = out_path.with_name("master_rows_vina.csv")
    stats = _write_master_rows_csv_from_vina(
        repo_root=repo_root,
        run_id=run_id,
        out_path=master_out,
        decoy_prefix=decoy_prefix,
        stage_required="stage3",
    )
    _write_vina_fallback_fdr(master_out, out_path.with_name("target_fdr_summary.csv"))
    _write_heatmap_input_csv(
        repo_root=repo_root,
        run_id=run_id,
        out_path=out_path,
        top_k=top_k,
        fda_mapping_csv=fda_mapping_csv,
        filter_invalid=filter_invalid,
        master_csv_override=master_out,
    )
    return stats

def _write_vina_heatmap_html_report(
    run_id: str,
    out_path: Path,
    heatmap_html: str,
    repo_root: Path,
) -> None:
    title = f"{run_id} Vina pilot heatmap"
    statistical_html = _render_statistical_analysis_section(repo_root, out_path.parent)
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
        f"{statistical_html}\n"
        "</body>\n"
        "</html>\n"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as handle:
        handle.write(body)


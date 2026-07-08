# ruff: noqa: F403, F405
from analysis.reporting.run_report.constants import *
def _ligand_side_effect_cache_hit_count(repo_root: Path) -> tuple[int, bool]:
    from analysis.reporting.ligand_side_effect_cache import cache_path

    path = cache_path(repo_root)
    if not path.exists():
        return 0, False
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except _REPORT_MAIN_ERRORS:
        return 0, True
    entries = payload.get("entries") if isinstance(payload, dict) else {}
    if not isinstance(entries, dict):
        return 0, True
    hit_count = 0
    for entry in entries.values():
        if isinstance(entry, dict) and entry.get("side_effects"):
            hit_count += 1
    return hit_count, True

def _target_specs_from_heatmap_csv(repo_root: Path, heatmap_csv: Optional[Path]) -> list[Any]:
    if heatmap_csv is None or not heatmap_csv.exists():
        return []
    from analysis.reporting.target_safety_evidence import (
        TargetSpec,
        extract_local_uniprots_from_pdb,
    )

    targets: Dict[str, Any] = {}
    try:
        with heatmap_csv.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                pdb_id = clean_report_text(row.get("pdb_id")).upper()
                if not pdb_id:
                    continue
                target_name = clean_report_text(row.get("target_name"))
                current = targets.get(pdb_id)
                if current is None:
                    targets[pdb_id] = TargetSpec(
                        pdb_id=pdb_id,
                        target_name=target_name,
                        uniprots=tuple(extract_local_uniprots_from_pdb(repo_root, pdb_id)),
                        gene_symbols=(),
                    )
                elif target_name and not current.target_name:
                    targets[pdb_id] = TargetSpec(
                        pdb_id=current.pdb_id,
                        target_name=target_name,
                        uniprots=current.uniprots,
                        gene_symbols=current.gene_symbols,
                    )
    except _REPORT_IO_ERRORS:
        return []
    return sorted(targets.values(), key=lambda item: str(item.pdb_id).lower())

def _target_side_effect_cache_missing_pdbs(
    repo_root: Path, heatmap_csv: Optional[Path]
) -> tuple[list[Any], list[str], bool]:
    targets = _target_specs_from_heatmap_csv(repo_root, heatmap_csv)
    if not targets:
        return [], [], False
    path = repo_root / "pathways" / "cache" / "target_safety_aggregated.json"
    if not path.exists():
        return targets, [str(target.pdb_id) for target in targets], False
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except _REPORT_MAIN_ERRORS:
        return targets, [str(target.pdb_id) for target in targets], True
    entries = (payload.get("entries") or []) if isinstance(payload, dict) else []
    cached_pdbs = {
        clean_report_text(entry.get("pdb_id")).upper()
        for entry in entries
        if isinstance(entry, dict) and clean_report_text(entry.get("pdb_id"))
    }
    missing = [str(target.pdb_id) for target in targets if str(target.pdb_id).upper() not in cached_pdbs]
    return targets, missing, True

def _ensure_target_side_effect_cache(
    *,
    repo_root: Path,
    heatmap_csv: Optional[Path],
    mode: str,
    logger: logging.Logger,
) -> None:
    if mode == "never":
        return
    targets, missing_pdbs, cache_exists = _target_side_effect_cache_missing_pdbs(
        repo_root, heatmap_csv
    )
    if not targets:
        return
    if mode == "auto" and cache_exists and not missing_pdbs:
        logger.info(
            "%s action=target_side_effect_cache status=ready targets=%d",
            COMPONENT,
            len(targets),
        )
        return
    reason = "missing" if not cache_exists else f"missing_pdbs:{len(missing_pdbs)}"
    logger.warning(
        "%s action=target_side_effect_cache status=refresh reason=%s heatmap=%s",
        COMPONENT,
        reason,
        heatmap_csv,
    )
    try:
        from analysis.reporting.target_safety_drift import (
            build_target_safety_drift_summary,
            write_target_safety_drift_summary,
        )
        from analysis.reporting.target_safety_evidence import (
            build_target_safety_cache,
            write_target_safety_cache,
        )

        previous_payload: Dict[str, Any] = {}
        default_path = repo_root / "pathways" / "cache" / "target_safety_aggregated.json"
        if default_path.exists():
            try:
                previous_payload = json.loads(default_path.read_text(encoding="utf-8"))
            except _REPORT_MAIN_ERRORS:
                previous_payload = {}
        payload = build_target_safety_cache(repo_root, targets)
        out_path = write_target_safety_cache(repo_root, payload)
        drift_summary = build_target_safety_drift_summary(previous_payload, payload)
        drift_path = write_target_safety_drift_summary(repo_root, drift_summary)
        logger.info(
            "%s action=target_side_effect_cache status=refreshed targets=%d out=%s drift=%s",
            COMPONENT,
            len(payload.get("entries") or []),
            out_path,
            drift_path,
        )
    except Exception as exc:  # pragma: no cover - optional remote enrichment
        logger.warning(
            "%s action=target_side_effect_cache status=refresh_failed error=%s",
            COMPONENT,
            exc,
            exc_info=True,
        )

def _ensure_ligand_side_effect_cache(
    *,
    repo_root: Path,
    run_id: str,
    heatmap_csv: Optional[Path],
    fda_mapping_csv: Optional[str],
    mode: str,
    min_hits: int,
    limit: int,
    event_limit: int,
    label_limit: int,
    max_side_effects: int,
    sleep_sec: float,
    fetch_label_exposure: bool,
    logger: logging.Logger,
) -> None:
    if mode == "never" or heatmap_csv is None or not heatmap_csv.exists():
        return
    hit_count, cache_exists = _ligand_side_effect_cache_hit_count(repo_root)
    if mode == "auto" and cache_exists and hit_count >= max(0, min_hits):
        logger.info(
            "%s action=ligand_side_effect_cache status=ready hits=%d",
            COMPONENT,
            hit_count,
        )
        return
    reason = "missing" if not cache_exists else f"hits_below_min:{hit_count}<{min_hits}"
    logger.warning(
        "%s action=ligand_side_effect_cache status=refresh reason=%s heatmap=%s",
        COMPONENT,
        reason,
        heatmap_csv,
    )
    try:
        from analysis.cli.refresh_ligand_side_effect_cache import refresh_cache

        refresh_cache(
            repo_root=repo_root,
            run_id=run_id,
            heatmap_csv=heatmap_csv,
            fda_mapping_csv=fda_mapping_csv or "",
            limit=max(0, limit),
            event_limit=max(1, event_limit),
            label_limit=max(1, label_limit),
            max_side_effects=max(1, max_side_effects),
            sleep_sec=max(0.0, sleep_sec),
            reuse_existing=True,
            fda_only=True,
            include_controls=False,
            rdk_only=False,
            fetch_label_exposure=fetch_label_exposure,
            spd_supplement_xlsx="",
        )
    except Exception as exc:  # pragma: no cover - optional remote enrichment
        logger.warning(
            "%s action=ligand_side_effect_cache status=refresh_failed error=%s",
            COMPONENT,
            exc,
            exc_info=True,
        )

def _ensure_side_effect_caches(
    *,
    repo_root: Path,
    run_id: str,
    heatmap_csv: Optional[Path],
    fda_mapping_csv: Optional[str],
    mode: str,
    min_hits: int,
    limit: int,
    event_limit: int,
    label_limit: int,
    max_side_effects: int,
    sleep_sec: float,
    fetch_label_exposure: bool,
    logger: logging.Logger,
) -> None:
    _ensure_ligand_side_effect_cache(
        repo_root=repo_root,
        run_id=run_id,
        heatmap_csv=heatmap_csv,
        fda_mapping_csv=fda_mapping_csv,
        mode=mode,
        min_hits=min_hits,
        limit=limit,
        event_limit=event_limit,
        label_limit=label_limit,
        max_side_effects=max_side_effects,
        sleep_sec=sleep_sec,
        fetch_label_exposure=fetch_label_exposure,
        logger=logger,
    )
    _ensure_target_side_effect_cache(
        repo_root=repo_root,
        heatmap_csv=heatmap_csv,
        mode=mode,
        logger=logger,
    )


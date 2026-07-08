# ruff: noqa: F403, F405
from analysis.reporting.run_report.constants import *
from analysis.reporting.run_report.export_io import _resolve_heatmap_source
from analysis.reporting.run_report.paths import _run_data_dir
def _extract_pdb_id_from_target_id(target_id: Any) -> str:
    return pdb_id_from_target_id(target_id)

def _collect_pdb_ids_from_manifest(repo_root: Path, run_id: str) -> Set[str]:
    manifest, _manifest_path = load_run_manifest(repo_root, run_id)
    proteins = (
        (manifest or {}).get("proteins", {}) if isinstance(manifest, dict) else {}
    )
    if not isinstance(proteins, dict):
        return set()

    pdb_ids: Set[str] = set()
    for protein_key, entry in proteins.items():
        token = _extract_pdb_id_from_target_id(protein_key)
        if token:
            pdb_ids.add(token)
        if isinstance(entry, dict):
            entry_pdb_id = str(entry.get("pdb_id") or "").strip().upper()
            if entry_pdb_id:
                pdb_ids.add(entry_pdb_id)
    return pdb_ids

def _collect_pdb_ids_from_heatmap_csv(heatmap_csv: Path) -> Set[str]:
    if not heatmap_csv.exists():
        return set()
    pdb_ids: Set[str] = set()
    with heatmap_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        has_pdb_id = "pdb_id" in fieldnames
        for row in reader:
            token = str(row.get("pdb_id") or "").strip().upper() if has_pdb_id else ""
            if not token:
                token = _extract_pdb_id_from_target_id(row.get("target_id"))
            if token:
                pdb_ids.add(token)
    return pdb_ids

def _collect_pdb_ids_for_pathway_reports(
    repo_root: Path, run_id: str, heatmap_source: Path
) -> Set[str]:
    from_manifest = _collect_pdb_ids_from_manifest(repo_root, run_id)
    if from_manifest:
        return from_manifest
    if heatmap_source.suffix.lower() == ".csv":
        return _collect_pdb_ids_from_heatmap_csv(heatmap_source)
    fallback_csv = _run_data_dir(repo_root, run_id) / "heatmap_input.csv"
    return _collect_pdb_ids_from_heatmap_csv(fallback_csv)

def _render_pathway_heatmap_document(
    run_id: str, pathway_slug: str, body_html: str
) -> str:
    title = f"{run_id} pathway heatmap: {pathway_slug}"
    return (
        "<!doctype html>\n"
        "<html>\n"
        "<head>\n"
        '  <meta charset="utf-8">\n'
        '  <meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f"  <title>{html.escape(title)}</title>\n"
        "</head>\n"
        "<body>\n"
        f"{body_html}\n"
        "</body>\n"
        "</html>\n"
    )

def _write_pathway_specific_heatmap_reports(
    report: Dict[str, Any],
    html_path: Path,
    repo_root: Path,
    heatmap_source_override: Optional[Path] = None,
    report_asset_mode: str = "auto",
) -> List[Path]:
    run_id = str(report.get("run_id") or "").strip()
    if not run_id:
        return []

    heatmap_source = heatmap_source_override or _resolve_heatmap_source(
        repo_root, run_id
    )
    if heatmap_source is None:
        return []

    pdb_ids = sorted(
        _collect_pdb_ids_for_pathway_reports(repo_root, run_id, heatmap_source)
    )
    if not pdb_ids:
        return []

    logger = logging.getLogger("run-report")
    cache = pathway_resolver.Cache(
        cache_dir=repo_root / "pathways" / "cache",
        refresh=False,
        logger=logger,
    )
    http_client = pathway_resolver.HttpClient()
    pathway_by_pdb = pathway_resolver.assign_single_pathway_per_pdb(
        pdb_ids, cache, http_client
    )
    if not pathway_by_pdb:
        return []

    grouped: Dict[str, Set[str]] = {}
    for pdb_id in pdb_ids:
        slug = pathway_by_pdb.get(pdb_id) or pathway_resolver.slugify_pathway_name(
            "pathway"
        )
        grouped.setdefault(slug, set()).add(pdb_id)

    # Emit pathway reports even for single-target groups so pathway coverage
    # remains visible for sparse or highly specific benchmark runs.
    eligible_groups = {slug: pdb_set for slug, pdb_set in grouped.items() if pdb_set}
    if not eligible_groups:
        return []

    written_paths: List[Path] = []
    pathway_prefix = html_path.stem
    for pathway_slug in sorted(eligible_groups.keys()):
        pdb_filter = eligible_groups[pathway_slug]
        if not pdb_filter:
            continue
        pathway_path = html_path.with_name(f"{pathway_prefix}_{pathway_slug}.html")
        try:
            heatmap_html = render_interactive_heatmap_html(
                repo_root,
                run_id,
                heatmap_source,
                allowed_pdb_ids=pdb_filter,
                report_asset_mode=report_asset_mode,
            )
            html_payload = _render_pathway_heatmap_document(
                run_id, pathway_slug, heatmap_html
            )
        except _REPORT_HTML_BUILD_ERRORS as exc:
            html_payload = _render_pathway_heatmap_document(
                run_id,
                pathway_slug,
                ("<div>Heatmap unavailable: " f"{html.escape(str(exc))}" "</div>"),
            )
        pathway_path.write_text(html_payload, encoding="utf-8")
        written_paths.append(pathway_path)
    return written_paths

def _write_html_reports(
    report: Dict[str, Any],
    out_path: Path,
    highlight_queries: List[str],
    repo_root: Path,
    heatmap_source_override: Optional[Path] = None,
    report_asset_mode: str = "auto",
    report_html_warn_bytes: int = REPORT_HTML_WARN_BYTES_DEFAULT,
) -> List[Path]:
    _write_html_report(
        report,
        out_path,
        highlight_queries,
        repo_root,
        heatmap_source_override=heatmap_source_override,
        report_asset_mode=report_asset_mode,
        report_html_warn_bytes=report_html_warn_bytes,
    )
    return _write_pathway_specific_heatmap_reports(
        report,
        out_path,
        repo_root,
        heatmap_source_override=heatmap_source_override,
        report_asset_mode=report_asset_mode,
    )


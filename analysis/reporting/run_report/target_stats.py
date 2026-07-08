# ruff: noqa: F403, F405
from analysis.reporting.run_report.constants import *
from analysis.reporting.run_report.paths import _run_post_docked_dir
def _resolve_decoy_prefix(
    repo_root: Path, run_id: str, cli_value: Optional[str] = None
) -> str:
    if cli_value is not None:
        candidate = strip_quotes(str(cli_value))
        if candidate:
            return candidate

    value = resolve_config_key(repo_root, run_id, (DECOY_PREFIX_KEY, DUD_PREFIX_KEY))
    if value:
        return value
    return DECOY_PREFIX_DEFAULT

def _as_int(x: Any) -> int:
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return 0

def _pick_repeated_value(
    rows: List[Dict[str, Any]], col: str, is_int: bool = False
) -> Optional[Any]:
    for r in rows:
        val = r.get(col)
        if val is None:
            continue
        s = str(val).strip()
        if not s or s.lower() == "nan":
            continue
        try:
            f = float(s)
            if math.isfinite(f):
                return int(f) if is_int else f
        except (TypeError, ValueError):
            continue
    return None

def _load_scorch_stats(
    repo_root: Path,
    run_id: str,
    pdb_id: str,
    variant: str,
    ph: str,
    decoy_prefix: str = DECOY_PREFIX_DEFAULT,
) -> Dict[str, Optional[Any]]:
    stats: Dict[str, Optional[float]] = {"mu": None, "sigma": None, "n": None}
    combo_dir = _run_post_docked_dir(repo_root, run_id) / pdb_id / variant / ph

    candidates = [
        combo_dir / f"{decoy_prefix}_scorch_scores_all.csv",
        combo_dir / "dud_scorch_scores_all.csv",
        combo_dir / "scorch_scores_all.csv",
    ]

    for path in candidates:
        if path.exists() and path.stat().st_size > 0:
            try:
                with path.open("r", encoding="utf-8") as f:
                    reader = csv.DictReader(f)
                    # Just need first row with data
                    for row in reader:
                        mu = as_float(row.get("scorch_mu_decoy"))
                        sigma = as_float(row.get("scorch_sigma_decoy"))
                        n = _as_int(row.get("scorch_n_decoys"))
                        if mu is not None or sigma is not None or n > 0:
                            stats["mu"] = mu
                            stats["sigma"] = sigma
                            stats["n"] = n if n > 0 else None
                            return stats
            except _REPORT_IO_ERRORS:
                continue
    return stats

def _resolve_target_name(pdb_id: str, cfg: Dict[str, Any]) -> str:
    normalized_pdb = clean_report_text(pdb_id).upper()
    if not normalized_pdb:
        return ""
    cached = _TARGET_NAME_CACHE.get(normalized_pdb)
    if cached is not None:
        return cached

    prefer_raw = clean_report_text(
        cfg.get("target_name_prefer")
        or cfg.get("TARGET_NAME_PREFER")
        or cfg.get("target_name_preference")
        or cfg.get("TARGET_NAME_PREFERENCE")
        or "auto"
    ).lower()
    prefer = prefer_raw if prefer_raw in {"auto", "compnd", "uniprot"} else "auto"
    repo_root = (
        Path(clean_report_text(cfg.get("repo_root"))).resolve()
        if clean_report_text(cfg.get("repo_root"))
        else None
    )
    try:
        if repo_root and repo_root.exists():
            cwd = Path.cwd()
            os.chdir(repo_root)
            try:
                target_name = clean_report_text(
                    derive_target_name(
                        normalized_pdb,
                        prefer=prefer,
                        pdb_root_override=repo_root / "input_pdbs",
                        cfg=cfg,
                    )
                )
            finally:
                os.chdir(cwd)
        else:
            target_name = clean_report_text(
                derive_target_name(
                    normalized_pdb,
                    prefer=prefer,
                    pdb_root_override=Path("input_pdbs"),
                    cfg=cfg,
                )
            )
    except (OSError, RuntimeError, ValueError, KeyError) as exc:
        logging.getLogger(__name__).debug(
            "%s derive_target_name_failed pdb=%s err=%s",
            COMPONENT,
            normalized_pdb,
            exc,
            exc_info=True,
        )
        target_name = ""

    _TARGET_NAME_CACHE[normalized_pdb] = target_name
    return target_name

def _ligand_identity_for_library(row: Dict[str, Any]) -> str:
    lig_base = clean_report_text(row.get("ligand_base"))
    if lig_base:
        return lig_base
    for key in ("ligand", "ligand_file"):
        raw = clean_report_text(row.get(key))
        if not raw:
            continue
        basename = os.path.basename(raw)
        stem = Path(basename).stem if basename else ""
        token = clean_report_text(stem or basename)
        if token:
            return token
    return ""

def _total_library_n(group_rows: List[Dict[str, Any]]) -> int:
    unique_ligands: Set[str] = set()
    for row in group_rows:
        if as_report_row_bool(row.get("is_decoy")):
            continue
        ligand_id = _ligand_identity_for_library(row)
        if ligand_id:
            unique_ligands.add(ligand_id)
    return max(1, len(unique_ligands))

def _row_z_selected_rank_sort_key(item: tuple[Dict[str, Any], float]) -> tuple[float, str, str]:
    row, t_val = item
    lig_base = clean_report_text(row.get("ligand_base"))
    lig_file = clean_report_text(row.get("ligand_file") or row.get("ligand") or "")
    return (-t_val, lig_base, lig_file)

def _target_display_name(pdb_id: str, target_name: str) -> str:
    pdb = clean_report_text(pdb_id)
    name = clean_report_text(target_name)
    if name and pdb:
        return f"{name} ({pdb})"
    if pdb:
        return pdb
    return name

def _format_pct_display(
    value: Any, decimals: int = 1, min_nonzero_pct: Optional[float] = None
) -> str:
    pct = as_float(value)
    if pct is None or not math.isfinite(pct):
        return ""
    pct_display = pct * 100.0
    if (
        min_nonzero_pct is not None
        and min_nonzero_pct > 0
        and pct_display > 0
        and pct_display < min_nonzero_pct
    ):
        pct_display = min_nonzero_pct
    digits = max(0, int(decimals))
    return f"{pct_display:.{digits}f}%"

def _compute_fdr_stats(group_rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    non_decoys = [r for r in group_rows if not as_report_row_bool(r.get("is_decoy"))]
    if not non_decoys:
        return {}

    score_field = None
    for r in non_decoys:
        val = (r.get("fdr_score_field") or "").strip()
        if val:
            score_field = val
            break

    n_decoys = None
    unique_decoys = None
    for r in non_decoys:
        dec = as_float(r.get("fdr_n_decoys"))
        if dec is not None:
            n_decoys = int(dec)
            break
    for r in non_decoys:
        uniq = as_float(r.get("fdr_unique_decoy_scores"))
        if uniq is not None:
            unique_decoys = int(uniq)
            break

    n_tested = 0
    n_hits_q05 = 0
    n_hits_q10 = 0
    best_q = None

    for r in non_decoys:
        q_val = as_float(r.get("fdr_q_target"))
        if q_val is None:
            continue
        n_tested += 1
        best_q = q_val if best_q is None else min(best_q, q_val)
        if q_val <= 0.05:
            n_hits_q05 += 1
        if q_val <= 0.10:
            n_hits_q10 += 1

    if best_q is None:
        best_q = 1.0

    reliable = (n_decoys is not None and n_decoys >= MIN_DECOYS_FOR_FDR) and (
        unique_decoys is not None and unique_decoys >= MIN_UNIQUE_DECOY_SCORES
    )
    return {
        "score_field": score_field,
        "n_decoys": n_decoys,
        "n_tested": n_tested,
        "n_hits_q05": n_hits_q05,
        "n_hits_q10": n_hits_q10,
        "best_q": best_q,
        "unique_decoy_scores": unique_decoys,
        "reliable": reliable,
    }

def _vector_from_row(row: Dict[str, Any], prefix: str) -> Optional[List[float]]:
    vx = as_float(row.get(f"{prefix}_x"))
    vy = as_float(row.get(f"{prefix}_y"))
    vz = as_float(row.get(f"{prefix}_z"))
    if (
        vx is not None
        and vy is not None
        and vz is not None
        and math.isfinite(vx)
        and math.isfinite(vy)
        and math.isfinite(vz)
    ):
        return [float(vx), float(vy), float(vz)]
    return None

def _resolve_manifest_pocket(
    manifest: Optional[Dict[str, Any]], pdb: str, variant: str, ph: str
) -> Dict[str, Any]:
    pocket = extract_pocket(manifest, pdb, variant, ph)
    center = pocket.get("center")
    box = pocket.get("box")

    normalized_center = None
    normalized_box = None
    if isinstance(center, (list, tuple)) and len(center) == 3:
        cx = as_float(center[0])
        cy = as_float(center[1])
        cz = as_float(center[2])
        if (
            cx is not None
            and cy is not None
            and cz is not None
            and math.isfinite(cx)
            and math.isfinite(cy)
            and math.isfinite(cz)
        ):
            normalized_center = [cx, cy, cz]

    if isinstance(box, (list, tuple)) and len(box) == 3:
        bx = as_float(box[0])
        by = as_float(box[1])
        bz = as_float(box[2])
        if (
            bx is not None
            and by is not None
            and bz is not None
            and math.isfinite(bx)
            and math.isfinite(by)
            and math.isfinite(bz)
        ):
            normalized_box = [bx, by, bz]

    return {
        "method": pocket.get("method"),
        "center": normalized_center,
        "box": normalized_box,
    }


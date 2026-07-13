# ruff: noqa: F403, F405
from analysis.reporting.run_report.constants import *
from analysis.reporting.run_report.paths import _run_data_dir
def _infer_decoy_prefix_from_tokens(decoy_prefix: str, tokens: List[str]) -> str:
    if decoy_prefix and decoy_prefix != DECOY_PREFIX_DEFAULT:
        return decoy_prefix
    if "dud" in tokens:
        return "dud"
    for tok in tokens:
        if "dud" in tok:
            return tok
    return decoy_prefix

def _sanitize_token(token: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_]+", "_", str(token))
    cleaned = cleaned.strip("_")
    return cleaned or "custom"

def _resolve_master_csvs(
    repo_root: Path, run_id: str, tokens: List[str], decoy_prefix: str
) -> List[Path]:
    data_dir = _run_data_dir(repo_root, run_id)
    primary = data_dir / "master_rows.csv"
    if primary.exists():
        return [primary]

    candidates: List[Path] = []
    for tok in tokens:
        if tok == "dud" or tok == decoy_prefix:
            continue
        prefix = _sanitize_token(tok)
        candidates.append(data_dir / f"{prefix}_master_rows.csv")
        candidates.append(data_dir / f"master_rows_{prefix}.csv")
    return [path for path in candidates if path.exists()]

def _safe_relpath(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)

def _is_decoy_filename(name: str, decoy_prefix: str) -> bool:
    token = str(name or "").strip().lower()
    if not token:
        return False
    if token.startswith("decoy") or token.startswith("dud") or "_dud_" in token:
        return True
    decoy = str(decoy_prefix or "").strip().lower()
    if decoy and (
        token.startswith(f"{decoy}_")
        or f"_{decoy}_" in token
        or token.startswith(decoy)
    ):
        return True
    return False

def _row_is_decoy_for_split(row: Dict[str, Any], decoy_prefix: str) -> bool:
    if row_is_control(row):
        return False
    explicit_decoy = str(row.get("is_decoy", "")).strip()
    if explicit_decoy:
        return row_has_explicit_decoy(row)
    if row_is_decoy_by_role(row):
        return True
    ligand_file = clean_report_text(ligand_filename_from_row(row))
    return _is_decoy_filename(ligand_file, decoy_prefix)

def _normalize_library_slug(value: Any) -> str:
    text = clean_report_text(value).lower()
    if not text:
        return ""
    text = re.sub(r"[^a-z0-9]+", "_", text).strip("_")
    if not text:
        return ""
    if text in {"na", "none", "null", "unknown"}:
        return ""
    return text

def _infer_library_from_source_csv(source_csv: Any, decoy_prefix: str) -> str:
    source = clean_report_text(source_csv)
    if not source:
        return ""
    basename = Path(source).name.lower()
    if basename == "consensus_reranked_scorch.csv":
        return "fda"
    suffix = "_consensus_reranked_scorch.csv"
    if not basename.endswith(suffix):
        return ""
    prefix = basename[: -len(suffix)]
    if not prefix:
        return ""
    decoy = _normalize_library_slug(decoy_prefix)
    if prefix in {"dud", "decoy"} or (decoy and prefix == decoy):
        return ""
    return _normalize_library_slug(prefix)

def _resolve_row_library_slug(row: Dict[str, Any], decoy_prefix: str) -> str:
    if _row_is_decoy_for_split(row, decoy_prefix):
        return "decoy"
    lib = _normalize_library_slug(row.get("library"))
    if lib and lib not in {"decoy", "dud"}:
        return lib
    mode = _normalize_library_slug(row.get("run_mode"))
    if mode and mode not in {"decoy", "dud"}:
        return mode
    source_lib = _infer_library_from_source_csv(row.get("source_csv"), decoy_prefix)
    if source_lib:
        return source_lib
    return "unknown"

def _write_split_master_rows(
    repo_root: Path,
    run_id: str,
    decoy_prefix: str,
    master_csvs: List[Path],
) -> Dict[str, Path]:
    data_dir = _run_data_dir(repo_root, run_id)
    rows: List[Dict[str, Any]] = []
    fieldnames: List[str] = []
    field_set: Set[str] = set()
    for csv_path in master_csvs:
        with csv_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for name in reader.fieldnames or []:
                if name not in field_set:
                    field_set.add(name)
                    fieldnames.append(name)
            for row in reader:
                rows.append(row)
    if not rows:
        return {}

    decoy_rows: List[Dict[str, Any]] = []
    by_library: Dict[str, List[Dict[str, Any]]] = {}
    for row in rows:
        lib_slug = _resolve_row_library_slug(row, decoy_prefix)
        if lib_slug == "decoy":
            decoy_rows.append(row)
            continue
        by_library.setdefault(lib_slug or "unknown", []).append(row)

    written: Dict[str, Path] = {}
    for lib_slug in sorted(by_library.keys()):
        subset_rows = list(by_library[lib_slug]) + list(decoy_rows)
        out_path = data_dir / f"{lib_slug}_master_rows.csv"
        with out_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(subset_rows)
        written[lib_slug] = out_path
    return written


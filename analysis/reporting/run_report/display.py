# ruff: noqa: F403, F405
from analysis.reporting.run_report.constants import *
def _norm_text(value: Any) -> str:
    text = clean_report_text(value).lower()
    if not text:
        return ""
    text = re.sub(r"[_-]+", " ", text)
    text = re.sub(r"\s+", " ", text).strip()
    text = re.sub(r"^[\W_]+|[\W_]+$", "", text)
    return re.sub(r"\s+", " ", text).strip()

def _is_iupac_like(text: str) -> bool:
    value = str(text or "").strip()
    if not value:
        return False
    length = len(value)
    if length < 25:
        return False
    digits = sum(ch.isdigit() for ch in value)
    punct = sum(1 for ch in value if not ch.isalnum() and not ch.isspace())
    if length >= 50:
        return True
    if length >= 35 and (punct / length) >= 0.12:
        return True
    if length >= 35 and (digits / length) >= 0.2:
        return True
    if length >= 30 and punct >= 6 and (digits / length) >= 0.1:
        return True
    return False

def _normalize_inchikey_candidate(text: str) -> str:
    value = str(text or "").strip()
    if not value:
        return ""
    value = re.sub(r"\s+", "-", value)
    value = re.sub(r"-+", "-", value.upper())
    return value

def _is_inchikey_like(text: str) -> bool:
    candidate = _normalize_inchikey_candidate(text)
    return bool(candidate and _INCHIKEY_RE.match(candidate))

def _is_relaxed_iupac_like(text: str) -> bool:
    value = str(text or "").strip()
    if not value:
        return False
    if value != value.lower():
        return False
    if len(value) < 18:
        return False
    digits = sum(ch.isdigit() for ch in value)
    if digits < 2:
        return False
    if value.count(" ") < 4:
        return False
    return True

def _is_unfriendly_display(text: str) -> bool:
    value = clean_report_text(text)
    if not value:
        return True
    if _is_inchikey_like(value):
        return True
    if _RDK_PLACEHOLDER_RE.search(value):
        return True
    if _is_iupac_like(value) or _is_relaxed_iupac_like(value):
        return True
    return False

def _parse_highlight_queries(value: Optional[str]) -> List[str]:
    if not value:
        return []
    queries = []
    for part in str(value).split(","):
        q = part.strip()
        if q:
            queries.append(q)
    return queries

def _matches_query(
    query: str,
    ligand_display: str,
    ligand_base: str,
    ligand_raw: str,
    mode: str,
    fda_index: Optional[Any],
) -> bool:
    q_norm = _norm_text(query)
    if not q_norm:
        return False
    disp_norm = _norm_text(ligand_display)
    base_norm = _norm_text(ligand_base)
    raw_norm = _norm_text(ligand_raw)

    if mode == "display_exact":
        if disp_norm and disp_norm == q_norm:
            return True
    elif mode == "display_contains":
        if disp_norm and q_norm in disp_norm:
            return True
    else:
        if (q_norm in disp_norm) or (q_norm in base_norm) or (q_norm in raw_norm):
            return True

    if fda_index is None:
        return False

    query_alias = clean_report_text(resolve_ligand_display_name(query, "", fda_index))
    query_alias_norm = _norm_text(query_alias)
    if not query_alias_norm:
        return False

    for text in (ligand_display, ligand_base, ligand_raw):
        alias = clean_report_text(resolve_ligand_display_name(text, "", fda_index))
        if alias and _norm_text(alias) == query_alias_norm:
            return True

    return False

def _resolve_ligand_display(
    row: Dict[str, Any], has_ligand_display: bool, fda_index: Optional[Any]
) -> str:
    display, _meta = _resolve_ligand_display_internal(
        row, has_ligand_display, fda_index
    )
    return display

def _resolve_ligand_display_internal(
    row: Dict[str, Any], has_ligand_display: bool, fda_index: Optional[Any]
) -> Tuple[str, Dict[str, Any]]:
    base = clean_report_text(row.get("ligand_base"))
    lig_file = clean_report_text(row.get("ligand_file") or row.get("ligand") or "")
    display = clean_report_text(row.get("ligand_display")) if has_ligand_display else ""
    is_control = as_report_row_bool(row.get("is_control"))
    mapping_used = False

    def build_meta(final_display: str) -> Dict[str, Any]:
        original_unfriendly = _is_unfriendly_display(display)
        final_unfriendly = _is_unfriendly_display(final_display)
        mapping_override = bool(
            has_ligand_display
            and not is_control
            and mapping_used
            and final_display
            and final_display != display
        )
        return {
            "is_control": is_control,
            "has_ligand_display": has_ligand_display,
            "original_display": display,
            "final_display": final_display,
            "original_unfriendly": original_unfriendly,
            "final_unfriendly": final_unfriendly,
            "mapping_used": mapping_used,
            "mapping_override": mapping_override,
        }

    if is_control:
        if display:
            return display, build_meta(display)
        fallback = clean_report_text(row.get("ligand"))
        final_display = fallback or base
        return final_display, build_meta(final_display)

    authoritative_display = clean_report_text(
        resolve_authoritative_ligand_display_name(base, lig_file, fda_index)
    )
    if authoritative_display:
        mapping_used = authoritative_display != display
        return authoritative_display, build_meta(authoritative_display)

    if display:
        if fda_index is not None and _is_unfriendly_display(display):
            alt = clean_report_text(resolve_ligand_display_name(base, lig_file, fda_index))
            if alt and alt != display and not _is_unfriendly_display(alt):
                mapping_used = True
                return alt, build_meta(alt)
        if not _is_unfriendly_display(display):
            return display, build_meta(display)
        if base and not _is_unfriendly_display(base):
            return base, build_meta(base)
        return display, build_meta(display)

    if fda_index is not None:
        alt = clean_report_text(resolve_ligand_display_name(base, lig_file, fda_index))
        if alt and not _is_unfriendly_display(alt):
            mapping_used = True
            return alt, build_meta(alt)

    fallback = clean_report_text(row.get("ligand"))
    final_display = fallback or base
    return final_display, build_meta(final_display)

def _summarize_display_resolution(
    rows: List[Dict[str, Any]],
    has_ligand_display: bool,
    fda_index: Optional[Any],
    mapping_csv: Optional[Path],
    logger: logging.Logger,
    top_n: int = 15,
) -> None:
    if not rows:
        return
    mapping_csv_value = str(mapping_csv) if mapping_csv else "missing"
    mapping_loaded = bool(fda_index)
    n_overridden = 0
    n_unfriendly_remaining = 0
    examples: List[Tuple[str, str, str, str, str]] = []
    seen: Set[Tuple[str, str, str, str, str]] = set()

    for row in rows:
        final_display, meta = _resolve_ligand_display_internal(
            row, has_ligand_display, fda_index
        )
        if meta["is_control"]:
            continue
        if meta["mapping_override"]:
            n_overridden += 1
        if meta["final_unfriendly"]:
            n_unfriendly_remaining += 1
            if len(examples) < top_n:
                lig_base = clean_report_text(row.get("ligand_base"))
                lig_file = clean_report_text(
                    row.get("ligand_file") or row.get("ligand") or ""
                )
                lig_file = os.path.basename(lig_file)
                library = clean_report_text(row.get("library"))
                key = (
                    lig_base,
                    lig_file,
                    library,
                    meta["original_display"],
                    final_display,
                )
                if key not in seen:
                    seen.add(key)
                    examples.append(key)

    logger.info(
        "%s action=fda_name_override overridden=%d unfriendly_remaining=%d mapping_csv=%s mapping_loaded=%s",
        COMPONENT,
        n_overridden,
        n_unfriendly_remaining,
        mapping_csv_value,
        mapping_loaded,
    )
    for lig_base, lig_file, library, original_display, final_display in examples:
        logger.info(
            "%s action=fda_name_unresolved ligand_base=%s ligand_file=%s library=%s original_display=%s final_display=%s mapping_csv=%s mapping_loaded=%s",
            COMPONENT,
            lig_base,
            lig_file,
            library,
            original_display,
            final_display,
            mapping_csv_value,
            mapping_loaded,
        )

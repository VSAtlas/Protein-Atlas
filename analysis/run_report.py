# -*- coding: utf-8 -*-
import argparse
import csv
import datetime
import html
import logging
import math
import os
import re
import shutil
import sys
import yaml  # type: ignore[import-untyped]
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

try:
    from analysis.fda_name_map import (
        resolve_ligand_display_name,
        resolve_mapping_csv_path,
        try_load_fda_index,
    )
    from analysis.dud_eval_discovery import derive_target_name
    from analysis.heatmap_html import render_interactive_heatmap_html
    from analysis.manifest_utils import extract_pocket, load_run_manifest
    import pathway_resolver
except ModuleNotFoundError:
    REPO_ROOT = Path(__file__).resolve().parents[1]
    if str(REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(REPO_ROOT))
    SRC_ROOT = REPO_ROOT / "src"
    if str(SRC_ROOT) not in sys.path:
        sys.path.insert(0, str(SRC_ROOT))
    from analysis.fda_name_map import (
        resolve_ligand_display_name,
        resolve_mapping_csv_path,
        try_load_fda_index,
    )
    from analysis.dud_eval_discovery import derive_target_name
    from analysis.heatmap_html import render_interactive_heatmap_html
    from analysis.manifest_utils import extract_pocket, load_run_manifest
    import pathway_resolver

COMPONENT = "[run-report]"
DECOY_PREFIX_KEY = "DECOY_PREFIX"
DUD_PREFIX_KEY = "DUD_PREFIX"
DECOY_PREFIX_DEFAULT = "dud"
_TEST_MODE_OFF_VALUES = {"", "0", "false", "no", "off", "none", "null"}
_TEST_MODE_ON_VALUES = {"true", "yes", "on", "1"}
_TEST_MODE_BOTH_VALUES = {
    "both",
    "fda_dud",
    "dud_fda",
    "fda+dud",
    "dud+fda",
    "fda-dud",
    "dud-fda",
}
MIN_DECOYS_FOR_FDR = 200
MIN_UNIQUE_DECOY_SCORES = 10
_INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_RDK_PLACEHOLDER_RE = re.compile(r"rdk[_-]?\d+", re.IGNORECASE)
_TARGET_NAME_CACHE: Dict[str, str] = {}


def _configure_logging(verbose: bool) -> logging.Logger:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(message)s")
    return logging.getLogger("run-report")


def _strip_quotes(value: str) -> str:
    stripped = str(value or "").strip()
    if (
        len(stripped) >= 2
        and stripped[0] == stripped[-1]
        and stripped[0] in ("'", '"')
    ):
        return stripped[1:-1].strip()
    return stripped


def _read_test_mode_from_file(path: Path) -> Optional[str]:
    if not path.exists():
        return None
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("#") or "=" not in stripped:
                    continue
                key, raw_value = stripped.split("=", 1)
                if key.strip() != "TEST_MODE_ENABLE":
                    continue
                return _strip_quotes(raw_value)
    except Exception:
        return None
    return None


def _parse_percent_searched(raw: str) -> Optional[float]:
    value = _strip_quotes(str(raw or "")).strip()
    if not value:
        return None
    scale = 1.0
    if value.endswith("%"):
        value = value[:-1].strip()
        scale = 0.01
    try:
        parsed = float(value) * scale
    except ValueError:
        return None
    if not math.isfinite(parsed) or parsed < 0:
        return None
    if parsed > 1 and parsed <= 100 and scale == 1.0:
        parsed = parsed / 100.0
    if parsed > 1:
        return None
    return parsed


def _read_percent_searched_from_file(path: Path) -> Optional[float]:
    if not path.exists():
        return None
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("#") or "=" not in stripped:
                    continue
                key, raw_value = stripped.split("=", 1)
                if key.strip() not in {"PERCENT_SEARCHED", "percent_searched"}:
                    continue
                parsed = _parse_percent_searched(raw_value)
                if parsed is not None:
                    return parsed
    except Exception:
        return None
    return None


def _resolve_percent_searched(repo_root: Path) -> Optional[float]:
    return _read_percent_searched_from_file(repo_root / "config.txt")


def _parse_filter_invalid(raw: object) -> Optional[bool]:
    if isinstance(raw, bool):
        return raw
    value = _strip_quotes(str(raw or "")).strip()
    if not value:
        return None
    lowered = value.lower()
    if lowered in _TEST_MODE_ON_VALUES:
        return True
    if lowered in _TEST_MODE_OFF_VALUES:
        return False
    return None


def _read_filter_invalid_from_file(path: Path) -> Optional[bool]:
    if not path.exists():
        return None
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("#") or "=" not in stripped:
                    continue
                key, raw_value = stripped.split("=", 1)
                if key.strip().lower() != "filter_invalid":
                    continue
                parsed = _parse_filter_invalid(raw_value)
                if parsed is not None:
                    return parsed
    except Exception:
        return None
    return None


def _resolve_filter_invalid(repo_root: Path, run_id: str) -> Optional[bool]:
    env_raw = os.environ.get("FILTER_INVALID")
    if env_raw is not None:
        parsed = _parse_filter_invalid(env_raw)
        if parsed is not None:
            return parsed
    candidates = [
        repo_root / "data" / run_id / "config.txt",
        repo_root / "data" / run_id / "config_snapshot.txt",
        repo_root / "manifests" / run_id / "config.txt",
        repo_root / "config.txt",
    ]
    for path in candidates:
        parsed = _read_filter_invalid_from_file(path)
        if parsed is not None:
            return parsed
    return None


def _read_target_name_prefer_from_file(path: Path) -> Optional[str]:
    if not path.exists():
        return None
    prefer_keys = {
        "TARGET_NAME_PREFER",
        "TARGET_NAME_PREFERENCE",
        "DUD_EVAL_TARGET_NAME_PREFER",
        "DUD_EVAL_TARGET_NAME_PREFERENCE",
    }
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("#") or "=" not in stripped:
                    continue
                key, raw_value = stripped.split("=", 1)
                key = key.strip().upper()
                if key not in prefer_keys:
                    continue
                value = _strip_quotes(raw_value).strip().lower()
                if value:
                    return value
    except Exception:
        return None
    return None


def _resolve_target_name_cfg(repo_root: Path, run_id: str) -> Dict[str, str]:
    env_value = _strip_quotes(os.environ.get("TARGET_NAME_PREFER", "")).strip().lower()
    if env_value:
        return {"target_name_prefer": env_value, "repo_root": str(repo_root)}

    candidates = [
        repo_root / "data" / run_id / "config.txt",
        repo_root / "data" / run_id / "config_snapshot.txt",
        repo_root / "manifests" / run_id / "config.txt",
        repo_root / "config.txt",
    ]
    for path in candidates:
        value = _read_target_name_prefer_from_file(path)
        if value:
            return {"target_name_prefer": value, "repo_root": str(repo_root)}
    return {"repo_root": str(repo_root)}


def _parse_report_limit(raw: object) -> Optional[int]:
    if isinstance(raw, bool):
        return None
    if isinstance(raw, int):
        return raw if raw >= 0 else None
    value = _strip_quotes(str(raw or "")).strip()
    if not value:
        return None
    try:
        parsed = int(value)
    except ValueError:
        try:
            parsed_float = float(value)
        except ValueError:
            return None
        if not math.isfinite(parsed_float) or parsed_float < 0:
            return None
        if not float(parsed_float).is_integer():
            return None
        parsed = int(parsed_float)
    if parsed < 0:
        return None
    return parsed


def _read_report_limit_from_file(path: Path, key: str) -> Optional[int]:
    if not path.exists():
        return None
    key_norm = key.strip().lower()
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("#") or "=" not in stripped:
                    continue
                k, raw_value = stripped.split("=", 1)
                if k.strip().lower() != key_norm:
                    continue
                parsed = _parse_report_limit(raw_value)
                if parsed is not None:
                    return parsed
    except Exception:
        return None
    return None


def _resolve_report_limit(repo_root: Path, run_id: str, key: str) -> Optional[int]:
    candidates = [
        repo_root / "data" / run_id / "config.txt",
        repo_root / "data" / run_id / "config_snapshot.txt",
        repo_root / "manifests" / run_id / "config.txt",
        repo_root / "config.txt",
    ]
    for path in candidates:
        parsed = _read_report_limit_from_file(path, key)
        if parsed is not None:
            return parsed
    return None


def _read_report_float_from_file(
    path: Path, key: str, allow_percent: bool = False
) -> Optional[float]:
    if not path.exists():
        return None
    key_norm = key.strip().lower()
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("#") or "=" not in stripped:
                    continue
                k, raw_value = stripped.split("=", 1)
                if k.strip().lower() != key_norm:
                    continue
                if allow_percent:
                    parsed = _parse_percent_searched(raw_value)
                    if parsed is not None:
                        return parsed
                raw = _strip_quotes(raw_value).strip()
                if not raw:
                    continue
                try:
                    parsed = float(raw)
                except ValueError:
                    continue
                if math.isfinite(parsed):
                    return parsed
    except Exception:
        return None
    return None


def _resolve_report_float(
    repo_root: Path, run_id: str, key: str, allow_percent: bool = False
) -> Optional[float]:
    candidates = [
        repo_root / "data" / run_id / "config.txt",
        repo_root / "data" / run_id / "config_snapshot.txt",
        repo_root / "manifests" / run_id / "config.txt",
        repo_root / "config.txt",
    ]
    for path in candidates:
        parsed = _read_report_float_from_file(path, key, allow_percent=allow_percent)
        if parsed is not None:
            return parsed
    return None


def _read_report_bool_from_file(path: Path, key: str) -> Optional[bool]:
    if not path.exists():
        return None
    key_norm = key.strip().lower()
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("#") or "=" not in stripped:
                    continue
                k, raw_value = stripped.split("=", 1)
                if k.strip().lower() != key_norm:
                    continue
                parsed = _parse_filter_invalid(raw_value)
                if parsed is not None:
                    return parsed
    except Exception:
        return None
    return None


def _resolve_report_bool(repo_root: Path, run_id: str, key: str) -> Optional[bool]:
    candidates = [
        repo_root / "data" / run_id / "config.txt",
        repo_root / "data" / run_id / "config_snapshot.txt",
        repo_root / "manifests" / run_id / "config.txt",
        repo_root / "config.txt",
    ]
    for path in candidates:
        parsed = _read_report_bool_from_file(path, key)
        if parsed is not None:
            return parsed
    return None


def _parse_test_mode_value(raw: object) -> List[str]:
    if isinstance(raw, bool):
        return ["dud", "fda"] if raw else ["fda"]
    s = str(raw).strip()
    if not s:
        return ["fda"]
    lowered = s.lower()
    if lowered in _TEST_MODE_OFF_VALUES:
        return ["fda"]
    if lowered in _TEST_MODE_ON_VALUES:
        return ["dud", "fda"]
    if lowered == "default":
        return ["fda"]
    if lowered in _TEST_MODE_BOTH_VALUES:
        return ["dud", "fda"]

    tokens = [tok for tok in re.split(r"[+,\s]+", lowered) if tok]
    normalized: List[str] = []
    for tok in tokens:
        if tok in {"and", "off", "none", "null"}:
            continue
        if tok == "default":
            tok = "fda"
        normalized.append(tok)
    if not normalized:
        return ["fda"]
    deduped: List[str] = []
    seen: set[str] = set()
    for tok in normalized:
        if tok not in seen:
            seen.add(tok)
            deduped.append(tok)
    return deduped


def _resolve_test_mode_tokens(repo_root: Path, run_id: str) -> List[str]:
    env_raw = os.environ.get("TEST_MODE_ENABLE")
    if env_raw is not None:
        return _parse_test_mode_value(env_raw)
    candidates = [
        repo_root / "data" / run_id / "config.txt",
        repo_root / "data" / run_id / "config_snapshot.txt",
        repo_root / "manifests" / run_id / "config.txt",
        repo_root / "config.txt",
    ]
    for path in candidates:
        raw = _read_test_mode_from_file(path)
        if raw is not None:
            return _parse_test_mode_value(raw)
    return ["fda"]


def _infer_decoy_prefix_from_tokens(
    decoy_prefix: str, tokens: List[str]
) -> str:
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
    data_dir = repo_root / "data" / run_id
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
    except Exception:
        return str(path)


def _is_decoy_filename(name: str, decoy_prefix: str) -> bool:
    token = str(name or "").strip().lower()
    if not token:
        return False
    if token.startswith("decoy") or token.startswith("dud") or "_dud_" in token:
        return True
    decoy = str(decoy_prefix or "").strip().lower()
    if decoy and (
        token.startswith(f"{decoy}_") or f"_{decoy}_" in token or token.startswith(decoy)
    ):
        return True
    return False


def _row_is_decoy_for_split(row: Dict[str, Any], decoy_prefix: str) -> bool:
    if _as_bool(row.get("is_decoy")):
        return True
    library = _clean_text(row.get("library")).lower()
    if library in {"decoy", "dud"}:
        return True
    run_mode = _clean_text(row.get("run_mode")).lower()
    if run_mode in {"dud", "decoy"}:
        return True
    ligand_file = _clean_text(
        row.get("ligand_file") or row.get("ligand") or row.get("Ligand_ID")
    )
    return _is_decoy_filename(ligand_file, decoy_prefix)


def _normalize_library_slug(value: Any) -> str:
    text = _clean_text(value).lower()
    if not text:
        return ""
    text = re.sub(r"[^a-z0-9]+", "_", text).strip("_")
    if not text:
        return ""
    if text in {"na", "none", "null", "unknown"}:
        return ""
    return text


def _infer_library_from_source_csv(source_csv: Any, decoy_prefix: str) -> str:
    source = _clean_text(source_csv)
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
    data_dir = repo_root / "data" / run_id
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


def _clean_text(value: Any) -> str:
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return ""
    return text


def _norm_text(value: Any) -> str:
    text = _clean_text(value).lower()
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
    value = _clean_text(text)
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

    query_alias = _clean_text(resolve_ligand_display_name(query, "", fda_index))
    query_alias_norm = _norm_text(query_alias)
    if not query_alias_norm:
        return False

    for text in (ligand_display, ligand_base, ligand_raw):
        alias = _clean_text(resolve_ligand_display_name(text, "", fda_index))
        if alias and _norm_text(alias) == query_alias_norm:
            return True

    return False


def _resolve_ligand_display(
    row: Dict[str, Any], has_ligand_display: bool, fda_index: Optional[Any]
) -> str:
    display, _meta = _resolve_ligand_display_internal(row, has_ligand_display, fda_index)
    return display


def _resolve_ligand_display_internal(
    row: Dict[str, Any], has_ligand_display: bool, fda_index: Optional[Any]
) -> Tuple[str, Dict[str, Any]]:
    base = _clean_text(row.get("ligand_base"))
    lig_file = _clean_text(row.get("ligand_file") or row.get("ligand") or "")
    display = _clean_text(row.get("ligand_display")) if has_ligand_display else ""
    is_control = _as_bool(row.get("is_control"))
    mapping_used = False

    def build_meta(final_display: str) -> Dict[str, Any]:
        original_unfriendly = _is_unfriendly_display(display)
        final_unfriendly = _is_unfriendly_display(final_display)
        mapping_override = bool(
            has_ligand_display
            and not is_control
            and mapping_used
            and original_unfriendly
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
        fallback = _clean_text(row.get("ligand"))
        final_display = fallback or base
        return final_display, build_meta(final_display)

    if display:
        if fda_index is not None and _is_unfriendly_display(display):
            alt = _clean_text(resolve_ligand_display_name(base, lig_file, fda_index))
            if alt and alt != display and not _is_unfriendly_display(alt):
                mapping_used = True
                return alt, build_meta(alt)
        if not _is_unfriendly_display(display):
            return display, build_meta(display)
        if base and not _is_unfriendly_display(base):
            return base, build_meta(base)
        return display, build_meta(display)

    if fda_index is not None:
        alt = _clean_text(resolve_ligand_display_name(base, lig_file, fda_index))
        if alt and not _is_unfriendly_display(alt):
            mapping_used = True
            return alt, build_meta(alt)

    fallback = _clean_text(row.get("ligand"))
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
                lig_base = _clean_text(row.get("ligand_base"))
                lig_file = _clean_text(row.get("ligand_file") or row.get("ligand") or "")
                lig_file = os.path.basename(lig_file)
                library = _clean_text(row.get("library"))
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


def _read_decoy_prefix_from_file(path: Path) -> Optional[str]:
    if not path.exists():
        return None
    decoy_value = None
    dud_value = None
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("#") or "=" not in stripped:
                    continue
                key, raw_value = stripped.split("=", 1)
                key = key.strip()
                value = _strip_quotes(raw_value)
                if not value:
                    continue
                if key == DECOY_PREFIX_KEY:
                    decoy_value = value
                elif key == DUD_PREFIX_KEY:
                    dud_value = value
    except Exception:
        return None
    return decoy_value or dud_value


def _resolve_decoy_prefix(
    repo_root: Path, run_id: str, cli_value: Optional[str] = None
) -> str:
    if cli_value is not None:
        candidate = _strip_quotes(str(cli_value))
        if candidate:
            return candidate

    candidates = [
        repo_root / "data" / run_id / "config.txt",
        repo_root / "data" / run_id / "config_snapshot.txt",
        repo_root / "manifests" / run_id / "config.txt",
        repo_root / "config.txt",
    ]
    for path in candidates:
        value = _read_decoy_prefix_from_file(path)
        if value:
            return value
    return DECOY_PREFIX_DEFAULT


def _as_float(x: Any) -> Optional[float]:
    try:
        if x is None:
            return None
        s = str(x).strip()
        if not s:
            return None
        v = float(s)
        if math.isnan(v):
            return None
        return v
    except Exception:
        return None


def _as_int(x: Any) -> int:
    try:
        return int(float(x))
    except Exception:
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
        except Exception:
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
    combo_dir = repo_root / "post_docked" / run_id / pdb_id / variant / ph

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
                        mu = _as_float(row.get("scorch_mu_decoy"))
                        sigma = _as_float(row.get("scorch_sigma_decoy"))
                        n = _as_int(row.get("scorch_n_decoys"))
                        if mu is not None or sigma is not None or n > 0:
                            stats["mu"] = mu
                            stats["sigma"] = sigma
                            stats["n"] = n if n > 0 else None
                            return stats
            except Exception:
                continue
    return stats


def _as_bool(x: Any) -> bool:
    if isinstance(x, bool):
        return x
    if isinstance(x, (int, float)):
        return bool(x)
    s = str(x).strip().lower()
    return s in ("1", "true", "yes", "on")


def _resolve_target_name(pdb_id: str, cfg: Dict[str, Any]) -> str:
    normalized_pdb = _clean_text(pdb_id).upper()
    if not normalized_pdb:
        return ""
    cached = _TARGET_NAME_CACHE.get(normalized_pdb)
    if cached is not None:
        return cached

    prefer_raw = _clean_text(
        cfg.get("target_name_prefer")
        or cfg.get("TARGET_NAME_PREFER")
        or cfg.get("target_name_preference")
        or cfg.get("TARGET_NAME_PREFERENCE")
        or "auto"
    ).lower()
    prefer = prefer_raw if prefer_raw in {"auto", "compnd", "uniprot"} else "auto"
    repo_root = Path(_clean_text(cfg.get("repo_root"))).resolve() if _clean_text(cfg.get("repo_root")) else None
    try:
        if repo_root and repo_root.exists():
            cwd = Path.cwd()
            os.chdir(repo_root)
            try:
                target_name = _clean_text(
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
            target_name = _clean_text(
                derive_target_name(
                    normalized_pdb,
                    prefer=prefer,
                    pdb_root_override=Path("input_pdbs"),
                    cfg=cfg,
                )
            )
    except Exception:
        target_name = ""

    _TARGET_NAME_CACHE[normalized_pdb] = target_name
    return target_name


def _ligand_identity_for_library(row: Dict[str, Any]) -> str:
    lig_base = _clean_text(row.get("ligand_base"))
    if lig_base:
        return lig_base
    for key in ("ligand", "ligand_file"):
        raw = _clean_text(row.get(key))
        if not raw:
            continue
        basename = os.path.basename(raw)
        stem = Path(basename).stem if basename else ""
        token = _clean_text(stem or basename)
        if token:
            return token
    return ""


def _total_library_n(group_rows: List[Dict[str, Any]]) -> int:
    unique_ligands: Set[str] = set()
    for row in group_rows:
        if _as_bool(row.get("is_decoy")):
            continue
        ligand_id = _ligand_identity_for_library(row)
        if ligand_id:
            unique_ligands.add(ligand_id)
    return max(1, len(unique_ligands))


def _target_display_name(pdb_id: str, target_name: str) -> str:
    pdb = _clean_text(pdb_id)
    name = _clean_text(target_name)
    if name and pdb:
        return f"{name} ({pdb})"
    if pdb:
        return pdb
    return name


def _format_pct_display(
    value: Any, decimals: int = 1, min_nonzero_pct: Optional[float] = None
) -> str:
    pct = _as_float(value)
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


def _median(values: List[float]) -> Optional[float]:
    if not values:
        return None
    sorted_values = sorted(float(v) for v in values)
    mid = len(sorted_values) // 2
    if len(sorted_values) % 2 == 1:
        return sorted_values[mid]
    return (sorted_values[mid - 1] + sorted_values[mid]) / 2.0


def _compute_fdr_stats(group_rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    non_decoys = [r for r in group_rows if not _as_bool(r.get("is_decoy"))]
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
        dec = _as_float(r.get("fdr_n_decoys"))
        if dec is not None:
            n_decoys = int(dec)
            break
    for r in non_decoys:
        uniq = _as_float(r.get("fdr_unique_decoy_scores"))
        if uniq is not None:
            unique_decoys = int(uniq)
            break

    n_tested = 0
    n_hits_q05 = 0
    n_hits_q10 = 0
    best_q = None

    for r in non_decoys:
        q_val = _as_float(r.get("fdr_q_target"))
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
    vals = [
        _as_float(row.get(f"{prefix}_x")),
        _as_float(row.get(f"{prefix}_y")),
        _as_float(row.get(f"{prefix}_z")),
    ]
    if all(v is not None and math.isfinite(v) for v in vals):
        return [float(v) for v in vals]  # type: ignore
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
        center_vals = [_as_float(v) for v in center]
        if all(v is not None for v in center_vals):
            normalized_center = [float(v) for v in center_vals]  # type: ignore

    if isinstance(box, (list, tuple)) and len(box) == 3:
        box_vals = [_as_float(v) for v in box]
        if all(v is not None for v in box_vals):
            normalized_box = [float(v) for v in box_vals]  # type: ignore

    return {
        "method": pocket.get("method"),
        "center": normalized_center,
        "box": normalized_box,
    }


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
    tokens = test_mode_tokens or _resolve_test_mode_tokens(repo_root, run_id)
    if master_csvs_override:
        master_csvs = [Path(p) for p in master_csvs_override]
    else:
        master_csvs = _resolve_master_csvs(repo_root, run_id, tokens, decoy_prefix)
    if not master_csvs:
        raise FileNotFoundError(
            f"Master CSV not found for run_id={run_id} tokens={tokens}"
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
    multi_target_max_pct = float(multi_target_max_pct) if multi_target_max_pct is not None else 0.0
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
        logger.warning("%s action=fda_mapping status=missing run_id=%s", COMPONENT, run_id)
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
        key = f"{pdb}|{variant}|{ph}"
        target_groups.setdefault(key, []).append(r)

        lig_base = _clean_text(r.get("ligand_base"))
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
        pdb, variant, ph = key.split("|")
        target_name = _resolve_target_name(pdb, target_name_cfg)
        target_display = _target_display_name(pdb, target_name)

        # QC Stats
        n_rows = len(group_rows)
        n_decoys = sum(1 for r in group_rows if _as_bool(r.get("is_decoy")))
        n_controls = sum(1 for r in group_rows if _as_bool(r.get("is_control")))

        # Extract metrics from first row (assuming they are repeated/consistent per target)
        first = group_rows[0]
        ef1 = _as_float(first.get("ef1"))
        roc_auc = _as_float(first.get("roc_auc"))
        roc_auc_adj = _as_float(first.get("roc_auc_adj")) or roc_auc
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

        # Rank non-decoys with valid t_selected
        rankable: List[tuple[Dict[str, Any], float]] = []
        for r in group_rows:
            if _as_bool(r.get("is_decoy")):
                continue
            if filter_invalid and not _as_bool(r.get("pose_valid_any")):
                continue
            t_val = _as_float(r.get("t_selected"))
            if t_val is None or not math.isfinite(t_val):
                continue
            rankable.append((r, t_val))

        def rank_key(item: tuple[Dict[str, Any], float]) -> tuple:
            row, t_val = item
            lig_base = _clean_text(row.get("ligand_base"))
            lig_file = _clean_text(row.get("ligand_file") or row.get("ligand") or "")
            return (-t_val, lig_base, lig_file)

        sorted_rankable = sorted(rankable, key=rank_key)
        ranked_entries = []
        ranked_rows: List[Dict[str, Any]] = []
        total_library_n = _total_library_n(group_rows)
        for idx, (r, t_sel) in enumerate(sorted_rankable, start=1):
            lig_base = _clean_text(r.get("ligand_base"))
            lig_display = _resolve_ligand_display(r, has_ligand_display, fda_index)
            t_src = r.get("t_selected_source") or None
            pct = idx / total_library_n

            entry = {
                "ligand_base": lig_base,
                "ligand_display": lig_display,
                "t_selected": t_sel,
                "t_selected_source": t_src,
                "pose_valid_any": _as_bool(r.get("pose_valid_any")),
                "pose_invalid_reason_top": _clean_text(
                    r.get("pose_invalid_reason_top")
                ),
                "rank": idx,
                "pct_rank": float(f"{pct:.6f}"),
                "is_control": _as_bool(r.get("is_control")),
                "is_decoy": _as_bool(r.get("is_decoy")),
            }
            ranked_entries.append(entry)
            ranked_rows.append({"row": r, "entry": entry})
            if lig_base and lig_base not in ligand_display_map:
                ligand_display_map[lig_base] = lig_display or lig_base

        compact_best_by_ligand: Dict[str, Dict[str, Any]] = {}
        for r in group_rows:
            if _as_bool(r.get("is_decoy")):
                continue
            t_val = _as_float(r.get("t_selected"))
            if t_val is None or not math.isfinite(t_val):
                continue
            if coverage_requires_pose_valid and not _as_bool(r.get("pose_valid_any")):
                continue
            lig_base = _ligand_identity_for_library(r)
            if not lig_base:
                continue
            lig_display = _resolve_ligand_display(r, has_ligand_display, fda_index)
            lig_file = _clean_text(r.get("ligand_file") or r.get("ligand") or "")
            candidate = {
                "ligand_base": lig_base,
                "ligand_display": lig_display or lig_base,
                "ligand_file": lig_file,
                "t_selected": t_val,
                "t_selected_source": r.get("t_selected_source") or None,
            }
            current = compact_best_by_ligand.get(lig_base)
            if current is None:
                compact_best_by_ligand[lig_base] = candidate
                continue
            current_t = _as_float(current.get("t_selected"))
            if current_t is None or t_val > current_t:
                compact_best_by_ligand[lig_base] = candidate
                continue
            if current_t is not None and t_val == current_t and lig_file < str(
                current.get("ligand_file") or ""
            ):
                compact_best_by_ligand[lig_base] = candidate

        def compact_rank_key(item: Dict[str, Any]) -> tuple:
            t_val = _as_float(item.get("t_selected"))
            return (
                -(t_val if t_val is not None else float("-inf")),
                _clean_text(item.get("ligand_base")),
                _clean_text(item.get("ligand_file")),
            )

        compact_rankable = sorted(compact_best_by_ligand.values(), key=compact_rank_key)
        for idx, item in enumerate(compact_rankable, start=1):
            lig_base = _clean_text(item.get("ligand_base"))
            if not lig_base:
                continue
            lig_display = _clean_text(item.get("ligand_display")) or lig_base
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
                    "t_selected": float(_as_float(item.get("t_selected")) or 0.0),
                    "t_selected_source": item.get("t_selected_source"),
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
                    ligand_raw = _clean_text(
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
                    t_val = entry["t_selected"]
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
                    "t_selected": entry["t_selected"],
                    "t_selected_source": entry["t_selected_source"],
                    "pose_valid_any": entry["pose_valid_any"],
                    "pose_invalid_reason_top": entry["pose_invalid_reason_top"],
                    "rank": entry["rank"],
                    "pct_rank": entry["pct_rank"],
                    "is_control": entry["is_control"],
                    "is_decoy": entry["is_decoy"],
                }
                library = _clean_text(row.get("library"))
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
                "t_selected": entry["t_selected"],
                "t_selected_source": entry["t_selected_source"],
                "ligand_display": entry["ligand_display"],
                "is_control": entry["is_control"],
            }
            library = _clean_text(row.get("library"))
            if library:
                detail["library"] = library
            multi_target_entries.setdefault(lig_base, []).append(detail)

        for entry in top_list:
            lig_base = _clean_text(entry.get("ligand_base"))
            if lig_base:
                ligand_top5_hits.setdefault(lig_base, []).append(
                    {
                        "target": key,
                        "t_selected": entry["t_selected"],
                        "t_selected_source": entry["t_selected_source"],
                    }
                )

        for entry in notable_list:
            lig_base = _clean_text(entry.get("ligand_base"))
            if lig_base:
                ligand_notable_hits.setdefault(lig_base, []).append(
                    {
                        "target": key,
                        "t_selected": entry["t_selected"],
                        "t_selected_source": entry["t_selected_source"],
                        "rank": entry["rank"],
                        "pct_rank": entry["pct_rank"],
                    }
                )

        fdr_stats = _compute_fdr_stats(group_rows)
        fdr_value: Optional[Any] = None
        if fdr_stats:
            best_q = _as_float(fdr_stats.get("best_q"))
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
    all_bases = sorted(
        set(ligand_top5_hits.keys()) | set(ligand_notable_hits.keys())
    )
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
            t_values = [float(e["t_selected"]) for e in qualifying]
            worst_pct = max(pct_values)
            mean_pct = sum(pct_values) / len(pct_values)
            best_rank = min(rank_values)
            best_t_selected = max(t_values)
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
                    "t_selected": entry["t_selected"],
                    "t_selected_source": entry["t_selected_source"],
                }
                library = _clean_text(entry.get("library"))
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
                    "best_t_selected": best_t_selected,
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
                float(_as_float(item.get("pct_rank")) or 1.0),
                -(float(_as_float(item.get("t_selected")) or float("-inf"))),
                _clean_text(item.get("target_id")),
            ),
        )
        best = entries_sorted[0]
        best_score = float(_as_float(best.get("t_selected")) or 0.0)
        best_pct_rank = float(_as_float(best.get("pct_rank")) or 1.0)
        best_target = _clean_text(best.get("target_display")) or _clean_text(
            best.get("pdb_id")
        )

        top_entries = entries_sorted[:top_targets_n]
        top_lines = []
        for idx, entry in enumerate(top_entries, start=1):
            target_display = _clean_text(entry.get("pdb_id"))
            z_text = _format_num(entry.get("t_selected"), ".6g")
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
            tail_median = _median(tail_pcts)
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
            if float(_as_float(entry.get("pct_rank")) or 1.0) <= breadth_pct_threshold
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
            float(_as_float(item.get("best_pct_rank")) or 1.0),
            -(float(_as_float(item.get("best_score")) or 0.0)),
            _clean_text(item.get("ligand_base")),
        )
    )
    highlights_top_count = len(ligand_highlights)
    if highlights_top_count and highlights_top_pct > 0:
        pct_fraction = min(float(highlights_top_pct), 1.0)
        highlights_top_count = max(1, int(math.ceil(highlights_top_count * pct_fraction)))
    filtered_highlights = ligand_highlights[:highlights_top_count]
    if highlights_max_ligands > 0:
        ligand_highlights = filtered_highlights[:highlights_max_ligands]
    else:
        ligand_highlights = filtered_highlights

    sources: Dict[str, Any] = {
        "master_rows_csv": _safe_relpath(primary_master_csv, repo_root),
        "manifest_yaml": str(canonical_manifest_rel),
    }
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


def write_yaml(report: Dict[str, Any], out_path: Path) -> None:
    # Use a compact representation for lists of numbers (center/box) if possible?
    # PyYAML default dump is okay.
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        yaml.dump(report, f, sort_keys=False, default_flow_style=False)


def _format_num(value: Any, fmt: str) -> str:
    try:
        num = float(value)
        if math.isfinite(num):
            return format(num, fmt)
    except Exception:
        return ""
    return ""


def _stage_artifacts(
    output_dir: Path, artifacts: List[Tuple[str, Path]]
) -> List[Tuple[str, str]]:
    staged: List[Tuple[str, str]] = []
    if not artifacts:
        return staged
    artifacts_dir = output_dir / "artifacts"
    logger = logging.getLogger("run-report")
    for label, src in artifacts:
        if not src.exists():
            continue
        dest = artifacts_dir / src.name
        if dest.exists():
            try:
                if dest.resolve() == src.resolve():
                    rel_path = dest.relative_to(output_dir).as_posix()
                    staged.append((label, rel_path))
                    continue
            except Exception:
                pass
            stem = src.stem
            suffix = src.suffix
            counter = 1
            while dest.exists():
                dest = artifacts_dir / f"{stem}-{counter}{suffix}"
                counter += 1
        try:
            artifacts_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dest)
        except Exception as exc:
            logger.warning(
                "%s action=stage_artifact_failed src=%s error=%s",
                COMPONENT,
                src,
                exc,
            )
            continue
        rel_path = dest.relative_to(output_dir).as_posix()
        staged.append((label, rel_path))
    return staged


def _resolve_heatmap_source(repo_root: Path, run_id: str) -> Optional[Path]:
    interactions_dataset = repo_root / "data" / run_id / "dataset" / "interactions"
    heatmap_csv = repo_root / "data" / run_id / "heatmap_input.csv"
    master_csv = repo_root / "data" / run_id / "master_rows.csv"
    # Preserve existing report source preference.
    if heatmap_csv.exists():
        return heatmap_csv
    if interactions_dataset.exists():
        return interactions_dataset
    if master_csv.exists():
        return master_csv
    return None


def _extract_pdb_id_from_target_id(target_id: Any) -> str:
    token = str(target_id or "").strip().split("|", 1)[0].strip().upper()
    return token


def _collect_pdb_ids_from_manifest(repo_root: Path, run_id: str) -> Set[str]:
    manifest, _manifest_path = load_run_manifest(repo_root, run_id)
    proteins = (manifest or {}).get("proteins", {}) if isinstance(manifest, dict) else {}
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
    fallback_csv = repo_root / "data" / run_id / "heatmap_input.csv"
    return _collect_pdb_ids_from_heatmap_csv(fallback_csv)


def _render_pathway_heatmap_document(run_id: str, pathway_slug: str, body_html: str) -> str:
    title = f"{run_id} pathway heatmap: {pathway_slug}"
    return (
        "<!doctype html>\n"
        "<html>\n"
        "<head>\n"
        "  <meta charset=\"utf-8\">\n"
        "  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
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
) -> List[Path]:
    run_id = str(report.get("run_id") or "").strip()
    if not run_id:
        return []

    heatmap_source = heatmap_source_override or _resolve_heatmap_source(repo_root, run_id)
    if heatmap_source is None:
        return []

    pdb_ids = sorted(_collect_pdb_ids_for_pathway_reports(repo_root, run_id, heatmap_source))
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

    eligible_groups = {
        slug: pdb_set for slug, pdb_set in grouped.items() if len(pdb_set) >= 2
    }
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
            )
            html_payload = _render_pathway_heatmap_document(
                run_id, pathway_slug, heatmap_html
            )
        except Exception as exc:
            html_payload = _render_pathway_heatmap_document(
                run_id,
                pathway_slug,
                (
                    "<div>Heatmap unavailable: "
                    f"{html.escape(str(exc))}"
                    "</div>"
                ),
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
) -> List[Path]:
    _write_html_report(
        report,
        out_path,
        highlight_queries,
        repo_root,
        heatmap_source_override=heatmap_source_override,
    )
    return _write_pathway_specific_heatmap_reports(
        report,
        out_path,
        repo_root,
        heatmap_source_override=heatmap_source_override,
    )


def _write_html_report(
    report: Dict[str, Any],
    out_path: Path,
    highlight_queries: List[str],
    repo_root: Path,
    heatmap_source_override: Optional[Path] = None,
) -> None:
    run_id = report.get("run_id", "")
    generated_at = report.get("generated_at", "")

    compact_rows = []
    compact_highlights = report.get("ligand_highlights", []) or []
    compact_cfg = report.get("ligand_highlight_config", {}) or {}
    top_targets_n = int(compact_cfg.get("top_targets_n") or 3)
    breadth_threshold_raw = _as_float(compact_cfg.get("breadth_pct_threshold"))
    if breadth_threshold_raw is None or not math.isfinite(breadth_threshold_raw):
        breadth_threshold_raw = 0.01
    breadth_threshold_label = f"{(breadth_threshold_raw * 100.0):g}%"
    if compact_highlights:
        for row in compact_highlights:
            ligand_name = _clean_text(row.get("ligand")) or _clean_text(
                row.get("ligand_display")
            )
            if not ligand_name:
                ligand_name = _clean_text(row.get("ligand_base")) or "—"
            best_target = _clean_text(row.get("best_target")) or "—"
            best_score = _format_num(row.get("best_score"), ".6g") or "—"
            best_percentile = _clean_text(row.get("best_percentile"))
            if not best_percentile:
                best_percentile = (
                    _format_pct_display(
                        row.get("best_pct_rank"),
                        2,
                        min_nonzero_pct=0.01,
                    )
                    or "—"
                )
            top_targets = _clean_text(row.get("top_targets")) or "—"
            coverage = _clean_text(row.get("coverage")) or "—"
            breadth = _clean_text(row.get("breadth")) or "0"
            tail_summary = _clean_text(row.get("tail_summary")) or "—"
            compact_rows.append(
                "<tr>"
                f"<td>{html.escape(ligand_name)}</td>"
                f"<td>{html.escape(best_target)}</td>"
                f"<td class=\"num\">{html.escape(best_score)}</td>"
                f"<td class=\"num best-pct\">{html.escape(best_percentile)}</td>"
                f"<td class=\"top-targets-col top-targets\">{html.escape(top_targets)}</td>"
                f"<td class=\"num\">{html.escape(coverage)}</td>"
                f"<td class=\"num\">{html.escape(breadth)}</td>"
                f"<td class=\"num\">{html.escape(tail_summary)}</td>"
                "</tr>"
            )

    highlight_rows = []
    highlight_headers = "".join(
        f"<th class=\"num\">{html.escape(q)}</th>" for q in highlight_queries
    )
    if not compact_rows:
        targets = report.get("targets", {}) or {}
        for target_id in sorted(targets.keys()):
            target = targets[target_id]
            target_name = _clean_text(target.get("target_name")) or "—"
            highlights = target.get("highlights", []) or []
            highlight_map = {h.get("query"): h for h in highlights if h}
            cells = []
            for query in highlight_queries:
                entry = highlight_map.get(query)
                if entry and entry.get("found"):
                    rank = entry.get("rank")
                    pct = _format_num(entry.get("pct_rank"), ".6f")
                    t_sel = _format_num(entry.get("t_selected"), ".6g")
                    if rank:
                        cell = f"{rank} ({pct}) {t_sel}".strip()
                    else:
                        cell = "not found"
                else:
                    cell = "not found"
                cells.append(f"<td class=\"num\">{html.escape(cell)}</td>")
            row_html = (
                f"<tr><td>{html.escape(target_id)}</td><td>{html.escape(target_name)}</td>{''.join(cells)}</tr>"
            )
            highlight_rows.append(row_html)

    heatmap_html = ""
    heatmap_source = heatmap_source_override or _resolve_heatmap_source(repo_root, run_id)
    if heatmap_source is not None:
        try:
            heatmap_html = render_interactive_heatmap_html(
                repo_root, run_id, heatmap_source
            )
        except Exception as exc:
            heatmap_html = (
                f"<div class=\"meta\">Heatmap unavailable: {html.escape(str(exc))}</div>"
            )
    else:
        heatmap_html = "<div class=\"meta\">Heatmap data not available.</div>"

    heatmap_section = (
        "<section class=\"section\" id=\"heatmap\">"
        "<h2>Heatmap</h2>"
        f"{heatmap_html}"
        "</section>"
    )

    report_suffix = ""
    if out_path.stem.startswith("report_"):
        report_suffix = out_path.stem[len("report_") :]
    data_dir = repo_root / "data" / run_id
    if report_suffix:
        report_yaml = data_dir / f"report_{report_suffix}.yaml"
        heatmap_csv_artifact = data_dir / f"heatmap_input_{report_suffix}.csv"
        heatmap_png = data_dir / f"heatmap_{report_suffix}.png"
    else:
        report_yaml = data_dir / "report.yaml"
        heatmap_csv_artifact = data_dir / "heatmap_input.csv"
        heatmap_png = data_dir / "heatmap.png"
    artifacts = [
        ("report.yaml", report_yaml),
        ("heatmap_input.csv", heatmap_csv_artifact),
        ("heatmap.png", heatmap_png),
    ]
    staged = _stage_artifacts(out_path.parent, artifacts)
    link_items = [
        f"<li><a href=\"{html.escape(href)}\">{html.escape(label)}</a></li>"
        for label, href in staged
    ]
    if link_items:
        artifacts_html = f"<ul class=\"artifact-list\">{''.join(link_items)}</ul>"
    else:
        artifacts_html = "<div class=\"meta\">No artifacts available.</div>"

    html_body = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(str(run_id))} report</title>
  <style>
    :root {{
      color-scheme: light;
      --bg: #f6f7fb;
      --card-bg: #ffffff;
      --text: #1f2937;
      --muted: #6b7280;
      --border: #e5e7eb;
      --accent: #0f172a;
      --highlight: #eef2ff;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      font-family: "Trebuchet MS", "Lucida Grande", "Lucida Sans Unicode", sans-serif;
      margin: 0;
      background: var(--bg);
      color: var(--text);
      line-height: 1.5;
    }}
    h1, h2 {{ margin: 0; color: var(--accent); }}
    .container {{ max-width: 1200px; margin: 0 auto; padding: 28px 20px 48px; }}
    .page-header {{
      display: flex;
      flex-wrap: wrap;
      justify-content: space-between;
      align-items: baseline;
      gap: 12px;
      margin-bottom: 16px;
    }}
    .page-title {{
      font-family: "Palatino Linotype", "Book Antiqua", Palatino, serif;
      font-size: 2.1rem;
    }}
    .meta {{ color: var(--muted); font-size: 0.95rem; }}
    .toc {{
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      margin: 16px 0 20px;
      padding: 10px 12px;
      background: var(--card-bg);
      border: 1px solid var(--border);
      border-radius: 12px;
    }}
    .toc a {{
      text-decoration: none;
      color: var(--accent);
      background: #f0f3f8;
      border: 1px solid #e2e8f0;
      padding: 4px 10px;
      border-radius: 999px;
      font-size: 0.9rem;
    }}
    .toc a:hover {{ background: #e8ecf4; }}
    .section {{
      background: var(--card-bg);
      border: 1px solid var(--border);
      border-radius: 14px;
      padding: 16px 18px;
      margin-bottom: 18px;
      box-shadow: 0 1px 2px rgba(15, 23, 42, 0.04);
    }}
    .section h2 {{ margin-bottom: 12px; font-size: 1.3rem; }}
    .table-wrap {{
      overflow-x: auto;
      border: 1px solid var(--border);
      border-radius: 10px;
      background: #ffffff;
    }}
    table {{
      border-collapse: collapse;
      width: 100%;
      font-size: 0.95rem;
      font-variant-numeric: tabular-nums;
    }}
    th, td {{ padding: 8px 10px; text-align: left; }}
    thead th {{
      position: sticky;
      top: 0;
      background: #f4f6fa;
      border-bottom: 1px solid var(--border);
      z-index: 1;
    }}
    tbody td {{ border-bottom: 1px solid var(--border); }}
    tbody tr:nth-child(even) {{ background: #f9fafb; }}
    tbody tr:hover {{ background: var(--highlight); }}
    th.num, td.num {{ text-align: right; }}
    .highlights-table th, .highlights-table td {{ vertical-align: top; }}
    .highlights-table .best-pct {{
      white-space: nowrap;
      min-width: 7.5rem;
      padding-right: 14px;
    }}
    .highlights-table .top-targets-col {{
      min-width: 22rem;
      padding-left: 14px;
    }}
    .top-targets {{
      white-space: pre;
      font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace;
      font-size: 0.9rem;
      line-height: 1.4;
      text-align: left;
    }}
    .artifact-list {{
      list-style: none;
      padding-left: 0;
      margin: 0;
      display: grid;
      gap: 6px;
    }}
    .artifact-list a {{
      text-decoration: none;
      color: var(--accent);
      background: #f8fafc;
      border: 1px solid var(--border);
      padding: 6px 10px;
      border-radius: 8px;
      display: inline-block;
    }}
    .artifact-list a:hover {{ background: #eef2f7; }}
  </style>
</head>
<body>
  <div class="container">
    <header class="page-header">
      <h1 class="page-title">{html.escape(str(run_id))} report</h1>
      <div class="meta">Generated at {html.escape(str(generated_at))}</div>
    </header>

    <nav class="toc" aria-label="Table of contents">
      <a href="#highlights">Highlights</a>
      <a href="#heatmap">Heatmap</a>
      <a href="#artifacts">Artifacts</a>
    </nav>

    <section class="section" id="highlights">
      <h2>Highlights</h2>
      <div class="table-wrap">
        <table class="highlights-table">
          <thead>
            <tr>
              {'<th>Ligand</th><th>Best target</th><th class="num">Best score</th><th class="num best-pct">Best percentile</th>'
                + f'<th class="top-targets-col">Top targets (N={top_targets_n})</th><th class="num">Coverage</th><th class="num">Breadth {breadth_threshold_label}</th><th class="num">Tail summary</th>'
                if compact_rows else '<th>target_id</th><th>target_name</th>' + highlight_headers}
            </tr>
          </thead>
          <tbody>
            {''.join(compact_rows) if compact_rows else ''.join(highlight_rows)}
          </tbody>
        </table>
      </div>
    </section>

    {heatmap_section}

    <section class="section" id="artifacts">
      <h2>Artifacts</h2>
      {artifacts_html}
    </section>
  </div>
</body>
</html>
"""

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as handle:
        handle.write(html_body)


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
    master_csv = master_csv_override or (repo_root / "data" / run_id / "master_rows.csv")
    if not master_csv.exists():
        raise FileNotFoundError(f"Master CSV not found: {master_csv}")

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
        logger.warning("%s action=fda_mapping status=missing run_id=%s", COMPONENT, run_id)
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
        pdb = _clean_text(row.get("pdb_id"))
        variant = _clean_text(row.get("variant"))
        ph = _clean_text(row.get("ph_label"))
        target_id = f"{pdb}|{variant}|{ph}"
        rows_by_target.setdefault(target_id, []).append(row)

    for target_id, group_rows in rows_by_target.items():
        rankable: List[tuple[Dict[str, Any], float]] = []
        for row in group_rows:
            if _as_bool(row.get("is_decoy")):
                continue
            if filter_invalid and not _as_bool(row.get("pose_valid_any")):
                continue
            t_val = _as_float(row.get("t_selected"))
            if t_val is None or not math.isfinite(t_val):
                continue
            rankable.append((row, t_val))

        def rank_key(item: tuple[Dict[str, Any], float]) -> tuple:
            row, t_val = item
            lig_base = _clean_text(row.get("ligand_base"))
            lig_file = _clean_text(row.get("ligand_file") or row.get("ligand") or "")
            return (-t_val, lig_base, lig_file)

        sorted_rankable = sorted(rankable, key=rank_key)
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
                if filter_invalid and not _as_bool(row.get("pose_valid_any")):
                    continue
                rank_info = row_rank.get(id(row))
                is_control = _as_bool(row.get("is_control"))
                if is_control:
                    filtered_rows.append(row)
                elif rank_info and rank_info["rank"] <= top_k:
                    filtered_rows.append(row)
        rows = filtered_rows

    if filter_invalid and (top_k is None or top_k <= 0):
        rows = [row for row in rows if _as_bool(row.get("pose_valid_any"))]

    def sort_key(row: Dict[str, Any]) -> tuple:
        rank_info = row_rank.get(id(row)) or {}
        rank_val = rank_info.get("rank")
        sort_rank = rank_val if isinstance(rank_val, int) else 1_000_000
        target_id = rank_info.get("target_id")
        if not target_id:
            pdb = _clean_text(row.get("pdb_id"))
            variant = _clean_text(row.get("variant"))
            ph = _clean_text(row.get("ph_label"))
            target_id = f"{pdb}|{variant}|{ph}"
        lig_base = _clean_text(row.get("ligand_base"))
        lig_file = _clean_text(row.get("ligand_file") or row.get("ligand") or "")
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
            pdb = _clean_text(row.get("pdb_id"))
            variant = _clean_text(row.get("variant"))
            ph = _clean_text(row.get("ph_label"))
            target_id = f"{pdb}|{variant}|{ph}"
            display = _resolve_ligand_display(row, has_ligand_display, fda_index)
            writer.writerow(
                {
                    "target_id": target_id,
                    "target_name": _resolve_target_name(pdb, target_name_cfg),
                    "pdb_id": pdb,
                    "variant": variant,
                    "ph_label": ph,
                    "ligand_base": _clean_text(row.get("ligand_base")),
                    "ligand_display": display,
                    "library": _clean_text(row.get("library")),
                    "z_selected": _clean_text(row.get("t_selected")),
                    "pose_valid_any": _clean_text(row.get("pose_valid_any")),
                    "pose_invalid_reason_top": _clean_text(
                        row.get("pose_invalid_reason_top")
                    ),
                    "rank": rank_info.get("rank", ""),
                    "pct_rank": rank_info.get("pct_rank", ""),
                    "is_decoy": _clean_text(row.get("is_decoy")),
                    "is_control": _clean_text(row.get("is_control")),
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate run report YAML")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--overwrite", action="store_true", default=True)
    parser.add_argument("--no-overwrite", action="store_false", dest="overwrite")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--decoy-prefix", default=None)
    parser.add_argument("--top-n", type=int, default=5)
    parser.add_argument("--extended-top-n", type=int, default=5)
    parser.add_argument("--notable-pct", type=float, default=0.01)
    parser.add_argument("--notable-max", type=int, default=25)
    parser.add_argument("--fda-mapping-csv", default=None)
    parser.add_argument("--highlight-ligands", default="imatinib")
    parser.add_argument(
        "--highlight-match",
        choices=["display_contains", "display_exact", "any_contains"],
        default="any_contains",
    )
    parser.add_argument("--highlight-max-per-target", type=int, default=5)
    parser.add_argument("--multi-target-min-targets", type=int, default=2)
    parser.add_argument("--multi-target-max-pct", type=float, default=0.05)
    parser.add_argument("--multi-target-max-hits", type=int, default=50)
    parser.add_argument(
        "--multi-target-sort",
        choices=["worst_pct_then_mean", "mean_pct", "best_rank_sum"],
        default="worst_pct_then_mean",
    )
    parser.add_argument("--emit-html", action="store_true")
    parser.add_argument("--html-path", default=None)
    parser.add_argument(
        "--emit-heatmap-csv",
        action="store_true",
        default=True,
        dest="emit_heatmap_csv",
    )
    parser.add_argument(
        "--no-emit-heatmap-csv",
        action="store_false",
        dest="emit_heatmap_csv",
    )
    parser.add_argument("--heatmap-csv-path", default=None)
    parser.add_argument("--heatmap-top-k", type=int, default=200)
    parser.add_argument("--top-targets-n", type=int, default=3)
    parser.add_argument("--breadth-pct-threshold", type=float, default=0.01)
    parser.add_argument(
        "--coverage-requires-pose-valid",
        action="store_true",
        default=False,
    )
    parser.add_argument("--pct-display-decimals", type=int, default=1)
    parser.add_argument("--tail-median-decimals", type=int, default=0)
    parser.add_argument("--highlights-top-pct", type=float, default=0.01)
    parser.add_argument("--highlights-max-ligands", type=int, default=25)
    parser.add_argument(
        "-combined",
        "--combined",
        action="store_true",
        default=False,
        help="Emit a single combined report instead of per-library reports.",
    )
    args = parser.parse_args()

    logger = _configure_logging(args.verbose)
    repo_root = Path(args.repo_root).resolve()
    run_id = args.run_id
    decoy_prefix = _resolve_decoy_prefix(repo_root, run_id, args.decoy_prefix)
    tokens = _resolve_test_mode_tokens(repo_root, run_id)
    decoy_prefix = _infer_decoy_prefix_from_tokens(decoy_prefix, tokens)
    logger.info(
        "%s action=preflight decoy_prefix=%s tokens=%s",
        COMPONENT,
        decoy_prefix,
        ",".join(tokens),
    )
    percent_searched = _resolve_percent_searched(repo_root)
    notable_pct = args.notable_pct
    multi_target_max_pct = args.multi_target_max_pct
    if percent_searched is not None:
        if "--notable-pct" not in sys.argv:
            notable_pct = percent_searched
        if "--multi-target-max-pct" not in sys.argv:
            multi_target_max_pct = percent_searched
        logger.info(
            "%s action=percent_searched value=%.6f",
            COMPONENT,
            percent_searched,
        )

    filter_invalid = _resolve_filter_invalid(repo_root, run_id)
    if filter_invalid is None:
        filter_invalid = False

    top_n = args.top_n
    extended_top_n = args.extended_top_n
    if "--top-n" not in sys.argv:
        cfg_top_n = _resolve_report_limit(repo_root, run_id, "REPORT_TOP_N")
        if cfg_top_n is not None:
            top_n = cfg_top_n
    if "--extended-top-n" not in sys.argv:
        cfg_extended = _resolve_report_limit(
            repo_root, run_id, "REPORT_EXTENDED_TOP_N"
        )
        if cfg_extended is not None:
            extended_top_n = cfg_extended

    top_targets_n = args.top_targets_n
    if "--top-targets-n" not in sys.argv:
        cfg_top_targets_n = _resolve_report_limit(repo_root, run_id, "TOP_TARGETS_N")
        if cfg_top_targets_n is not None:
            top_targets_n = cfg_top_targets_n

    breadth_pct_threshold = args.breadth_pct_threshold
    if "--breadth-pct-threshold" not in sys.argv:
        cfg_breadth_pct = _resolve_report_float(
            repo_root, run_id, "BREADTH_PCT_THRESHOLD", allow_percent=True
        )
        if cfg_breadth_pct is not None:
            breadth_pct_threshold = cfg_breadth_pct

    coverage_requires_pose_valid = args.coverage_requires_pose_valid
    if "--coverage-requires-pose-valid" not in sys.argv:
        cfg_coverage_requires = _resolve_report_bool(
            repo_root, run_id, "COVERAGE_REQUIRES_POSE_VALID"
        )
        if cfg_coverage_requires is not None:
            coverage_requires_pose_valid = cfg_coverage_requires

    pct_display_decimals = args.pct_display_decimals
    if "--pct-display-decimals" not in sys.argv:
        cfg_pct_display_decimals = _resolve_report_limit(
            repo_root, run_id, "PCT_DISPLAY_DECIMALS"
        )
        if cfg_pct_display_decimals is not None:
            pct_display_decimals = cfg_pct_display_decimals

    tail_median_decimals = args.tail_median_decimals
    if "--tail-median-decimals" not in sys.argv:
        cfg_tail_median_decimals = _resolve_report_limit(
            repo_root, run_id, "TAIL_MEDIAN_DECIMALS"
        )
        if cfg_tail_median_decimals is not None:
            tail_median_decimals = cfg_tail_median_decimals

    highlights_top_pct = args.highlights_top_pct
    highlights_max_ligands = args.highlights_max_ligands

    master_csvs = _resolve_master_csvs(repo_root, run_id, tokens, decoy_prefix)
    if not master_csvs:
        logger.error(
            "%s action=fail error=%s",
            COMPONENT,
            f"Master CSV not found for run_id={run_id} tokens={tokens}",
        )
        return 1

    jobs: List[Tuple[Optional[str], List[Path]]] = []
    if args.combined:
        jobs.append((None, master_csvs))
    else:
        split_map = _write_split_master_rows(repo_root, run_id, decoy_prefix, master_csvs)
        split_keys = sorted(split_map.keys())
        if len(split_keys) > 1:
            for lib_slug in split_keys:
                jobs.append((lib_slug, [split_map[lib_slug]]))
        elif len(split_keys) == 1:
            only = split_keys[0]
            jobs.append((None, [split_map[only]]))
        else:
            jobs.append((None, master_csvs))

    try:
        highlight_queries = _parse_highlight_queries(args.highlight_ligands)
        if args.highlight_max_per_target > 0:
            highlight_queries = highlight_queries[: args.highlight_max_per_target]
        else:
            highlight_queries = []
        top_k = args.heatmap_top_k if args.heatmap_top_k > 0 else None

        for lib_suffix, job_master_csvs in jobs:
            suffix = _normalize_library_slug(lib_suffix) if lib_suffix else ""
            data_dir = repo_root / "data" / run_id
            if suffix:
                yaml_path = data_dir / f"report_{suffix}.yaml"
                default_heatmap_path = data_dir / f"heatmap_input_{suffix}.csv"
                default_html_path = data_dir / f"report_{suffix}.html"
            else:
                yaml_path = data_dir / "report.yaml"
                default_heatmap_path = data_dir / "heatmap_input.csv"
                default_html_path = data_dir / "report.html"

            if yaml_path.exists() and not args.overwrite:
                logger.info("%s action=skip reason=exists path=%s", COMPONENT, yaml_path)
                continue

            report = build_report(
                run_id,
                repo_root,
                top_n=top_n,
                extended_top_n=extended_top_n,
                notable_pct=notable_pct,
                notable_max=args.notable_max,
                decoy_prefix=decoy_prefix,
                fda_mapping_csv=args.fda_mapping_csv,
                highlight_ligands=args.highlight_ligands,
                highlight_match=args.highlight_match,
                highlight_max_per_target=args.highlight_max_per_target,
                multi_target_min_targets=args.multi_target_min_targets,
                multi_target_max_pct=multi_target_max_pct,
                multi_target_max_hits=args.multi_target_max_hits,
                multi_target_sort=args.multi_target_sort,
                filter_invalid=filter_invalid,
                top_targets_n=top_targets_n,
                breadth_pct_threshold=breadth_pct_threshold,
                coverage_requires_pose_valid=coverage_requires_pose_valid,
                pct_display_decimals=pct_display_decimals,
                tail_median_decimals=tail_median_decimals,
                highlights_top_pct=highlights_top_pct,
                highlights_max_ligands=highlights_max_ligands,
                test_mode_tokens=tokens,
                master_csvs_override=job_master_csvs,
            )
            write_yaml(report, yaml_path)
            logger.info(
                "%s action=write status=ok path=%s library=%s",
                COMPONENT,
                yaml_path,
                suffix or "combined",
            )

            heatmap_path: Optional[Path] = None
            if args.emit_heatmap_csv:
                if args.heatmap_csv_path and len(jobs) == 1:
                    heatmap_path = Path(args.heatmap_csv_path)
                else:
                    heatmap_path = default_heatmap_path
                try:
                    _write_heatmap_input_csv(
                        repo_root,
                        run_id,
                        heatmap_path,
                        top_k,
                        args.fda_mapping_csv,
                        filter_invalid,
                        master_csv_override=job_master_csvs[0] if job_master_csvs else None,
                    )
                    logger.info(
                        "%s action=write_heatmap_csv status=ok path=%s library=%s",
                        COMPONENT,
                        heatmap_path,
                        suffix or "combined",
                    )
                except Exception as exc:
                    logger.warning(
                        "%s action=write_heatmap_csv status=failed path=%s error=%s library=%s",
                        COMPONENT,
                        heatmap_path,
                        exc,
                        suffix or "combined",
                    )
                    heatmap_path = None

            if args.emit_html:
                if args.html_path and len(jobs) == 1:
                    html_path = Path(args.html_path)
                else:
                    html_path = default_html_path
                try:
                    heatmap_source_for_html = heatmap_path
                    if heatmap_source_for_html is None and job_master_csvs:
                        heatmap_source_for_html = job_master_csvs[0]
                    pathway_htmls = _write_html_reports(
                        report,
                        html_path,
                        highlight_queries,
                        repo_root,
                        heatmap_source_override=heatmap_source_for_html,
                    )
                    logger.info(
                        "%s action=write_html status=ok path=%s library=%s",
                        COMPONENT,
                        html_path,
                        suffix or "combined",
                    )
                    if pathway_htmls:
                        logger.info(
                            "%s action=write_html_pathways status=ok count=%d library=%s",
                            COMPONENT,
                            len(pathway_htmls),
                            suffix or "combined",
                        )
                except Exception as exc:
                    logger.warning(
                        "%s action=write_html status=failed path=%s error=%s library=%s",
                        COMPONENT,
                        html_path,
                        exc,
                        suffix or "combined",
                    )
    except Exception as e:
        logger.error("%s action=fail error=%s", COMPONENT, e)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

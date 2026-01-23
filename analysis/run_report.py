# -*- coding: utf-8 -*-
import argparse
import csv
import datetime
import html
import logging
import math
import os
import re
import sys
import yaml  # type: ignore[import-untyped]
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    from analysis.fda_name_map import (
        resolve_ligand_display_name,
        resolve_mapping_csv_path,
        try_load_fda_index,
    )
    from analysis.manifest_utils import extract_pocket, load_run_manifest
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
    from analysis.manifest_utils import extract_pocket, load_run_manifest

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
    base = _clean_text(row.get("ligand_base"))
    lig_file = _clean_text(row.get("ligand_file") or row.get("ligand") or "")
    if has_ligand_display:
        display = _clean_text(row.get("ligand_display"))
        if display:
            return display
    if fda_index is not None:
        display = _clean_text(resolve_ligand_display_name(base, lig_file, fda_index))
        if display:
            return display
    display = _clean_text(row.get("ligand"))
    if display:
        return display
    return base


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
    test_mode_tokens: Optional[List[str]] = None,
) -> Dict[str, Any]:
    tokens = test_mode_tokens or _resolve_test_mode_tokens(repo_root, run_id)
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
    if not has_ligand_display:
        mapping_csv = resolve_mapping_csv_path(
            repo_root, run_id, cli_value=fda_mapping_csv
        )
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


    # Group by target combo
    # Key: "<pdb_id>|<variant>|<ph_label>"
    target_groups: Dict[str, List[Dict[str, Any]]] = {}
    ligand_top5_hits: Dict[str, List[Dict[str, Any]]] = {}
    ligand_notable_hits: Dict[str, List[Dict[str, Any]]] = {}
    ligand_display_map: Dict[str, str] = {}
    multi_target_entries: Dict[str, List[Dict[str, Any]]] = {}

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
        n_rankable = len(sorted_rankable)
        for idx, (r, t_sel) in enumerate(sorted_rankable, start=1):
            lig_base = _clean_text(r.get("ligand_base"))
            lig_display = _resolve_ligand_display(r, has_ligand_display, fda_index)
            t_src = r.get("t_selected_source") or None
            pct = idx / n_rankable if n_rankable else 0.0

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

    sources: Dict[str, Any] = {
        "master_rows_csv": str(primary_master_csv.relative_to(repo_root)),
        "manifest_yaml": str(canonical_manifest_rel),
    }
    if len(master_csvs) > 1:
        sources["master_rows_csvs"] = [
            str(path.relative_to(repo_root)) for path in master_csvs
        ]

    report = {
        "run_id": run_id,
        "generated_at": datetime.datetime.utcnow().isoformat() + "Z",
        "sources": sources,
        "summary": summary,
        "targets": targets_out,
        "ligands": ligands_out,
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


def _write_html_report(
    report: Dict[str, Any],
    out_path: Path,
    highlight_queries: List[str],
    repo_root: Path,
) -> None:
    run_id = report.get("run_id", "")
    generated_at = report.get("generated_at", "")

    highlight_headers = "".join(
        f"<th>{html.escape(q)}</th>" for q in highlight_queries
    )
    highlight_rows = []
    targets = report.get("targets", {}) or {}
    for target_id in sorted(targets.keys()):
        target = targets[target_id]
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
            cells.append(f"<td>{html.escape(cell)}</td>")
        row_html = (
            f"<tr><td>{html.escape(target_id)}</td>{''.join(cells)}</tr>"
        )
        highlight_rows.append(row_html)

    multi_rows = []
    for hit in report.get("multi_target_hits", []) or []:
        multi_rows.append(
            "<tr>"
            f"<td>{html.escape(str(hit.get('ligand_display', '')))}</td>"
            f"<td>{html.escape(str(hit.get('ligand_base', '')))}</td>"
            f"<td>{html.escape(str(hit.get('targets_qualified', '')))}</td>"
            f"<td>{html.escape(_format_num(hit.get('worst_pct'), '.6f'))}</td>"
            f"<td>{html.escape(_format_num(hit.get('mean_pct'), '.6f'))}</td>"
            f"<td>{html.escape(str(hit.get('best_rank', '')))}</td>"
            f"<td>{html.escape(_format_num(hit.get('best_t_selected'), '.6g'))}</td>"
            "</tr>"
        )

    links = []
    report_yaml = repo_root / "data" / run_id / "report.yaml"
    if report_yaml.exists():
        links.append(f"<li><a href=\"{html.escape(str(report_yaml))}\">report.yaml</a></li>")
    heatmap_csv = repo_root / "data" / run_id / "heatmap_input.csv"
    if heatmap_csv.exists():
        links.append(
            f"<li><a href=\"{html.escape(str(heatmap_csv))}\">heatmap_input.csv</a></li>"
        )
    heatmap_png = repo_root / "data" / run_id / "heatmap.png"
    if heatmap_png.exists():
        links.append(
            f"<li><a href=\"{html.escape(str(heatmap_png))}\">heatmap.png</a></li>"
        )

    html_body = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>{html.escape(str(run_id))} report</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 24px; }}
    h1, h2 {{ margin-bottom: 8px; }}
    table {{ border-collapse: collapse; width: 100%; margin-bottom: 24px; }}
    th, td {{ border: 1px solid #ddd; padding: 6px 8px; text-align: left; }}
    th {{ background: #f3f3f3; }}
    .meta {{ color: #666; font-size: 0.9em; }}
  </style>
</head>
<body>
  <h1>{html.escape(str(run_id))} report</h1>
  <div class="meta">Generated at {html.escape(str(generated_at))}</div>

  <h2>Highlights</h2>
  <table>
    <thead>
      <tr>
        <th>target_id</th>
        {highlight_headers}
      </tr>
    </thead>
    <tbody>
      {''.join(highlight_rows)}
    </tbody>
  </table>

  <h2>Multi-target hits</h2>
  <table>
    <thead>
      <tr>
        <th>ligand_display</th>
        <th>ligand_base</th>
        <th>targets_qualified</th>
        <th>worst_pct</th>
        <th>mean_pct</th>
        <th>best_rank</th>
        <th>best_t_selected</th>
      </tr>
    </thead>
    <tbody>
      {''.join(multi_rows)}
    </tbody>
  </table>

  <h2>Artifacts</h2>
  <ul>
    {''.join(links)}
  </ul>
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
) -> None:
    master_csv = repo_root / "data" / run_id / "master_rows.csv"
    if not master_csv.exists():
        raise FileNotFoundError(f"Master CSV not found: {master_csv}")

    rows = []
    with master_csv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        has_ligand_display = "ligand_display" in (reader.fieldnames or [])
        for r in reader:
            rows.append(r)

    fda_index = None
    if not has_ligand_display:
        mapping_csv = resolve_mapping_csv_path(
            repo_root, run_id, cli_value=fda_mapping_csv
        )
        if mapping_csv is not None:
            fda_index = try_load_fda_index(mapping_csv)

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
        n_rankable = len(sorted_rankable)
        for idx, (row, _t_val) in enumerate(sorted_rankable, start=1):
            pct = idx / n_rankable if n_rankable else 0.0
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
        "pdb_id",
        "variant",
        "ph_label",
        "ligand_base",
        "ligand_display",
        "library",
        "t_selected",
        "pose_valid_any",
        "pose_invalid_reason_top",
        "rank",
        "pct_rank",
        "is_decoy",
        "is_control",
    ]

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
                    "pdb_id": pdb,
                    "variant": variant,
                    "ph_label": ph,
                    "ligand_base": _clean_text(row.get("ligand_base")),
                    "ligand_display": display,
                    "library": _clean_text(row.get("library")),
                    "t_selected": _clean_text(row.get("t_selected")),
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
    parser.add_argument("--extended-top-n", type=int, default=15)
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

    out_path = repo_root / "data" / run_id / "report.yaml"
    if out_path.exists() and not args.overwrite:
        logger.info("%s action=skip reason=exists path=%s", COMPONENT, out_path)
        return 0

    try:
        report = build_report(
            run_id,
            repo_root,
            top_n=args.top_n,
            extended_top_n=args.extended_top_n,
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
            test_mode_tokens=tokens,
        )
        write_yaml(report, out_path)
        logger.info("%s action=write status=ok path=%s", COMPONENT, out_path)

        if args.emit_heatmap_csv:
            heatmap_path = (
                Path(args.heatmap_csv_path)
                if args.heatmap_csv_path
                else repo_root / "data" / run_id / "heatmap_input.csv"
            )
            top_k = args.heatmap_top_k if args.heatmap_top_k > 0 else None
            try:
                _write_heatmap_input_csv(
                    repo_root,
                    run_id,
                    heatmap_path,
                    top_k,
                    args.fda_mapping_csv,
                    filter_invalid,
                )
                logger.info(
                    "%s action=write_heatmap_csv status=ok path=%s",
                    COMPONENT,
                    heatmap_path,
                )
            except Exception as exc:
                logger.warning(
                    "%s action=write_heatmap_csv status=failed path=%s error=%s",
                    COMPONENT,
                    heatmap_path,
                    exc,
                )

        if args.emit_html:
            html_path = (
                Path(args.html_path)
                if args.html_path
                else repo_root / "data" / run_id / "report.html"
            )
            highlight_queries = _parse_highlight_queries(args.highlight_ligands)
            if args.highlight_max_per_target > 0:
                highlight_queries = highlight_queries[: args.highlight_max_per_target]
            else:
                highlight_queries = []
            try:
                _write_html_report(report, html_path, highlight_queries, repo_root)
                logger.info(
                    "%s action=write_html status=ok path=%s", COMPONENT, html_path
                )
            except Exception as exc:
                logger.warning(
                    "%s action=write_html status=failed path=%s error=%s",
                    COMPONENT,
                    html_path,
                    exc,
                )
    except Exception as e:
        logger.error("%s action=fail error=%s", COMPONENT, e)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

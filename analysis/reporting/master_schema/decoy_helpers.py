# -*- coding: utf-8 -*-
import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

from analysis.reporting.report_config import resolve_config_key
from analysis.reporting.value_utils import strip_quotes
from analysis.reporting.master_schema.constants import (
    COMPONENT,
    DECOY_PREFIX_DEFAULT,
    DECOY_PREFIX_KEY,
    DUD_PREFIX_KEY,
    _DECOY_RE,
    canonical_ligand_base,
)
from chemdb.decoy_control import (
    ligand_filename_from_row,
    row_has_explicit_decoy,
    row_is_control,
    row_is_decoy_by_role,
)

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


def _normalize_library_name(value: Any) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    lowered = text.lower()
    if lowered in {"unknown", "none", "null", "na", "n/a"}:
        return ""
    return text


def _infer_library_from_source_csv(source_rel: str, decoy_prefix: str) -> str:
    basename = Path(str(source_rel or "")).name.lower()
    if not basename:
        return ""
    if basename == "consensus_reranked_scorch.csv":
        return "FDA"
    suffix = "_consensus_reranked_scorch.csv"
    if not basename.endswith(suffix):
        return ""
    token = basename[: -len(suffix)]
    if not token:
        return ""
    decoy = str(decoy_prefix or "").strip().lower()
    if token in {"dud", "decoy"} or (decoy and token == decoy):
        return "DECOY"
    return token.upper()


def _is_hidden_cache_path(path: Path) -> bool:
    return any(
        part.startswith(".") and part not in {".", ".."} for part in Path(path).parts
    )


def _resolve_library_name(
    row: Dict[str, Any], source_rel: str, decoy_prefix: str, is_decoy: bool
) -> str:
    if is_decoy:
        return "DECOY"
    raw_library = _normalize_library_name(row.get("library"))
    if raw_library:
        if raw_library.strip().lower() in {"decoy", "dud"}:
            return "DECOY"
        return raw_library.upper()
    source_library = _infer_library_from_source_csv(source_rel, decoy_prefix)
    if source_library and source_library != "DECOY":
        return source_library
    run_mode = _normalize_library_name(row.get("run_mode"))
    if run_mode and run_mode.strip().lower() not in {"decoy", "dud"}:
        return run_mode.upper()
    return "UNKNOWN"


def _discover_consensus_files(
    run_root: Path, tokens: List[str], decoy_prefix: str, logger: logging.Logger
) -> List[Path]:
    library_tokens = [tok for tok in tokens if tok not in {"dud", decoy_prefix}]
    if not library_tokens:
        library_tokens = ["fda"]

    patterns: List[str] = []
    for tok in library_tokens:
        if tok == "fda":
            patterns.append("consensus_reranked_scorch.csv")
        else:
            prefix = _sanitize_token(tok)
            patterns.append(f"{prefix}_consensus_reranked_scorch.csv")

    found: List[Path] = []
    for pattern in patterns:
        found.extend(run_root.rglob(pattern))

    if not found and "consensus_reranked_scorch.csv" not in patterns:
        found.extend(run_root.rglob("consensus_reranked_scorch.csv"))

    filtered: List[Path] = []
    decoy_names = {
        f"{decoy_prefix}_consensus_reranked_scorch.csv",
        "dud_consensus_reranked_scorch.csv",
    }
    for path in found:
        if _is_hidden_cache_path(path):
            continue
        if path.name in decoy_names:
            continue
        filtered.append(path)

    deduped: Dict[str, Path] = {}
    for path in filtered:
        deduped[str(path)] = path
    ordered = [deduped[key] for key in sorted(deduped.keys())]
    logger.info(
        "%s action=discover patterns=%s files=%d",
        COMPONENT,
        ",".join(patterns),
        len(ordered),
    )
    return ordered


def _is_decoy_file(filename: str, decoy_prefix: str) -> bool:
    name = str(filename or "")
    if _DECOY_RE.search(name):
        return True
    prefix = str(decoy_prefix or "").strip().lower()
    if not prefix:
        return False
    return name.lower().startswith(f"{prefix}_")


def _row_is_decoy(row: Dict[str, Any], decoy_prefix: str) -> bool:
    if row_is_control(row):
        return False
    if row_is_decoy_by_role(row):
        return True
    if row_has_explicit_decoy(row):
        return True
    lig_name = ligand_filename_from_row(row)
    return _is_decoy_file(lig_name, decoy_prefix)


def _canonical_ligand_base_with_prefix(lig: str, decoy_prefix: str) -> str:
    base = canonical_ligand_base(lig)
    prefix = re.escape(str(decoy_prefix or "").strip())
    if prefix:
        for pattern in (
            rf"(_{prefix}_gnina_stage\d+)$",
            rf"(_gnina_{prefix}_stage\d+)$",
            rf"(_dock6_{prefix}_stage\d+)$",
            rf"(_{prefix}_dock6_stage\d+)$",
            rf"(_{prefix}_stage\d+)$",
            rf"(__{prefix}_ledock_stage\d+)$",
            rf"(__ledock_{prefix}_stage\d+)$",
            rf"(__{prefix}_dock6_stage\d+)$",
            rf"(__dock6_{prefix}_stage\d+)$",
        ):
            base = re.sub(pattern, "", base)
    base = re.sub(r"(__ledock_stage\d+)$", "", base)
    base = re.sub(r"(__dock6_stage\d+)$", "", base)
    base = re.sub(r"(_gnina_stage\d+)$", "", base)
    base = re.sub(r"(_stage\d+)$", "", base)
    base = base.replace("__", "_")
    base = re.sub(r"_+$", "", base)
    return base

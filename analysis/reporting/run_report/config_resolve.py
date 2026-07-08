# ruff: noqa: F403, F405
from analysis.reporting.run_report.constants import *
def _configure_logging(verbose: bool) -> logging.Logger:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(message)s")
    return logging.getLogger("run-report")

def _parse_percent_searched(raw: str) -> Optional[float]:
    value = strip_quotes(str(raw or "")).strip()
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
    raw = read_config_key(path, {"PERCENT_SEARCHED", "percent_searched"})
    if raw is None:
        return None
    return _parse_percent_searched(raw)

def _resolve_percent_searched(repo_root: Path) -> Optional[float]:
    return _read_percent_searched_from_file(repo_root / "config.txt")

def _parse_filter_invalid(raw: object) -> Optional[bool]:
    if isinstance(raw, bool):
        return raw
    value = strip_quotes(str(raw or "")).strip()
    if not value:
        return None
    lowered = value.lower()
    if lowered in TEST_MODE_ON_VALUES:
        return True
    if lowered in TEST_MODE_OFF_VALUES:
        return False
    return None

def _read_filter_invalid_from_file(path: Path) -> Optional[bool]:
    raw = read_config_key(path, {"filter_invalid", "FILTER_INVALID"})
    if raw is None:
        return None
    return _parse_filter_invalid(raw)

def _resolve_filter_invalid(repo_root: Path, run_id: str) -> Optional[bool]:
    env_raw = os.environ.get("FILTER_INVALID")
    if env_raw is not None:
        parsed = _parse_filter_invalid(env_raw)
        if parsed is not None:
            return parsed
    candidates = report_config_candidates(repo_root, run_id)
    for path in candidates:
        parsed = _read_filter_invalid_from_file(path)
        if parsed is not None:
            return parsed
    return None

def _read_target_name_prefer_from_file(path: Path) -> Optional[str]:
    prefer_keys = (
        "TARGET_NAME_PREFER",
        "TARGET_NAME_PREFERENCE",
        "DUD_EVAL_TARGET_NAME_PREFER",
        "DUD_EVAL_TARGET_NAME_PREFERENCE",
    )
    raw = read_config_key(path, prefer_keys)
    if raw:
        return raw.strip().lower()
    return None

def _resolve_target_name_cfg(repo_root: Path, run_id: str) -> Dict[str, str]:
    env_value = strip_quotes(os.environ.get("TARGET_NAME_PREFER", "")).strip().lower()
    if env_value:
        return {"target_name_prefer": env_value, "repo_root": str(repo_root)}

    candidates = report_config_candidates(repo_root, run_id)
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
    value = strip_quotes(str(raw or "")).strip()
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
    value = read_config_key(path, key)
    return _parse_report_limit(value) if value is not None else None

def _resolve_report_limit(repo_root: Path, run_id: str, key: str) -> Optional[int]:
    candidates = report_config_candidates(repo_root, run_id)
    for path in candidates:
        parsed = _read_report_limit_from_file(path, key)
        if parsed is not None:
            return parsed
    return None

def _read_report_float_from_file(
    path: Path, key: str, allow_percent: bool = False
) -> Optional[float]:
    value = read_config_key(path, key)
    if value is None:
        return None
    if allow_percent:
        parsed = _parse_percent_searched(value)
        if parsed is not None:
            return parsed
    return as_float(value)

def _resolve_report_float(
    repo_root: Path, run_id: str, key: str, allow_percent: bool = False
) -> Optional[float]:
    candidates = report_config_candidates(repo_root, run_id)
    for path in candidates:
        parsed = _read_report_float_from_file(path, key, allow_percent=allow_percent)
        if parsed is not None:
            return parsed
    return None

def _read_report_bool_from_file(path: Path, key: str) -> Optional[bool]:
    value = read_config_key(path, key)
    return _parse_filter_invalid(value) if value is not None else None

def _resolve_report_bool(repo_root: Path, run_id: str, key: str) -> Optional[bool]:
    candidates = report_config_candidates(repo_root, run_id)
    for path in candidates:
        parsed = _read_report_bool_from_file(path, key)
        if parsed is not None:
            return parsed
    return None

def _read_report_choice_from_file(
    path: Path, key: str, allowed: Set[str]
) -> Optional[str]:
    candidate = str(read_config_key(path, key) or "").strip().lower()
    return candidate if candidate in allowed else None

def _resolve_report_asset_mode(
    repo_root: Path, run_id: str, cli_value: Optional[str] = None
) -> str:
    allowed = {"auto", "inline", "relative", "cdn"}
    if cli_value is not None:
        candidate = str(cli_value).strip().lower()
        if candidate in allowed:
            return candidate
    candidates = report_config_candidates(repo_root, run_id)
    for path in candidates:
        parsed = _read_report_choice_from_file(path, "REPORT_ASSET_MODE", allowed)
        if parsed is not None:
            return parsed
    return "auto"


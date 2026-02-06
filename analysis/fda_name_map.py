# -*- coding: utf-8 -*-
import re
import sys
from pathlib import Path
from typing import Any, Optional

try:
    from prep_ligands.metabolite_resolver import (  # type: ignore[import-not-found]
        load_library_index,
        resolve_corresponding_name_for_rdk,
        resolve_corresponding_name_from_text,
    )
except Exception:
    REPO_ROOT = Path(__file__).resolve().parents[1]
    if str(REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(REPO_ROOT))
    SRC_ROOT = REPO_ROOT / "src"
    if str(SRC_ROOT) not in sys.path:
        sys.path.insert(0, str(SRC_ROOT))
    try:
        from prep_ligands.metabolite_resolver import (  # type: ignore[import-not-found]
            load_library_index,
            resolve_corresponding_name_for_rdk,
            resolve_corresponding_name_from_text,
        )
    except Exception:
        # Keep report generation non-fatal when optional resolver deps (e.g., RDKit ABI)
        # are unavailable in the runtime environment.
        def load_library_index(_mapping_csv: str) -> Optional[Any]:
            return None

        def resolve_corresponding_name_for_rdk(
            _rdk_id: int, _fda_index: Any
        ) -> Optional[str]:
            return None

        def resolve_corresponding_name_from_text(
            _text: str, _fda_index: Any
        ) -> Optional[str]:
            return None

_RDK_RE = re.compile(r"rdk[_-]?(\d+)", re.IGNORECASE)


def _strip_quotes(value: str) -> str:
    stripped = str(value or "").strip()
    if (
        len(stripped) >= 2
        and stripped[0] == stripped[-1]
        and stripped[0] in ("'", '"')
    ):
        return stripped[1:-1].strip()
    return stripped


def _clean_text(value: Any) -> str:
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return ""
    return text


def _read_config_key(path: Path, key_name: str) -> Optional[str]:
    if not path.exists():
        return None
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("#") or "=" not in stripped:
                    continue
                key, raw_value = stripped.split("=", 1)
                if key.strip() != key_name:
                    continue
                value = _strip_quotes(raw_value)
                if value:
                    return value
    except Exception:
        return None
    return None


def _resolve_path(value: str, base_dir: Path) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    return path


def resolve_mapping_csv_path(
    repo_root: Path, run_id: str, cli_value: Optional[str] = None
) -> Optional[Path]:
    if cli_value is not None:
        candidate = _strip_quotes(str(cli_value))
        if candidate:
            path = _resolve_path(candidate, repo_root)
            if path.exists():
                return path

    config_candidates = [
        repo_root / "data" / run_id / "config.txt",
        repo_root / "data" / run_id / "config_snapshot.txt",
        repo_root / "manifests" / run_id / "config.txt",
        repo_root / "config.txt",
    ]

    # Priority for implicit defaults:
    # 1) FDA_MAPPING_CSV (current config key)
    # 2) backward-compatible aliases used by older flows
    key_candidates = [
        "FDA_MAPPING_CSV",
        "fda_mapping_csv",
        "fda_mapping_filled.csv",
        "IDENTIFY_MAPPING_CSV",
    ]
    for path in config_candidates:
        for key_name in key_candidates:
            value = _read_config_key(path, key_name)
            if not value:
                continue
            resolved = _resolve_path(value, path.parent)
            if resolved.exists():
                return resolved

    fallbacks = [
        repo_root / "fda_mapping_filled.csv",
        repo_root / "data" / "fda_mapping_from_pdbqt.csv",
    ]
    for path in fallbacks:
        if path.exists():
            return path
    return None


def try_load_fda_index(mapping_csv: Optional[Path]) -> Optional[Any]:
    if mapping_csv is None:
        return None
    try:
        return load_library_index(str(mapping_csv))
    except Exception:
        return None


def resolve_ligand_display_name(
    ligand_base: str, ligand_file: str, fda_index: Any
) -> str:
    base = _clean_text(ligand_base)
    lig_file = _clean_text(ligand_file)
    if fda_index is None:
        return base

    for text in (base, lig_file):
        match = _RDK_RE.search(text or "")
        if match:
            digits = match.group(1)
            rdk_id = f"rdk_{digits.zfill(7)}"
            name = _clean_text(resolve_corresponding_name_for_rdk(rdk_id, fda_index))
            if name:
                return name

    for text in (base, Path(lig_file).stem if lig_file else ""):
        name = _clean_text(resolve_corresponding_name_from_text(text, fda_index))
        if name:
            return name

    return base

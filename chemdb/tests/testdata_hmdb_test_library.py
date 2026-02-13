from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Mapping, Tuple

# Repository root: chemdb/tests -> chemdb -> repo root
REPO_ROOT = Path(__file__).resolve().parents[2]


def _root_from_cfg(
    cfg: Mapping[str, object], keys: tuple[str, ...], default_dir: Path
) -> Path:
    """
    Prefer a config dict key, then environment, else fall back to the repo-local default.
    """
    for key in keys:
        raw = cfg.get(key)
        if not raw:
            raw = os.environ.get(key)
        if raw:
            try:
                return Path(raw)
            except Exception:
                continue
    return default_dir


def infer_ligand_roots(cfg: Mapping[str, object] | None = None) -> Tuple[Path, Path]:
    """
    Return (extracted_root, prepped_root) using the same keys tests/configs already pass.
    """
    cfg_map: Mapping[str, object] = cfg or {}
    extracted_root = _root_from_cfg(
        cfg_map,
        ("EXTRACTED_LIGANDS_DIR", "LIGAND_EXTRACTED_DIR"),
        REPO_ROOT / "extracted_ligands",
    )
    prepped_root = _root_from_cfg(
        cfg_map,
        ("PREPPED_LIGANDS_DIR", "OUTPUT_LIGANDS_DIR", "PREPPED_LIGANDS_ROOT"),
        REPO_ROOT / "prepped_ligands",
    )
    return extracted_root, prepped_root


def ensure_hmdb_test_library(
    extracted_root: Path, prepped_root: Path
) -> dict[str, Path]:
    """
    Ensure hmdb_test_library_10 exists under both roots by mirroring test_library_10.

    The copy is idempotent and safe if the target directories are already present.
    """
    extracted_src = Path(extracted_root) / "test_library_10"
    prepped_src = Path(prepped_root) / "test_library_10"

    if not extracted_src.exists():
        fallback = REPO_ROOT / "extracted_ligands" / "test_library_10"
        if fallback.exists():
            extracted_src = fallback
    if not prepped_src.exists():
        fallback = REPO_ROOT / "prepped_ligands" / "test_library_10"
        if fallback.exists():
            prepped_src = fallback

    if not extracted_src.exists():
        raise FileNotFoundError(f"Missing extracted test_library_10 at {extracted_src}")
    if not prepped_src.exists():
        raise FileNotFoundError(f"Missing prepped test_library_10 at {prepped_src}")

    extracted_dst = Path(extracted_root) / "hmdb_test_library_10"
    prepped_dst = Path(prepped_root) / "hmdb_test_library_10"

    extracted_dst.parent.mkdir(parents=True, exist_ok=True)
    prepped_dst.parent.mkdir(parents=True, exist_ok=True)

    shutil.copytree(extracted_src, extracted_dst, dirs_exist_ok=True)
    shutil.copytree(prepped_src, prepped_dst, dirs_exist_ok=True)

    return {"extracted": extracted_dst, "prepped": prepped_dst}

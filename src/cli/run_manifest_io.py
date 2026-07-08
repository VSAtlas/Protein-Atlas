"""YAML load/save and manifest path helpers for run manifests."""

from __future__ import annotations

import logging
import os
import tempfile
import threading
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Tuple

from cli.run_manifest_support import require_yaml as _require_yaml
from config.output_paths import runtime_root

def load_manifest(manifest_path: Path) -> Optional[Dict[str, Any]]:
    yaml_mod = _require_yaml()
    if yaml_mod is None:
        return None

    if not manifest_path.exists():
        return {}

    try:
        with manifest_path.open("r", encoding="utf-8") as fh:
            data = yaml_mod.safe_load(fh) or {}
        if not isinstance(data, dict):
            return {}
        return data
    except Exception:
        logging.warning(
            "[run-manifest] Failed to load manifest path=%s",
            manifest_path,
            exc_info=True,
        )
        return {}


def write_manifest(manifest_path: Path, manifest: Dict[str, Any]) -> None:
    yaml_mod = _require_yaml()
    if yaml_mod is None:
        return

    tmp_path: Optional[Path] = None
    try:
        manifest_path.parent.mkdir(parents=True, exist_ok=True)

        # Unique temp per process/thread to avoid replace collisions under parallel writes
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=manifest_path.name + ".",
            suffix=f".{os.getpid()}.{threading.get_ident()}.tmp",
            dir=manifest_path.parent,
            delete=False,
        ) as fh:
            tmp_path = Path(fh.name)
            yaml_mod.safe_dump(manifest, fh, default_flow_style=False, sort_keys=False)
            fh.flush()
            os.fsync(fh.fileno())

        tmp_path.replace(manifest_path)
    except Exception:
        logging.warning(
            "[run-manifest] Failed to write manifest path=%s",
            manifest_path,
            exc_info=True,
        )
    finally:
        try:
            if tmp_path and tmp_path.exists():
                tmp_path.unlink()
        except Exception:
            pass


def manifest_lock_path(manifest_path: Path) -> Path:
    return manifest_path.with_suffix(manifest_path.suffix + ".lock")


def protein_key(
    pdb_id: str, variant_label: Optional[str], ph_tag: Optional[str]
) -> str:
    variant_norm = (variant_label or "legacy").strip().upper() or "LEGACY"
    ph_norm = (ph_tag or "base").strip()
    ph_token = ph_norm if ph_norm else "base"
    return f"{str(pdb_id).upper()}|{variant_norm}|{ph_token}"


def get_manifest_paths(cfg: Mapping[str, Any], run_id: str) -> Tuple[Path, Path]:
    """
    Return (manifest_dir, manifest_path) for the current run and ensure the
    directory exists. Never raises.
    """

    try:
        manifest_dir = runtime_root(cfg, "MANIFESTS_DIR", "manifests") / str(run_id)
        manifest_dir.mkdir(parents=True, exist_ok=True)
        manifest_path = manifest_dir / "run_manifest.yaml"
        return manifest_dir, manifest_path
    except Exception:
        logging.warning(
            "[run-manifest] Failed to prepare manifest paths run_id=%s",
            run_id,
            exc_info=True,
        )
        fallback_dir = runtime_root(cfg, "MANIFESTS_DIR", "manifests")
        return fallback_dir, fallback_dir / "run_manifest.yaml"

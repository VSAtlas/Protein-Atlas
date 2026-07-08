"""Shared config, path, and deferral state for pose-capture helpers."""

from __future__ import annotations

import os
import threading
from pathlib import Path
from typing import Optional, Sequence

from config.normalize import _to_bool, _to_float, _to_int
from config.runtime_config import load_config, validate_config
from path_router.path_router import make_paths

EXCLUDE_HET_IDS = {
    "HOH",
    "WAT",
    "NA",
    "K",
    "CL",
    "MG",
    "MN",
    "CA",
    "ZN",
    "FE",
    "CO",
    "CU",
    "NI",
    "MO",
    "SO4",
    "PO4",
    "ACT",
    "ACE",
    "IPH",
    "FMT",
    "BME",
    "MPD",
    "DMS",
    "IPA",
    "IMD",
    "DTT",
    "TRS",
    "MES",
    "HEP",
    "CIT",
    "TAR",
    "TLA",
    "GLY",
    "EDO",
    "GOL",
    "PEG",
    "GLC",
    "GAL",
    "MAN",
    "NAG",
    "BMA",
    "FUC",
    "TRE",
    "BGC",
    "BOG",
    "HEM",
    "FAD",
    "FMN",
    "NAD",
    "NAP",
    "NADH",
    "SAM",
    "SAH",
}

try:
    _cfg = load_config() or {}
    try:
        validate_config(_cfg)
    except Exception:
        pass
except Exception:
    _cfg = {}

def _cfg_value(name: str, *, env: str | None = None):
    value = _cfg.get(name, None)
    if value is None and env:
        value = os.environ.get(env)
    return value


def _cfg_bool(name: str, *, env: str | None = None, default: bool = False) -> bool:
    return _to_bool(_cfg_value(name, env=env), default=default)


def _cfg_int(name: str, *, env: str | None = None, default: int = 0) -> int:
    return _to_int(_cfg_value(name, env=env), default=default)


def _cfg_float(name: str, default: float) -> float:
    return _to_float(_cfg_value(name), default=default)

def _resolve_receptor_and_outprefix(
    cfg,
    pdb_id,
    variant=None,
    ph_token=None,
    tag="renders",
    *,
    receptor_kind: str = "pdbqt",
):
    p = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    if str(receptor_kind).lower() == "cleaned":
        receptor = p.receptor_cleaned_pdb(variant)
    else:
        receptor = p.receptor_pdbqt(variant, ph_token)
    outdir = p.docked_variant_root(variant) / tag
    outdir.mkdir(parents=True, exist_ok=True)
    return p, str(receptor), outdir

def _normalize_candidate_path(candidate, base_dir: Optional[Path]) -> Optional[Path]:
    if candidate is None:
        return None
    cand_str = str(candidate).strip()
    if not cand_str:
        return None
    cand_path = Path(cand_str).expanduser()
    if cand_path.is_absolute() or base_dir is None:
        return cand_path.resolve()
    return (base_dir / cand_path).resolve()

def _default_outprefix_name(candidates: Sequence[Optional[str]], fallback: str) -> str:
    for cand in candidates:
        if not cand:
            continue
        stem = Path(str(cand)).stem
        if stem:
            return stem
    return fallback

def _cfg_data_or_default(cfg: Optional[dict]) -> Optional[dict]:
    return cfg if cfg is not None else (_cfg or None)


def _resolve_output_prefix(
    outprefix,
    outdir: Optional[Path],
    *,
    fallback_candidates: Sequence[Optional[str]],
    fallback_name: str,
) -> Optional[Path]:
    resolved = _normalize_candidate_path(outprefix, outdir)
    if resolved is None and outdir is not None:
        resolved = (outdir / _default_outprefix_name(fallback_candidates, fallback_name)).resolve()
    if resolved is None and outprefix:
        resolved = Path(str(outprefix)).expanduser().resolve()
    if resolved is not None:
        resolved.parent.mkdir(parents=True, exist_ok=True)
    return resolved


def _resolve_receptor_path(
    receptor_path,
    default_receptor: str | None,
) -> Optional[Path]:
    router_parent = Path(default_receptor).parent if default_receptor else None
    resolved = (
        _normalize_candidate_path(receptor_path, router_parent) if receptor_path else None
    )
    if resolved is None and default_receptor:
        resolved = Path(default_receptor).resolve()
    if resolved is None and receptor_path:
        resolved = Path(str(receptor_path)).expanduser().resolve()
    return resolved


def _resolve_original_path(original_pdb, paths_obj: object | None) -> Optional[Path]:
    if not original_pdb:
        return None
    resolved = _normalize_candidate_path(original_pdb, None)
    if resolved is not None and resolved.is_file():
        return resolved
    if paths_obj is not None:
        for attr in ('root_pdb_dir', 'processed_root', 'input_root'):
            base = getattr(paths_obj, attr, None)
            if base is None:
                continue
            cand = _normalize_candidate_path(original_pdb, Path(base))
            if cand is not None and cand.is_file():
                return cand
    return Path(str(original_pdb)).expanduser().resolve()

_DEFER_MODE = _cfg_bool("DEFER_PYMOL", env="DEFER_PYMOL", default=False)
_QUEUE_PATH = (
    _cfg.get("PYMOL_DEFER_QUEUE")
    or os.environ.get("PYMOL_DEFER_QUEUE")
    or str(Path(os.environ.get("OVERALL_DIR", Path.cwd())) / "deferred_pymol_jobs.jsonl")
)
_Q_LOCK = threading.Lock()

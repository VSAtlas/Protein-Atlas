"""
APO/HOLO mode resolution and audit helpers.

This module centralizes:
- APO/HOLO mode normalization from config/env
- Variant-specific receptor path lookup
- SHA1-based APO vs HOLO dedup helpers
- Audit JSON writers for APO/HOLO decisions and usage
"""

from __future__ import annotations

import os
import json
import logging
from pathlib import Path
from typing import Dict, Optional

from path_router import make_paths, docked_dir


def _norm_path(p: str | Path) -> str:
    """Normalize path to a clean, forward-slash string for logs & keys."""
    return os.path.abspath(str(p)).replace("\\", "/")


def _clean_mode_token(s: str | None) -> str:
    s = (s or "").strip().lower()
    # normalize separators
    s = s.replace("-", "").replace("_", "")
    return s


def _debug_normalize_mode_token(tok: str | None) -> str:
    t = (tok or "").strip().lower()
    if t in {"", "none", "null", "false", "0", "legacy"}:
        return "legacy"
    if t in {"apo", "apo_only"}:
        return "apo"
    if t in {"holo", "holo_only"}:
        return "holo"
    if t in {"apo_vs_holo", "apo+holo", "both"}:
        return "apo_vs_holo"
    return f"unknown:{t}"


def resolve_apo_holo_mode(cfg: dict) -> tuple[str, list]:
    """Normalize APO/HOLO mode from config tokens."""

    import os  # if not already imported

    env_raw = os.environ.get("APO_HOLO_MODE")
    raw_value = env_raw if env_raw is not None else cfg.get("APO_HOLO_MODE")
    logging.info("[apo-holo.debug] env.APO_HOLO_MODE_raw=%r cfg.APO_HOLO_MODE_raw=%r", env_raw,
                 cfg.get("APO_HOLO_MODE"))

    token = _clean_mode_token(str(raw_value) if raw_value is not None else "")

    if token in {"", "none", "legacy", "null", "false", "0"}:
        mode, variants = "legacy", [None]
    elif token == "apo":
        mode, variants = "apo", ["APO"]
    elif token == "holo":
        mode, variants = "holo", ["HOLO"]
    elif token in {"apovsholo", "apoandholo", "both"}:
        mode, variants = "apo_vs_holo", ["APO", "HOLO"]
    else:
        logging.warning(
            "[apo-holo.debug] unknown_mode=%r defaulting=apo_vs_holo", raw_value
        )
        mode, variants = "apo_vs_holo", ["APO", "HOLO"]

    logging.info("[apo-holo] mode=%s expanded=%s", mode, variants)
    return mode, variants


# Put near other helpers
def _variant_receptor_path(pdb_id: str, variant: str | None, cfg: dict) -> str | None:
    # Return the cleaned receptor PDB path for a given variant if it exists, else None
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    rec = paths.receptor_cleaned_pdb(variant)
    return str(rec) if rec.exists() else None


def file_sha1(path: str) -> str:
    import hashlib
    h = hashlib.sha1()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def delete_variant_trees(pdb_id: str, variant: str, cfg: dict) -> None:
    # Delete processed receptor and docking trees for a specific variant (idempotent)
    import shutil
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    proc_variant_root = paths.receptor_dir(variant).parent      # processed_pdbs/<PDB>/<VARIANT>/
    dock_variant_root = paths.docked_variant_root(variant)      # docked/<PDB>/<VARIANT>/
    for d in (proc_variant_root, dock_variant_root):
        if d.exists():
            shutil.rmtree(d, ignore_errors=True)


def dedup_identical_variants(pdb_id: str, cfg: dict) -> None:
    """
    If HOLO and APO cleaned receptors are byte-identical, delete HOLO and keep APO.
    """
    holo = _variant_receptor_path(pdb_id, "HOLO", cfg)
    apo  = _variant_receptor_path(pdb_id, "APO",  cfg)

    # log actual resolved paths and existence flags up-front.
    apo_p = Path(apo) if apo else None
    holo_p = Path(holo) if holo else None
    apo_exists = apo_p.exists() if apo_p else False
    holo_exists = holo_p.exists() if holo_p else False
    logger = logging.getLogger()
    logger.info(
        "[apo-vs-holo] compare.dedup apo=%s exists=%s holo=%s exists=%s",
        (_norm_path(apo_p) if apo_p else "None"), ("T" if apo_exists else "F"),
        (_norm_path(holo_p) if holo_p else "None"), ("T" if holo_exists else "F"),
    )

    if not holo or not apo:
        logging.warning(
            "[apo-vs-holo] pdb_id=%s stage=dedup action=skip reason=missing_paths apo=%s holo=%s",
            pdb_id,
            apo,
            holo,
        )
        _record_apo_holo_decision(cfg, pdb_id, "HOLO", "missing_paths")
        return
    try:
        holo_sha = file_sha1(holo)
        apo_sha = file_sha1(apo)
    except Exception as e:
        logging.warning(
            "[apo-vs-holo] pdb_id=%s stage=dedup action=skip reason=sha_error err=%s",
            pdb_id,
            e,
        )
        _record_apo_holo_decision(cfg, pdb_id, "HOLO", "sha_error")
        return

    if holo_sha == apo_sha:
        logging.info(
            "[apo-vs-holo] identical receptors for %s; deleting HOLO (keeping APO) apo_sha=%s holo_sha=%s",
            pdb_id,
            apo_sha,
            holo_sha,
        )
        delete_variant_trees(pdb_id, "HOLO", cfg)
        _record_apo_holo_decision(cfg, pdb_id, "HOLO", "deleted_postrun")
    else:
        logging.info(
            "[apo-vs-holo] pdb_id=%s stage=dedup action=keep reason=not_identical apo_sha=%s holo_sha=%s",
            pdb_id,
            apo_sha,
            holo_sha,
        )
        _record_apo_holo_decision(cfg, pdb_id, "HOLO", "not_identical")


def _audit_variant_key(variant: str | None) -> str:
    return (variant or "LEGACY").upper()


def _ensure_apo_holo_variant_entry(cfg: Dict, pdb_id: str, variant: str | None):
    if not cfg.get("AUDIT_JSON", True):
        return None
    store = cfg.setdefault("_APO_HOLO_AUDIT", {})
    per_pdb = store.setdefault(pdb_id.upper(), {})
    key = _audit_variant_key(variant)
    entry = per_pdb.setdefault(
        key,
        {
            "variant": key,
            "env_token": (os.environ.get("APO_HOLO_VARIANT", "") or ""),
            "records": [],
            "dedup_decision": "pending",
        },
    )
    entry["env_token"] = (os.environ.get("APO_HOLO_VARIANT", "") or "")
    entry.setdefault("records", [])
    entry.setdefault("dedup_decision", "pending")
    return entry


def _flush_apo_holo_audit(cfg: Dict, pdb_id: str) -> None:
    if not cfg.get("AUDIT_JSON", True):
        return
    store = cfg.get("_APO_HOLO_AUDIT")
    if not store:
        return
    payload = store.get(pdb_id.upper())
    if not payload:
        return
    try:
        out_root = docked_dir(pdb_id, variant=None, ph_tag=None, legacy=False)
        out_root.mkdir(parents=True, exist_ok=True)
        (out_root / "apo_holo_audit.json").write_text(json.dumps(payload, indent=2))
    except Exception:
        pass


def _record_apo_holo_usage(
    cfg: Dict,
    pdb_id: str,
    variant: str | None,
    ph_label: str | None,
    receptor_path: str | Path,
) -> None:
    entry = _ensure_apo_holo_variant_entry(cfg, pdb_id, variant)
    if entry is None:
        return
    rec_path = Path(receptor_path)
    exists_flag = rec_path.exists()
    record = {
        "ph_label": None if ph_label is None else str(ph_label),
        "receptor_path": str(rec_path),
        "exists": exists_flag,
    }
    if exists_flag:
        try:
            record["sha1"] = file_sha1(str(rec_path))
        except Exception as err:
            record["sha1_error"] = str(err)
    entry.setdefault("records", []).append(record)
    _flush_apo_holo_audit(cfg, pdb_id)


def _record_apo_holo_decision(cfg: Dict, pdb_id: str, variant: str | None, decision: str) -> None:
    entry = _ensure_apo_holo_variant_entry(cfg, pdb_id, variant)
    if entry is None:
        return
    current = entry.get("dedup_decision")
    if current in {None, "", "pending"}:
        entry["dedup_decision"] = decision
    elif current != decision:
        history = entry.setdefault("decision_history", [])
        history.append(decision)
    _flush_apo_holo_audit(cfg, pdb_id)

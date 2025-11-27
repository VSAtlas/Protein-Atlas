# debug_fs.py
# Debug filesystem wrappers and helpers to collapse repeated '.sanitized' tokens.

from __future__ import annotations

import hashlib
import logging
import os
import re
import traceback
from pathlib import Path
from typing import Iterable, Set

# --- Debug wrappers to locate legacy/incorrect folder creation ---

_orig_mkdir = Path.mkdir
def _dbg_mkdir(self, *a, **k):
    path_str = str(self)
    # Log legacy cleaned_ligands creations
    if path_str.lower().endswith("_cleaned_ligands"):
        logging.error(
            "[DBG] Path.mkdir for legacy path: %s\n%s",
            path_str,
            "".join(traceback.format_stack(limit=6)),
        )
    # Log unwanted processed_pdbs/<PDB>_CLEANED creations (case-insensitive)
    if "/processed_pdbs/" in path_str and re.search(r"(?i)_cleaned/?$", path_str):
        logging.error(
            "[DBG] Path.mkdir for processed_pdbs CLEANED dir: %s\n%s",
            path_str,
            "".join(traceback.format_stack(limit=6)),
        )
    return _orig_mkdir(self, *a, **k)

_orig_makedirs = os.makedirs
def _dbg_makedirs(name, *a, **k):
    p = str(name)
    if p.lower().endswith("_cleaned_ligands"):
        logging.error(
            "[DBG] os.makedirs for legacy path: %s\n%s",
            p,
            "".join(traceback.format_stack(limit=6)),
        )
    if "/processed_pdbs/" in p and re.search(r"(?i)_cleaned/?$", p):
        logging.error(
            "[DBG] os.makedirs for processed_pdbs CLEANED dir: %s\n%s",
            p,
            "".join(traceback.format_stack(limit=6)),
        )
    return _orig_makedirs(name, *a, **k)

def install_debug_makedirs() -> None:
    """
    Install debug wrappers for Path.mkdir and os.makedirs so we can see
    any creation of legacy _cleaned/_cleaned_ligands paths.
    """
    Path.mkdir = _dbg_mkdir
    os.makedirs = _dbg_makedirs

# --- sanitize-collapse helpers ---

_SANITIZED_RUN = re.compile(r"(?:\.sanitized){2,}")

def _collapse_sanitized_token(fn: str) -> str:
    """Collapse any repeated '.sanitized' tokens anywhere in the stem."""
    stem, ext = os.path.splitext(fn)
    new_stem = _SANITIZED_RUN.sub(".sanitized", stem)
    return new_stem + ext

def _sha1(path: str, bufsize: int = 1 << 20) -> str:
    h = hashlib.sha1()
    with open(path, "rb") as f:
        while True:
            b = f.read(bufsize)
            if not b:
                break
            h.update(b)
    return h.hexdigest()

def _files_identical(a: str, b: str) -> bool:
    try:
        if os.path.getsize(a) != os.path.getsize(b):
            return False
        return _sha1(a) == _sha1(b)
    except Exception:
        return False

def collapse_sanitized_names(
    root_dirs: Iterable[object],
    exts: Set[str] = {".pdb", ".sdf", ".mol2", ".pdbqt"},
    logger=None,
) -> None:
    """
    Walk given roots and collapse repeated '.sanitized' in filenames.
    If the canonical name exists:
      - if byte-identical, delete the redundant file
      - if different, keep the canonical (shortest run) and delete the longer-run file; warn.
    """
    for root in root_dirs:
        if not root:
            continue
        root_str = os.fspath(root)
        if not os.path.isdir(root_str):
            continue
        for dirpath, _, files in os.walk(root_str):
            for fn in files:
                ext = os.path.splitext(fn)[1].lower()
                if ext not in exts:
                    continue
                new_fn = _collapse_sanitized_token(fn)
                if new_fn == fn:
                    continue
                src = os.path.join(dirpath, fn)
                dst = os.path.join(dirpath, new_fn)
                rel_src = os.path.relpath(src, root_str)
                rel_dst = os.path.relpath(dst, root_str)
                try:
                    if os.path.exists(dst):
                        if _files_identical(src, dst):
                            os.remove(src)
                            if logger:
                                logger.info(
                                    f"[sanitize-collapse] dedup: removed duplicate '{rel_src}' (kept '{rel_dst}')"
                                )
                        else:
                            # Prefer the shorter '.sanitized' run (i.e., dst). Remove the longer one.
                            os.remove(src)
                            if logger:
                                logger.warning(
                                    f"[sanitize-collapse] conflict: kept '{rel_dst}', removed longer-run '{rel_src}'"
                                )
                    else:
                        os.rename(src, dst)
                        if logger:
                            logger.info(
                                f"[sanitize-collapse] rename: '{rel_src}' -> '{rel_dst}'"
                            )
                except Exception as e:
                    if logger:
                        logger.error(
                            f"[sanitize-collapse] failed on '{rel_src}' -> '{rel_dst}': {e}"
                        )


def collapse_sanitized_names_for_cfg(cfg, logger=None) -> None:
    """
    Convenience wrapper used from main.py.

    Given the ConfigDict (or plain dict) `cfg`, gather all relevant
    directories that may contain '.sanitized' filenames and call
    collapse_sanitized_names(...) on them.

    This should be a thin wrapper: no new behavior beyond picking
    the right roots and delegating to collapse_sanitized_names.
    """
    candidate_keys = [
        "PREPPED_LIGANDS_DIR",
        "PREPPED_LIGANDS_ROOT",
        "OUTPUT_LIGANDS_DIR",
        "LIGAND_DIR",
        "LIGANDS_MOL2_DIR",
        "OUTPUT_DIR",
        "PDBQT_DIR",
    ]

    roots = []
    # Support both dict-style and attribute-style access.
    for key in candidate_keys:
        val = None
        if hasattr(cfg, "get"):
            try:
                val = cfg.get(key)  # ConfigDict path
            except Exception:
                val = None
        if val is None:
            # Fallback to attribute-style (cfg.PREPPED_LIGANDS_DIR, etc.)
            val = getattr(cfg, key, None)
        if val:
            roots.append(os.fspath(val))

    if not roots:
        return

    if logger:
        logger.info("[sanitize-collapse] scanning %d roots: %s", len(roots), roots)

    collapse_sanitized_names(roots, logger=logger)

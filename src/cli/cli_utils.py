# -*- coding: utf-8 -*-
from __future__ import annotations

import os
from typing import Iterable, Optional


def _normalize_flag_name(tok: str) -> str:
    """
    Normalize a CLI flag token by stripping leading dashes and
    lowercasing. This lets '-flag' and '--flag' be interchangeable.
    """
    return str(tok).lstrip("-").lower()


def _cli_val(argv, flag):
    """
    Return the value following a flag, treating '-flag' and '--flag'
    as equivalent. Example: _cli_val(sys.argv, "--run-id") will find
    values from either '-run-id' or '--run-id'.
    """
    try:
        target = _normalize_flag_name(flag)
        for i, tok in enumerate(argv):
            # Only consider tokens that look like flags
            if not tok.startswith("-"):
                continue
            if _normalize_flag_name(tok) == target:
                if i + 1 < len(argv) and not argv[i + 1].startswith("-"):
                    return argv[i + 1]
    except Exception:
        pass
    return None


def _cli_has(argv, flag):
    """
    Return True if a flag is present, treating '-flag' and '--flag'
    as equivalent.
    """
    try:
        target = _normalize_flag_name(flag)
        for tok in argv:
            if not tok.startswith("-"):
                continue
            if _normalize_flag_name(tok) == target:
                return True
        return False
    except Exception:
        return False


# --- Single-ligand ---
def _parse_single_from_cli(argv) -> str:
    """
    Minimal CLI parser for: --single <pattern>
    Returns the pattern string or "" if not provided.
    """
    try:
        if "--single" in argv:
            i = argv.index("--single")
            if i + 1 < len(argv) and not argv[i + 1].startswith("-"):
                return argv[i + 1]
    except Exception:
        pass
    return ""


def _parse_fast_flag(argv) -> bool:
    """Return True if argv includes fast/-fast/--fast (case-insensitive)."""
    try:
        return any(tok.lower().lstrip("-") == "fast" for tok in argv)
    except Exception:
        return False


def _norm_pdb_id(token: str) -> Optional[str]:
    """
    Normalize a user token to a 4-char PDB ID (uppercase).
    Accepts bare IDs (2HYY), QoL flags (--2HYY or -2HYY), and filenames (2HYY.pdb).
    Returns None if it cannot produce a 4-char alnum ID.
    """
    if not token:
        return None
    t = str(token).strip()
    # Strip any leading dashes (one or two)
    while t.startswith("-"):
        t = t[1:]
    t = os.path.basename(t)
    if t.lower().endswith(".pdb"):
        t = t[:-4]
    t = t.replace("_cleaned", "")
    t = t.upper()
    if len(t) >= 4:
        cand = t[:4]
        return cand if cand.isalnum() else None
    return None


def _split_ids(s: str) -> list[str]:
    """Split a comma/whitespace separated string into normalized 4-char IDs."""
    if not s:
        return []
    parts = s.replace(",", " ").split()
    out = []
    for p in parts:
        nid = _norm_pdb_id(p)
        if nid:
            out.append(nid)
    return out


def _dedupe_order(seq: Iterable[str]) -> list[str]:
    """De-duplicate while preserving first-seen order."""
    seen = set()
    out = []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def _parse_specified_proteins(argv, cfg) -> tuple[list[str], str]:
    """
    Resolve requested PDB IDs with precedence CLI > ENV > CFG.

    CLI:
      --pdb  XIAP        (repeatable)
      -pdb   XIAP        (short form, repeatable)
      --pdbs "XIAP,1BN1" (comma/space separated)
      -pdbs  "XIAP 1BN1" (short form)
      --XIAP, -XIAP      (QoL: any --<4char> alnum)
      XIAP, xiap.pdb     (bare tokens)

    ENV:
      ONLY_PDBS="XIAP 1BN1"

    CFG:
      SPECIFIED_PROTEINS: JSON list or string "XIAP, 1BN1"

    Returns: (normalized_ids, source or "")
    """
    cli_ids: list[str] = []
    consumed_value_idx: set[int] = set()

    # --- --pdb / -pdb (repeatable) ---
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok in ("--pdb", "-pdb"):
            if i + 1 < len(argv) and not argv[i + 1].startswith("-"):
                nid = _norm_pdb_id(argv[i + 1])
                if nid:
                    cli_ids.append(nid)
                    consumed_value_idx.add(i + 1)
            i += 2
            continue
        i += 1

    # --- --pdbs / -pdbs "XIAP,1BN1" ---
    for flag in ("--pdbs", "-pdbs"):
        try:
            if flag in argv:
                j = argv.index(flag)
                if j + 1 < len(argv) and not argv[j + 1].startswith("-"):
                    cli_ids.extend(_split_ids(argv[j + 1]))
                    consumed_value_idx.add(j + 1)
        except Exception:
            # be robust to weird argv
            pass

    # --- QoL: --XIAP / -XIAP style (exact length, 4-char alnum) ---
    for tok in argv:
        low = tok.lower()
        # don't treat fast/-fast/--fast as a PDB short-form token
        if low in ("fast", "-fast", "--fast"):
            continue
        if (tok.startswith("--") and len(tok) == 6) or (
            tok.startswith("-") and len(tok) == 5
        ):
            nid = _norm_pdb_id(tok)
            if nid:
                cli_ids.append(nid)

    # --- Bare tokens: XIAP, xiap.pdb, XIAP_cleaned.pdb ---
    #
    # We skip:
    #   - argv[0] (script name)
    #   - tokens we've already consumed as values to known flags
    #   - anything starting with '-' (flags)
    #
    for idx, tok in enumerate(argv[1:], start=1):
        if idx in consumed_value_idx:
            continue
        if tok.startswith("-"):
            continue

        base = os.path.basename(tok)

        # Strip common suffix patterns
        lower = base.lower()
        if lower.endswith("_cleaned.pdb"):
            core = base[: -len("_cleaned.pdb")]
        elif lower.endswith(".pdb"):
            core = base[:-4]
        else:
            core = base

        core = core.strip()
        if not core:
            continue

        # Require exactly 4 alnum chars to avoid grabbing e.g. run-id strings
        if len(core) == 4 and core.isalnum():
            nid = _norm_pdb_id(core)
            if nid:
                cli_ids.append(nid)

    if cli_ids:
        return _dedupe_order(cli_ids), "CLI"

    # --- ENV ---
    env_val = os.environ.get("ONLY_PDBS", "").strip()
    if env_val:
        return _dedupe_order(_split_ids(env_val)), "ENV"

    # --- CFG ---
    cfg_val = cfg.get("SPECIFIED_PROTEINS", "")
    cfg_list: list[str] = []
    if isinstance(cfg_val, str):
        cfg_list = _split_ids(cfg_val)
    else:
        try:
            cfg_list = [t for t in cfg_val]
        except Exception:
            cfg_list = []
    cfg_list = [_norm_pdb_id(t) for t in cfg_list if _norm_pdb_id(t)]
    if cfg_list:
        return _dedupe_order(cfg_list), "CFG"

    return [], ""

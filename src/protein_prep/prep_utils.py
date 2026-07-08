from __future__ import annotations

import json
import logging
import os
import shlex
import subprocess
from collections.abc import Mapping as MappingABC
from pathlib import Path
from typing import Iterable, List, Mapping, Optional, cast

from protein_prep.pdb_fixer_runtime import (
    get_atom_rules,
    load_canonical_cofactors,
    load_canonical_metals,
    load_canonical_waters,
)


ALIASES = get_atom_rules()
config: dict = {}
_CONFIG: dict = config
_IONS_CFG_CACHE: tuple[set[str], str] | None = None
_IONS_CFG_LOGGED = False


def set_runtime_config(cfg: dict | None) -> None:
    global config, _CONFIG, _IONS_CFG_CACHE, _IONS_CFG_LOGGED
    config = cfg if isinstance(cfg, dict) else {}
    _CONFIG = config
    _IONS_CFG_CACHE = None
    _IONS_CFG_LOGGED = False


def _to_upper_set(values: Iterable[str] | None) -> set[str]:
    out: set[str] = set()
    for val in values or []:
        if val is None:
            continue
        text = str(val).strip()
        if text:
            out.add(text.upper())
    return out


def _cfg(key: str, default: str = "", legacy_key: str | None = None) -> str:
    """env > config[key] > config[legacy_key] > default (all coerced to str)"""
    v = os.environ.get(key)
    if isinstance(v, str) and v != "":
        return v
    if key in _CONFIG and str(_CONFIG.get(key)) != "":
        return str(_CONFIG.get(key))
    if legacy_key and legacy_key in _CONFIG and str(_CONFIG.get(legacy_key)) != "":
        return str(_CONFIG.get(legacy_key))
    return str(default)


def _cfg_bool(key: str, default: bool = False) -> bool:
    v = (_cfg(key, str(int(default))) or "").strip().lower()
    return v in ("1", "true", "yes", "y", "on")


def _cfg_float(key: str, default: float) -> float:
    try:
        return float(_cfg(key, str(default)))
    except Exception:
        return float(default)


def _cfg_int(key: str, default: int) -> int:
    try:
        return int(float(_cfg(key, str(default))))
    except Exception:
        return int(default)


def _cfg_chain_keep_list() -> set[str]:
    raw = _cfg("CHAIN_KEEP_LIST", "") or ""
    toks = [t.strip() for t in raw.replace(";", ",").split(",") if t.strip()]
    return {t if len(t) == 1 else t[:1] for t in toks}


def _resolve_variant_token(
    cfg: Mapping[str, object] | None = None, override: Optional[str] = None
) -> Optional[str]:
    token = (override or "").strip().upper() if override else ""
    if token in {"APO", "HOLO"}:
        return token
    env_token = (os.environ.get("APO_HOLO_VARIANT") or "").strip().upper()
    if env_token in {"APO", "HOLO"}:
        return env_token
    if cfg is not None:
        cfg_token = str(cfg.get("_CURRENT_VARIANT", "")).strip().upper()
        if cfg_token in {"APO", "HOLO"}:
            return cfg_token
    return None


def _load_retain_allowlist(cfg: Optional[dict]) -> tuple[set[str], str]:
    global _IONS_CFG_CACHE, _IONS_CFG_LOGGED
    if _IONS_CFG_CACHE is not None:
        allow, source = _IONS_CFG_CACHE
        if not _IONS_CFG_LOGGED:
            logging.info(
                "[ions.cfg] retain_tokens=%s source=%s", ",".join(sorted(allow)), source
            )
            _IONS_CFG_LOGGED = True
        return allow, source

    def _tokens_from_attr(value: object) -> set[str]:
        if value is None:
            return set()
        if isinstance(value, MappingABC):
            return _to_upper_set(value.keys())
        if isinstance(value, str):
            return _to_upper_set([value])
        if isinstance(value, (list, tuple, set, frozenset)):
            return _to_upper_set(cast(Iterable[str], value))
        return set()

    def _load_canonical_attr(name: str, loader) -> set[str]:
        tokens = _tokens_from_attr(getattr(ALIASES, name, None))
        if tokens:
            return tokens
        try:
            loaded = loader(None)
        except Exception as exc:
            logging.warning(
                "[ions.cfg] canonical_load_failed attr=%s err=%s", name, exc
            )
            return set()
        return _to_upper_set(loaded)

    canonical_metals = _load_canonical_attr("canonical_metals", load_canonical_metals)
    canonical_cofactors = _load_canonical_attr(
        "canonical_cofactors", load_canonical_cofactors
    )
    canonical_waters = _load_canonical_attr("canonical_waters", load_canonical_waters)
    legacy_tokens = _to_upper_set(getattr(ALIASES, "retain_resnames", []))

    base_tokens = (
        canonical_metals | canonical_cofactors | canonical_waters | legacy_tokens
    )

    candidate = None
    if cfg and "retain_in_receptor_resnames" in cfg:
        candidate = cfg.get("retain_in_receptor_resnames")
    elif "retain_in_receptor_resnames" in config:
        candidate = config.get("retain_in_receptor_resnames")

    config_tokens: set[str] = set()
    if isinstance(candidate, (list, tuple, set)):
        config_tokens = _to_upper_set(candidate)
    elif isinstance(candidate, str) and candidate.strip():
        text = candidate.strip()
        parsed: Iterable[str] | None = None
        if text.startswith("[") and text.endswith("]"):
            try:
                loaded = json.loads(text)
                if isinstance(loaded, list):
                    parsed = loaded
            except Exception as exc:
                logging.warning(
                    "[ions.cfg] inline_json_parse_failed=%s err=%s", text[:40], exc
                )
        if parsed is None:
            tokens = [
                tok.strip() for tok in text.replace(";", ",").split(",") if tok.strip()
            ]
            if tokens:
                parsed = tokens
        if parsed is not None:
            config_tokens = _to_upper_set(parsed)

    allow_set = set(base_tokens) | set(config_tokens)
    source = "aliases"
    if config_tokens and base_tokens:
        source = "aliases+config"
    elif config_tokens and not base_tokens:
        source = "config"
    elif not base_tokens:
        source = "aliases"

    _IONS_CFG_CACHE = (allow_set, source)
    if not _IONS_CFG_LOGGED:
        logging.info(
            "[ions.cfg] retain_tokens=%s source=%s", ",".join(sorted(allow_set)), source
        )
        _IONS_CFG_LOGGED = True
    return allow_set, source


def _persist_subproc(
    tag: str,
    cmd: List[str],
    cp: "subprocess.CompletedProcess",
    outdir: Path,
    receptor_pdbqt: Path,
) -> None:
    """
    Save command, stdout, stderr to <work>/<tag>.* and emit a one-line summary with rc and receptor size.
    """
    try:
        outdir.mkdir(parents=True, exist_ok=True)
    except Exception:
        pass
    (outdir / f"{tag}.cmd.txt").write_text(
        " ".join(shlex.quote(x) for x in cmd), encoding="utf-8", errors="ignore"
    )
    (outdir / f"{tag}.stdout.txt").write_text(
        cp.stdout or "", encoding="utf-8", errors="ignore"
    )
    (outdir / f"{tag}.stderr.txt").write_text(
        cp.stderr or "", encoding="utf-8", errors="ignore"
    )
    exists = receptor_pdbqt.exists()
    size = receptor_pdbqt.stat().st_size if exists else 0
    logging.info(
        "[receptor-attempt %s] rc=%s exists=%s size=%d stderr=%s",
        tag,
        cp.returncode,
        exists,
        size,
        str(outdir / f"{tag}.stderr.txt"),
    )


def _short_path_for_log(path: Path) -> str:
    try:
        cwd = Path.cwd()
        return str(path.resolve(strict=False).relative_to(cwd))
    except Exception:
        try:
            return str(path.resolve(strict=False))
        except Exception:
            return str(path)

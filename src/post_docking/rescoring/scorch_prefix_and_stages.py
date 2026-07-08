from __future__ import annotations

import argparse
import logging
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

TEST_MODE_OFF_VALUES = {"", "0", "false", "no", "off", "none", "null"}
TEST_MODE_ON_VALUES = {"true", "yes", "on", "1"}
TEST_MODE_BOTH_VALUES = {
    "both",
    "fda_dud",
    "dud_fda",
    "fda+dud",
    "dud+fda",
    "fda-dud",
    "dud-fda",
}

VINA_STAGE_DIRS_BASE = ("stage1", "stage2", "stage3")
GNINA_STAGE_DIRS_BASE = ("gnina_stage1", "gnina_stage2", "gnina_stage3")
DOCK6_STAGE_DIRS_BASE = ("dock6_stage1", "dock6_stage2", "dock6_stage3")
LEDOCK_STAGE_DIRS_BASE = ("ledock_stage1", "ledock_stage2", "ledock_stage3")
POST_STAGE_DIRS_BASE = ("ledock_pdbqt", "dock6_pdbqt")


def normalize_decoy_prefix(value: Optional[object], default: str = "dud") -> str:
    prefix = str(value).strip().strip('"').strip("'") if value is not None else ""
    return prefix if prefix else default


def decoy_prefix_aliases(decoy_prefix: str) -> Tuple[str, ...]:
    prefix = normalize_decoy_prefix(decoy_prefix)
    aliases = [prefix]
    for alias in ("fda_dud", "dud"):
        if alias not in aliases:
            aliases.append(alias)
    return tuple(sorted(aliases, key=len, reverse=True))


def parse_test_mode_value(raw: object) -> list[str]:
    if isinstance(raw, bool):
        return ["dud", "fda"] if raw else ["fda"]

    s = str(raw).strip()
    if not s:
        return ["fda"]
    lowered = s.lower()
    if lowered in TEST_MODE_OFF_VALUES:
        return ["fda"]
    if lowered in TEST_MODE_ON_VALUES:
        return ["dud", "fda"]
    if lowered == "default":
        return ["fda"]
    if lowered in TEST_MODE_BOTH_VALUES:
        return ["dud", "fda"]

    tokens = [tok for tok in re.split(r"[+,\s]+", lowered) if tok]
    normalized: list[str] = []
    for tok in tokens:
        if tok in {"and", "off", "none", "null"}:
            continue
        if tok == "default":
            tok = "fda"
        normalized.append(tok)

    if not normalized:
        return ["fda"]

    deduped: list[str] = []
    seen: set[str] = set()
    for tok in normalized:
        if tok not in seen:
            seen.add(tok)
            deduped.append(tok)
    return deduped


def parse_test_mode_tokens(
    cfg: Dict[str, object],
    *,
    prefer_env: bool = False,
) -> list[str]:
    raw = cfg.get("TEST_MODE_ENABLE", "")
    if prefer_env and "TEST_MODE_ENABLE" in os.environ:
        raw = os.environ.get("TEST_MODE_ENABLE", "")
    return parse_test_mode_value(raw)


def infer_decoy_prefix_from_test_mode(
    cfg: Dict[str, object],
    *,
    default_prefix: str = "dud",
    prefer_env: bool = False,
) -> str:
    for tok in parse_test_mode_tokens(cfg, prefer_env=prefer_env):
        if tok == "fda":
            continue
        if tok == "dud":
            return default_prefix
        return normalize_decoy_prefix(tok, default_prefix)
    return default_prefix


def resolve_decoy_prefix_from_config(
    cfg: Dict[str, object],
    *,
    decoy_prefix_key: str = "DECOY_PREFIX",
    dud_prefix_key: str = "DUD_PREFIX",
    default_prefix: str = "dud",
) -> str:
    if decoy_prefix_key in cfg:
        raw = str(cfg.get(decoy_prefix_key, "")).strip()
        if raw:
            return normalize_decoy_prefix(cfg.get(decoy_prefix_key), default_prefix)
        return default_prefix
    if dud_prefix_key in cfg:
        raw = str(cfg.get(dud_prefix_key, "")).strip()
        if raw:
            return normalize_decoy_prefix(cfg.get(dud_prefix_key), default_prefix)
    return infer_decoy_prefix_from_test_mode(cfg, default_prefix=default_prefix)


def resolve_decoy_prefix_override(args: argparse.Namespace) -> Optional[str]:
    raw = getattr(args, "decoy_prefix", None)
    if raw is None:
        return None
    token = str(raw).strip()
    return token or None


def resolve_decoy_prefix_override_with_env(
    args: argparse.Namespace,
    *,
    env_keys: Sequence[str],
    default_prefix: str = "dud",
) -> Optional[str]:
    direct = resolve_decoy_prefix_override(args)
    if direct:
        return normalize_decoy_prefix(direct, default_prefix)
    for key in env_keys:
        raw = os.environ.get(key)
        if raw and str(raw).strip():
            return normalize_decoy_prefix(raw, default_prefix)
    return None


def decoy_prefixes_from_test_mode(
    cfg: Dict[str, object],
    override: Optional[str],
    *,
    default_prefix: str = "dud",
    prefer_env: bool = False,
) -> list[str]:
    if override:
        return [normalize_decoy_prefix(override, default_prefix)]

    tokens = parse_test_mode_tokens(cfg, prefer_env=prefer_env)
    fallback_prefix = resolve_decoy_prefix_from_config(
        cfg,
        default_prefix=default_prefix,
    )
    prefixes: list[str] = []
    for tok in tokens:
        if tok == "fda":
            continue
        if tok == "dud":
            prefixes.append(fallback_prefix)
        else:
            prefixes.append(normalize_decoy_prefix(tok, default_prefix))

    if not prefixes:
        prefixes = [normalize_decoy_prefix(fallback_prefix, default_prefix)]

    deduped: list[str] = []
    seen: set[str] = set()
    for prefix in prefixes:
        if prefix not in seen:
            seen.add(prefix)
            deduped.append(prefix)
    return deduped


def decoy_stage_dirs(order: Sequence[int], *, decoy_prefix: str) -> Tuple[str, ...]:
    prefix = normalize_decoy_prefix(decoy_prefix)
    return tuple(f"{prefix}_stage{int(idx)}" for idx in order)


def decoy_engine_stage_dirs(
    engine: str, order: Sequence[int], *, decoy_prefix: str
) -> Tuple[str, ...]:
    prefix = normalize_decoy_prefix(decoy_prefix)
    return tuple(f"{engine}_{prefix}_stage{int(idx)}" for idx in order)


def decoy_engine_stage_dirs_legacy(
    engine: str, order: Sequence[int], *, decoy_prefix: str
) -> Tuple[str, ...]:
    prefix = normalize_decoy_prefix(decoy_prefix)
    return tuple(f"{prefix}_{engine}_stage{int(idx)}" for idx in order)


def decoy_engine_stage_dirs_all(engine: str, *, decoy_prefix: str) -> Tuple[str, ...]:
    prefix = normalize_decoy_prefix(decoy_prefix)
    return (
        f"{engine}_{prefix}_stage1",
        f"{engine}_{prefix}_stage2",
        f"{engine}_{prefix}_stage3",
        f"{prefix}_{engine}_stage1",
        f"{prefix}_{engine}_stage2",
        f"{prefix}_{engine}_stage3",
    )


def decoy_post_stage_dirs(*, decoy_prefix: str) -> Tuple[str, ...]:
    prefix = normalize_decoy_prefix(decoy_prefix)
    return (f"{prefix}_ledock_pdbqt", f"{prefix}_dock6_pdbqt")


def vina_stage_dirs(*, decoy_prefix: str) -> Tuple[str, ...]:
    decoy_dirs = tuple(
        stage
        for alias in decoy_prefix_aliases(decoy_prefix)
        for stage in decoy_stage_dirs((1, 2, 3), decoy_prefix=alias)
    )
    return VINA_STAGE_DIRS_BASE + decoy_dirs


def gnina_stage_dirs(*, decoy_prefix: str) -> Tuple[str, ...]:
    return (
        GNINA_STAGE_DIRS_BASE
        + tuple(
            stage
            for alias in decoy_prefix_aliases(decoy_prefix)
            for stage in decoy_engine_stage_dirs("gnina", (1, 2, 3), decoy_prefix=alias)
        )
        + tuple(
            stage
            for alias in decoy_prefix_aliases(decoy_prefix)
            for stage in decoy_engine_stage_dirs_legacy(
                "gnina", (1, 2, 3), decoy_prefix=alias
            )
        )
    )


def dock6_stage_dirs(*, decoy_prefix: str) -> Tuple[str, ...]:
    return (
        DOCK6_STAGE_DIRS_BASE
        + tuple(
            stage
            for alias in decoy_prefix_aliases(decoy_prefix)
            for stage in decoy_engine_stage_dirs("dock6", (1, 2, 3), decoy_prefix=alias)
        )
        + tuple(
            stage
            for alias in decoy_prefix_aliases(decoy_prefix)
            for stage in decoy_engine_stage_dirs_legacy(
                "dock6", (1, 2, 3), decoy_prefix=alias
            )
        )
    )


def ledock_stage_dirs(*, decoy_prefix: str) -> Tuple[str, ...]:
    return (
        LEDOCK_STAGE_DIRS_BASE
        + tuple(
            stage
            for alias in decoy_prefix_aliases(decoy_prefix)
            for stage in decoy_engine_stage_dirs("ledock", (1, 2, 3), decoy_prefix=alias)
        )
        + tuple(
            stage
            for alias in decoy_prefix_aliases(decoy_prefix)
            for stage in decoy_engine_stage_dirs_legacy(
                "ledock", (1, 2, 3), decoy_prefix=alias
            )
        )
    )


def post_stage_dirs(*, decoy_prefix: str) -> Tuple[str, ...]:
    return POST_STAGE_DIRS_BASE + decoy_post_stage_dirs(decoy_prefix=decoy_prefix)


def discover_stage3_roots(variant_root: Path, *, decoy_prefix: str) -> Dict[str, Path]:
    roots: Dict[str, Path] = {}
    fda = variant_root / "stage3"
    dud = variant_root / f"{normalize_decoy_prefix(decoy_prefix)}_stage3"
    if fda.exists():
        roots["fda"] = fda
    if dud.exists():
        roots["dud"] = dud
    logger = logging.getLogger("rescoring_scorch")
    logger.info(
        "[scorch.discover] variant_root=%s fda=%s dud=%s",
        str(variant_root),
        str(fda) if fda.exists() else "None",
        str(dud) if dud.exists() else "None",
    )
    return roots


def _prefer_existing(
    root: Optional[Path], primary: Sequence[str], legacy: Sequence[str]
) -> Tuple[str, ...]:
    if root is None:
        return tuple(primary)
    if any((root / name).exists() for name in primary):
        return tuple(primary)
    if any((root / name).exists() for name in legacy):
        return tuple(legacy)
    return tuple(primary)


def stage_dir_candidates(
    source: str,
    mode: str,
    root: Optional[Path] = None,
    *,
    decoy_prefix: str,
) -> Tuple[str, ...]:
    src = (source or "").lower()
    mode_norm = (mode or "").lower()
    prefix = normalize_decoy_prefix(decoy_prefix)
    if src == "vina":
        return (
            ("stage3", "stage2", "stage1")
            if mode_norm != "dud"
            else tuple(
                stage
                for alias in decoy_prefix_aliases(prefix)
                for stage in decoy_stage_dirs((3, 2, 1), decoy_prefix=alias)
            )
        )
    if src == "gnina":
        if mode_norm != "dud":
            return ("gnina_stage3", "gnina_stage2", "gnina_stage1")
        primary = tuple(
            stage
            for alias in decoy_prefix_aliases(prefix)
            for stage in decoy_engine_stage_dirs("gnina", (3, 2, 1), decoy_prefix=alias)
        )
        legacy = tuple(
            stage
            for alias in decoy_prefix_aliases(prefix)
            for stage in decoy_engine_stage_dirs_legacy(
                "gnina", (3, 2, 1), decoy_prefix=alias
            )
        )
        return _prefer_existing(root, primary, legacy)
    if src == "dock6":
        if mode_norm != "dud":
            return ("dock6_pdbqt",)
        post_dirs = tuple(
            f"{alias}_dock6_pdbqt" for alias in decoy_prefix_aliases(prefix)
        )
        return _prefer_existing(root, post_dirs, post_dirs)
    if src == "ledock":
        if mode_norm != "dud":
            return ("ledock_pdbqt",)
        post_dirs = tuple(
            f"{alias}_ledock_pdbqt" for alias in decoy_prefix_aliases(prefix)
        )
        return _prefer_existing(root, post_dirs, post_dirs)
    return ()


def _collect_combo_from_rel(parts: Sequence[str]) -> Optional[Tuple[str, str, str]]:
    if len(parts) < 3:
        return None
    return parts[0], parts[1], parts[2]


def _collect_combo_from_stage_rel(parts: Sequence[str]) -> Optional[Tuple[str, str, str]]:
    if len(parts) < 2:
        return None
    if len(parts) == 2:
        return parts[0], "", ""
    if len(parts) == 3:
        middle = str(parts[1]).strip()
        if middle.lower().startswith("ph"):
            return parts[0], "", middle
        return parts[0], parts[1], ""
    return parts[0], parts[1], parts[2]


def discover_combos(
    run_root: Path,
    post_root: Path,
    *,
    decoy_prefix: str,
    pdb_id: Optional[str] = None,
    variant: Optional[str] = None,
    ph: Optional[str] = None,
) -> Set[Tuple[str, str, str]]:
    combos: Set[Tuple[str, str, str]] = set()
    pdb_filter = normalized_filter_token(pdb_id)
    run_search_roots = _scoped_pdb_roots(run_root, pdb_filter)
    post_search_roots = _scoped_pdb_roots(post_root, pdb_filter)
    for stage in (
        vina_stage_dirs(decoy_prefix=decoy_prefix)
        + gnina_stage_dirs(decoy_prefix=decoy_prefix)
        + dock6_stage_dirs(decoy_prefix=decoy_prefix)
        + ledock_stage_dirs(decoy_prefix=decoy_prefix)
    ):
        for search_root in run_search_roots:
            for path in search_root.rglob(stage):
                try:
                    rel = path.relative_to(run_root).parts
                except ValueError:
                    continue
                combo = _collect_combo_from_stage_rel(rel)
                if combo:
                    combos.add(combo)

    for stage_dir in post_stage_dirs(decoy_prefix=decoy_prefix):
        for search_root in post_search_roots:
            for path in search_root.rglob(stage_dir):
                try:
                    rel = path.relative_to(post_root).parts
                except ValueError:
                    continue
                combo = _collect_combo_from_stage_rel(rel)
                if combo:
                    combos.add(combo)
    return filter_combos(combos, pdb_id=pdb_id, variant=variant, ph=ph)


def _scoped_pdb_roots(root: Path, pdb_filter: Optional[str]) -> Tuple[Path, ...]:
    if not pdb_filter:
        return (root,)
    direct = root / pdb_filter
    roots: List[Path] = []
    if direct.exists():
        roots.append(direct)
    try:
        for child in root.iterdir():
            if (
                child.is_dir()
                and child.name.lower() == pdb_filter.lower()
                and child not in roots
            ):
                roots.append(child)
    except Exception:
        pass
    return tuple(roots or (direct,))


def normalized_filter_token(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    token = str(value).strip()
    return token if token else None


def combo_matches_filters(
    combo: Tuple[str, str, str],
    *,
    pdb_id: Optional[str] = None,
    variant: Optional[str] = None,
    ph: Optional[str] = None,
) -> bool:
    combo_pdb, combo_variant, combo_ph = combo
    pdb_filter = normalized_filter_token(pdb_id)
    variant_filter = normalized_filter_token(variant)
    ph_filter = normalized_filter_token(ph)
    if pdb_filter and combo_pdb.lower() != pdb_filter.lower():
        return False
    if variant_filter and combo_variant.lower() != variant_filter.lower():
        return False
    if ph_filter and combo_ph.lower() != ph_filter.lower():
        return False
    return True


def filter_combos(
    combos: Set[Tuple[str, str, str]],
    *,
    pdb_id: Optional[str] = None,
    variant: Optional[str] = None,
    ph: Optional[str] = None,
) -> Set[Tuple[str, str, str]]:
    if not (pdb_id or variant or ph):
        return combos
    return {
        combo
        for combo in combos
        if combo_matches_filters(combo, pdb_id=pdb_id, variant=variant, ph=ph)
    }


def done_sentinel_path(
    post_run_root: Path, combo: Tuple[str, str, str], *, done_dir: str, done_sentinel: str
) -> Path:
    pdb_id, variant, ph = combo
    return post_run_root / pdb_id / variant / ph / done_dir / done_sentinel


def filter_done_combos_with_stats(
    combos: Set[Tuple[str, str, str]],
    post_run_root: Path,
    overwrite: bool,
    logger: logging.Logger,
    *,
    component: str,
    done_dir: str,
    done_sentinel: str,
) -> Tuple[Set[Tuple[str, str, str]], Dict[str, int]]:
    if overwrite or not combos:
        return combos, {"combos_input": len(combos), "combos_skipped": 0}

    grouped: Dict[Tuple[str, str], List[Tuple[str, str, str]]] = defaultdict(list)
    for combo in combos:
        grouped[(combo[0], combo[1])].append(combo)

    pending: Set[Tuple[str, str, str]] = set()
    skipped = 0
    for (pdb_id, variant), ph_combos in grouped.items():
        sorted_combos = sorted(ph_combos, key=lambda item: item[2])
        done_flags = [
            done_sentinel_path(
                post_run_root,
                c,
                done_dir=done_dir,
                done_sentinel=done_sentinel,
            ).exists()
            for c in sorted_combos
        ]
        if all(done_flags):
            skipped += len(sorted_combos)
            logger.info(
                "%s action=idempotency status=skip pdb_id=%s variant=%s reason=all_ph_done n_ph=%d",
                component,
                pdb_id,
                variant,
                len(sorted_combos),
            )
            continue
        for combo, done in zip(sorted_combos, done_flags):
            if done:
                skipped += 1
                logger.info(
                    "%s action=idempotency status=skip pdb_id=%s variant=%s ph=%s reason=done_sentinel",
                    component,
                    combo[0],
                    combo[1],
                    combo[2],
                )
            else:
                pending.add(combo)

    if skipped:
        logger.info(
            "%s action=idempotency status=ok combos_pending=%d combos_skipped=%d",
            component,
            len(pending),
            skipped,
        )
    return pending, {"combos_input": len(combos), "combos_skipped": skipped}


def filter_done_combos(
    combos: Set[Tuple[str, str, str]],
    post_run_root: Path,
    overwrite: bool,
    logger: logging.Logger,
    *,
    component: str,
    done_dir: str,
    done_sentinel: str,
) -> Set[Tuple[str, str, str]]:
    pending, _stats = filter_done_combos_with_stats(
        combos,
        post_run_root,
        overwrite,
        logger,
        component=component,
        done_dir=done_dir,
        done_sentinel=done_sentinel,
    )
    return pending


def log_no_combos_after_filters(
    logger: logging.Logger,
    *,
    component: str,
    run_id: str,
    decoy_prefix: str,
    pdb_id_filter: str,
    variant_filter: str,
    ph_filter: str,
    idempotent_only: bool,
    combos_input: int = 0,
    combos_skipped: int = 0,
    filtered_scope: bool = False,
) -> None:
    if idempotent_only:
        logger.info(
            "%s action=discover status=skip reason=no_combos_after_filters_idempotent run_id=%s decoy_prefix=%s combos_input=%d combos_skipped=%d filters[pdb_id=%s variant=%s ph=%s]",
            component,
            run_id,
            decoy_prefix,
            int(combos_input),
            int(combos_skipped),
            pdb_id_filter,
            variant_filter,
            ph_filter,
        )
        return
    if filtered_scope:
        logger.info(
            "%s action=discover status=skip reason=no_combos_after_filters_filtered_scope run_id=%s decoy_prefix=%s combos_input=%d combos_skipped=%d filters[pdb_id=%s variant=%s ph=%s]",
            component,
            run_id,
            decoy_prefix,
            int(combos_input),
            int(combos_skipped),
            pdb_id_filter,
            variant_filter,
            ph_filter,
        )
        return
    logger.warning(
        "%s action=discover status=skip reason=no_combos_after_filters run_id=%s decoy_prefix=%s filters[pdb_id=%s variant=%s ph=%s]",
        component,
        run_id,
        decoy_prefix,
        pdb_id_filter,
        variant_filter,
        ph_filter,
    )


def mark_combo_done(
    post_run_root: Path,
    combo: Tuple[str, str, str],
    logger: logging.Logger,
    *,
    component: str,
    done_dir: str,
    done_sentinel: str,
) -> None:
    done_path = done_sentinel_path(
        post_run_root,
        combo,
        done_dir=done_dir,
        done_sentinel=done_sentinel,
    )
    try:
        done_path.parent.mkdir(parents=True, exist_ok=True)
        done_path.write_text("ok\n", encoding="utf-8")
    except Exception as exc:
        logger.warning(
            "%s action=idempotency status=warn reason=write_done_failed combo=%s path=%s error=%s",
            component,
            combo,
            done_path,
            exc,
        )
        return
    logger.info(
        "%s action=idempotency status=ok reason=write_done combo=%s path=%s",
        component,
        combo,
        done_path,
    )

from __future__ import annotations

import logging
import math
import re
from pathlib import Path
from typing import Any, Dict, List, Optional


def _selection_schedule(
    cfg: Dict[str, Any],
    docking_mode: str,
    n_stages: Optional[int],
    logger: Optional[logging.Logger] = None,
) -> List[float]:
    """
    Selection percentages per stage.
    Uses config overrides when present:
      - DISCOVERY_SELECTION_PCTS
      - POLYPHARM_SELECTION_PCTS
    """
    default_disc = [1.0, 0.1, 0.01, 0.001, 0.001]
    default_poly = [1.0, 0.05, 0.005]
    default_sched = {
        "discovery": default_disc,
        "polypharmacology": default_poly,
    }.get(docking_mode, [1.0] * n_stages if n_stages is not None else [1.0])

    def _fmt_selection_err(config_key: str, details: str) -> str:
        return (
            f"{config_key}: {details}. "
            "Use comma-separated values as either fractions (e.g. 1.0,0.15,0.10) "
            "or percents (e.g. 100,15,10 or 100%,15%,10%). "
            "Do not mix unsuffixed fraction and percent scales in one list."
        )

    def _tokenize_schedule(raw_val: Any, config_key: str) -> List[str]:
        if isinstance(raw_val, str):
            text = raw_val.strip()
            if not text:
                return []
            toks = [tok.strip() for tok in text.split(",")]
        elif isinstance(raw_val, (list, tuple)):
            if not raw_val:
                return []
            toks = [str(tok).strip() for tok in raw_val]
        else:
            raise ValueError(
                _fmt_selection_err(
                    config_key,
                    f"expected comma-separated string/list, got {type(raw_val).__name__}",
                )
            )

        if any(tok == "" for tok in toks):
            raise ValueError(_fmt_selection_err(config_key, "contains empty token"))
        return toks

    def _parse_schedule_values(raw_val: Any, config_key: str) -> List[float]:
        tokens = _tokenize_schedule(raw_val, config_key)
        if not tokens:
            return []

        parsed: List[float] = []
        unsuffixed_scales: set[str] = set()
        saw_explicit_percent = False
        saw_unsuffixed_fraction = False

        for idx, token in enumerate(tokens, start=1):
            explicit_percent = token.endswith("%")
            num_text = token[:-1].strip() if explicit_percent else token
            try:
                value = float(num_text)
            except Exception as e:
                raise ValueError(
                    _fmt_selection_err(
                        config_key, f"token #{idx} '{token}' is not numeric ({e})"
                    )
                ) from e

            if not math.isfinite(value):
                raise ValueError(
                    _fmt_selection_err(config_key, f"token #{idx} '{token}' is NaN/inf")
                )
            if value <= 0.0:
                raise ValueError(
                    _fmt_selection_err(
                        config_key, f"token #{idx} '{token}' must be > 0"
                    )
                )

            if explicit_percent:
                saw_explicit_percent = True
                if value > 100.0:
                    raise ValueError(
                        _fmt_selection_err(
                            config_key, f"token #{idx} '{token}' percent exceeds 100"
                        )
                    )
                parsed.append(value / 100.0)
                continue

            if value <= 1.0:
                unsuffixed_scales.add("fraction")
                saw_unsuffixed_fraction = True
                parsed.append(value)
            elif value <= 100.0:
                unsuffixed_scales.add("percent")
                parsed.append(value / 100.0)
            else:
                raise ValueError(
                    _fmt_selection_err(
                        config_key, f"token #{idx} '{token}' exceeds maximum of 100"
                    )
                )

        if len(unsuffixed_scales) > 1:
            raise ValueError(
                _fmt_selection_err(
                    config_key,
                    f"mixed unsuffixed scales are ambiguous (raw={raw_val!r})",
                )
            )
        if saw_explicit_percent and saw_unsuffixed_fraction:
            raise ValueError(
                _fmt_selection_err(
                    config_key,
                    "cannot mix explicit % tokens with unsuffixed fraction tokens (e.g. 100%,0.1 is ambiguous; use 0.1% or 10%)",
                )
            )

        return parsed

    key = {
        "discovery": "DISCOVERY_SELECTION_PCTS",
        "polypharmacology": "POLYPHARM_SELECTION_PCTS",
    }.get(docking_mode)

    sched = list(default_sched)
    if key:
        raw = cfg.get(key)
        if raw not in (None, "", [], ()):
            parsed = _parse_schedule_values(raw, key)
            if not parsed:
                raise ValueError(_fmt_selection_err(key, "no values provided"))
            if not math.isclose(parsed[0], 1.0, rel_tol=0.0, abs_tol=1e-12):
                raise ValueError(
                    _fmt_selection_err(
                        key,
                        f"first value must resolve to 1.0 (100%) because stage1 always runs on the full ligand pool (got {parsed[0]:.12g})",
                    )
                )
            sched = parsed

    if n_stages is None:
        return list(sched)
    if len(sched) < n_stages:
        sched = sched + [sched[-1]] * (n_stages - len(sched))
    elif len(sched) > n_stages:
        sched = sched[:n_stages]
    return sched


def _iter_pdbqt_dirfirst(root: Path, allowed_subdirs: Optional[set[str]] = None):
    try:
        if not root or not root.exists():
            return
        for p in root.glob("*.pdbqt"):
            yield p
        for d in root.iterdir():
            if not d.is_dir():
                continue
            if allowed_subdirs is not None and d.name not in allowed_subdirs:
                continue
            for p in d.glob("*.pdbqt"):
                yield p
    except Exception:
        for p in root.rglob("*.pdbqt"):
            yield p


def _parse_ph_dir_value(dir_name: str) -> Optional[float]:
    token = str(dir_name).strip()
    if not token:
        return None

    lowered = token.lower()
    num_pattern = re.compile(r"[+-]?\d+(?:\.\d+)?$")

    def _coerce(candidate: str) -> Optional[float]:
        try:
            val = float(candidate)
        except Exception:
            return None
        return val if 0.0 <= val <= 14.0 else None

    if lowered.startswith("ph"):
        remainder = lowered[2:].lstrip("_- ")
        if remainder and remainder[0] in "+-.0123456789":
            for candidate in (remainder, remainder.replace("_", ".")):
                if num_pattern.fullmatch(candidate):
                    parsed = _coerce(candidate)
                    if parsed is not None:
                        return parsed
            match = re.search(r"([+-]?\d+(?:\.\d+)?)", remainder)
            if match:
                parsed = _coerce(match.group(1))
                if parsed is not None:
                    return parsed

    cleaned = lowered.strip("_- ")
    cleaned_alt = cleaned.replace("_", ".")
    for candidate in (cleaned, cleaned_alt):
        if candidate and num_pattern.fullmatch(candidate):
            parsed = _coerce(candidate)
            if parsed is not None:
                return parsed

    return None


def _is_ph_dir_name(dir_name: str) -> bool:
    return _parse_ph_dir_value(dir_name) is not None


def _is_under_ph_subdir(p: Path, roots: list[Path]) -> bool:
    for root in roots:
        try:
            rel = p.resolve().relative_to(Path(root).resolve())
        except Exception:
            continue
        if not rel.parts:
            continue
        if _is_ph_dir_name(rel.parts[0]):
            return True
    return False


def select_ligands_for_next(
    cfg: Dict[str, Any],
    docking_mode: str,
    i: int,
    stages: List[Dict[str, Any]],
    scores: Dict[str, float],
    logger: logging.Logger,
    base_pool_n: Optional[int] = None,
    force_include: Optional[set[str]] = None,
) -> List[str]:
    if not scores:
        return sorted(force_include) if force_include else []

    schedule = _selection_schedule(cfg, docking_mode, len(stages), logger=logger)
    pct = schedule[i + 1] if i + 1 < len(schedule) else 0.01
    pool_n = base_pool_n if (base_pool_n is not None) else len(scores)
    k_target = max(1, int(pool_n * pct))
    k = max(1, min(k_target, len(scores)))

    next_list = [
        ligand_id for ligand_id, _ in sorted(scores.items(), key=lambda kv: kv[1])[:k]
    ]

    if force_include:
        in_set = set(next_list)
        forced_add = [
            ligand_id for ligand_id in sorted(force_include) if ligand_id not in in_set
        ]
        next_list.extend(forced_add)
        if forced_add:
            logger.info(
                "[Force-carry] Added %d extracted ligands to next stage.",
                len(forced_add),
            )

    logger.info(
        "Selected %d by score (+%d forced) = %d total (%.5f%% of base=%d).",
        k,
        len(force_include or []),
        len(next_list),
        pct * 100.0,
        pool_n,
    )
    return next_list


def compute_stage_membership_from_scores(
    cfg: Dict[str, Any],
    docking_mode: str,
    scores: Dict[str, float],
    *,
    higher_is_better: bool,
    n_stages: Optional[int] = None,
) -> Dict[int, List[str]]:
    if not scores:
        return {}

    valid_scores = {
        ligand_id: score
        for ligand_id, score in scores.items()
        if isinstance(score, (int, float)) and math.isfinite(score)
    }
    if not valid_scores:
        return {}

    order = [
        ligand_id
        for ligand_id, _ in sorted(
            valid_scores.items(), key=lambda kv: kv[1], reverse=higher_is_better
        )
    ]
    total = len(order)
    schedule = _selection_schedule(cfg, docking_mode, n_stages)

    remaining = list(order)
    stage_membership: Dict[int, List[str]] = {}
    for stage_idx in range(len(schedule), 0, -1):
        pct = schedule[stage_idx - 1]
        target = max(1, int(math.ceil(total * pct)))
        assign_count = len(remaining) if stage_idx == 1 else min(target, len(remaining))
        stage_membership[stage_idx] = remaining[:assign_count]
        remaining = remaining[assign_count:]

    return {k: stage_membership[k] for k in sorted(stage_membership)}


__all__ = [
    "_selection_schedule",
    "_iter_pdbqt_dirfirst",
    "_parse_ph_dir_value",
    "_is_ph_dir_name",
    "_is_under_ph_subdir",
    "select_ligands_for_next",
    "compute_stage_membership_from_scores",
]

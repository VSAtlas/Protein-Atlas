from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

try:
    import yaml  # type: ignore
except Exception:  # pragma: no cover - optional dependency
    yaml = None  # type: ignore

_REPO_ROOT = Path(__file__).resolve().parents[2]
_PROTEIN_CLASS_RULES: Optional[Dict[str, Any]] = None


def evaluate_druggability_for_active_site(
    cfg: Dict[str, Any],
    *,
    pdb_id: str,
    variant: Optional[str],
    ph_label: Optional[str],
    center: Tuple[float, float, float],
    box_size: Tuple[float, float, float],
    logger: logging.Logger,
) -> None:
    """
    Compatibility no-op for callers that still pass active-site geometry.

    The previous external protein-druggability evaluator is no longer part of
    the supported project surface.
    """
    logger.debug(
        "[druggability.skip] pdb_id=%s variant=%s ph=%s reason=unsupported_public_surface",
        pdb_id,
        variant,
        ph_label,
    )


def _load_protein_class_rules(logger: logging.Logger) -> Dict[str, Any]:
    global _PROTEIN_CLASS_RULES
    if _PROTEIN_CLASS_RULES is not None:
        return _PROTEIN_CLASS_RULES

    rules: Dict[str, Any] = {}
    if yaml is None:
        logger.debug("[protein-family.yaml.skip] reason=missing_yaml_module")
        _PROTEIN_CLASS_RULES = rules
        return rules

    yaml_path = _REPO_ROOT / "chemdb" / "protein_names.yaml"
    try:
        if not yaml_path.is_file():
            logger.warning("[protein-family.yaml.missing] path=%s", yaml_path)
            _PROTEIN_CLASS_RULES = rules
            return rules
        with yaml_path.open("r", encoding="utf-8") as handle:
            loaded = yaml.safe_load(handle) or {}
        if isinstance(loaded, dict):
            rules = loaded
        else:
            logger.warning(
                "[protein-family.yaml.invalid] path=%s type=%s",
                yaml_path,
                type(loaded),
            )
    except Exception as exc:
        logger.warning("[protein-family.yaml.error] err=%s", exc)
        rules = {}

    _PROTEIN_CLASS_RULES = rules
    return rules


def _infer_family_prior_from_metadata(
    uniprot_id: Optional[str],
    protein_name: Optional[str],
    protein_family: Optional[str],
    ec_numbers: List[str],
    pdb_header_text: str,
    logger: logging.Logger,
) -> Tuple[Optional[str], List[str], Optional[str]]:
    """
    Use YAML-driven canonical classes to infer a family-based difficulty prior.

    Returns (family_prior_tier, triggers, protein_class_name).
    """
    rules = _load_protein_class_rules(logger)
    classes = rules.get("classes") or {}
    if not isinstance(classes, dict) or not classes:
        return None, [], None

    parts: List[str] = []
    for val in (pdb_header_text, protein_name, protein_family, uniprot_id):
        if val:
            parts.append(str(val))
    if ec_numbers:
        parts.extend(ec_numbers)
    text = " ".join(parts).lower()

    best_class: Optional[str] = None
    best_tier: Optional[str] = None
    triggers: List[str] = []
    best_val: Optional[int] = None

    def tier_to_val(tier: str) -> int:
        return {"A": 0, "B": 1, "C": 2}.get(tier.upper(), 1)

    for class_name, info in classes.items():
        if not isinstance(info, dict):
            continue
        cname = str(class_name).strip()
        tier_prior = str(info.get("tier_prior", "B")).upper()
        if tier_prior not in {"A", "B", "C"}:
            tier_prior = "B"
        class_val = tier_to_val(tier_prior)

        class_triggers: List[str] = []
        for raw_kw in info.get("text_keywords") or []:
            kw = str(raw_kw).strip().lower()
            if kw and kw in text:
                class_triggers.append(f"{cname}:kw:{kw}")

        for raw_prefix in info.get("ec_prefixes") or []:
            prefix = str(raw_prefix).strip().lower()
            if not prefix:
                continue
            for ec in ec_numbers:
                if ec.lower().startswith(prefix):
                    class_triggers.append(f"{cname}:ec:{prefix}")
                    break

        if not class_triggers:
            continue
        if best_val is None or class_val > best_val:
            best_val = class_val
            best_class = cname
            best_tier = tier_prior
            triggers = class_triggers

    if best_class and best_tier:
        logger.info(
            "[protein-family.prior] uniprot=%s class=%s tier=%s triggers=%s",
            uniprot_id or "",
            best_class,
            best_tier,
            ",".join(triggers),
        )

    return best_tier, triggers, best_class

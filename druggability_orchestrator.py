from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from path_router import make_paths
from druggability_evaluation import get_fpocket_output_root, _summarize_fpocket_info


@dataclass
class EnginePolicy:
    tier: str
    use_gnina: bool
    use_ledock: bool
    use_dock6: bool
    reason: str = ""


def _normalize_bool(raw: Any, default: bool = False) -> bool:
    if isinstance(raw, bool):
        return raw
    if raw is None:
        return default
    s = str(raw).strip().lower()
    if s in {"1", "true", "yes", "on"}:
        return True
    if s in {"0", "false", "no", "off"}:
        return False
    return default


def _use_protein_druggability(cfg: Dict[str, Any]) -> bool:
    raw = cfg.get("USE_PROTEIN_DRUGGABILITY", cfg.get("use_protein_druggability", False))
    return _normalize_bool(raw, default=False)


def _use_gnina_flag(cfg: Dict[str, Any]) -> bool:
    raw = cfg.get("USE_GNINA", cfg.get("use_gnina", False))
    return _normalize_bool(raw, default=False)

def _use_family_prior(cfg: Dict[str, Any]) -> bool:
    raw = cfg.get("USE_PROTEIN_FAMILY_PRIOR", cfg.get("use_protein_family_prior", False))
    return _normalize_bool(raw, default=False)


def _load_druggability_metrics(
    cfg: Dict[str, Any],
    pdb_id: str,
    variant: Optional[str],
    ph_label: Optional[str],
    center: Optional[Tuple[float, float, float]],
    logger: logging.Logger,
):
    if center is None:
        logger.info(
            "[druggability.orchestrator.skip] pdb_id=%s variant=%s ph=%s reason=no_center",
            pdb_id,
            variant,
            ph_label,
        )
        return None

    try:
        paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    except Exception as exc:
        logger.warning(
            "[druggability.orchestrator.skip] pdb_id=%s variant=%s ph=%s reason=paths err=%s",
            pdb_id,
            variant,
            ph_label,
            exc,
        )
        return None

    preferred = getattr(paths, "input_pdb_path", None)
    receptor_pdb = Path(preferred) if preferred and Path(preferred).exists() else Path(
        paths.receptor_cleaned_pdb(variant)
    )
    stem = receptor_pdb.stem

    root = get_fpocket_output_root(cfg)
    info_path = root / f"{stem}_out" / f"{stem}_info.txt"
    if not info_path.exists():
        logger.info(
            "[druggability.orchestrator.skip] pdb_id=%s variant=%s ph=%s reason=missing_info path=%s",
            pdb_id,
            variant,
            ph_label,
            info_path,
        )
        return None

    return _summarize_fpocket_info(
        info_path=info_path,
        receptor_pdb=receptor_pdb,
        pdb_id=pdb_id,
        variant=variant,
        ph_label=ph_label,
        center=center,
        logger=logger,
    )


def decide_engine_policy(
    cfg: Dict[str, Any],
    pdb_id: str,
    variant: Optional[str],
    ph_label: Optional[str],
    center: Optional[Tuple[float, float, float]],
    ledock_enabled: bool,
    dock6_enabled: bool,
    logger: logging.Logger,
) -> EnginePolicy:
    """
    Decide GNINA/LeDock/DOCK6 usage for this (pdb_id, variant, ph).
    """
    gnina_global = _use_gnina_flag(cfg)
    use_pd = _use_protein_druggability(cfg)

    if not use_pd:
        policy = EnginePolicy(
            tier="legacy",
            use_gnina=bool(gnina_global),
            use_ledock=bool(ledock_enabled),
            use_dock6=bool(dock6_enabled),
            reason="legacy_no_protein_druggability",
        )
        logger.info(
            "[druggability.orchestrator.policy] pdb_id=%s variant=%s ph=%s tier=%s use_pd=0 use_gnina=%d use_ledock=%d use_dock6=%d reason=%s",
            pdb_id,
            variant,
            ph_label,
            policy.tier,
            int(policy.use_gnina),
            int(policy.use_ledock),
            int(policy.use_dock6),
            policy.reason,
        )
        return policy

    metrics = _load_druggability_metrics(
        cfg=cfg,
        pdb_id=pdb_id,
        variant=variant,
        ph_label=ph_label,
        center=center,
        logger=logger,
    )

    if metrics is None:
        policy = EnginePolicy(
            tier="fallback",
            use_gnina=bool(gnina_global),
            use_ledock=bool(ledock_enabled),
            use_dock6=bool(dock6_enabled),
            reason="missing_metrics_fallback_to_legacy",
        )
        logger.info(
            "[druggability.orchestrator.policy] pdb_id=%s variant=%s ph=%s tier=%s use_pd=1 use_gnina=%d use_ledock=%d use_dock6=%d reason=%s",
            pdb_id,
            variant,
            ph_label,
            policy.tier,
            int(policy.use_gnina),
            int(policy.use_ledock),
            int(policy.use_dock6),
            policy.reason,
        )
        return policy

    use_family = _use_family_prior(cfg)
    base_tier = metrics.tier
    family_prior = getattr(metrics, "family_prior_tier", None)

    effective_tier = base_tier

    if use_family and family_prior in {"A", "B", "C"}:
        if family_prior == "C":
            effective_tier = "C"
        elif family_prior == "A":
            if base_tier == "C":
                effective_tier = "B"
            else:
                effective_tier = base_tier
        else:
            effective_tier = base_tier
    else:
        effective_tier = base_tier

    if effective_tier == "A":
        base_tier_use_gnina = False
        base_tier_use_ledock = True
        base_tier_use_dock6 = False
    elif effective_tier == "B":
        base_tier_use_gnina = False
        base_tier_use_ledock = True
        base_tier_use_dock6 = True
    else:
        base_tier_use_gnina = True
        base_tier_use_ledock = True
        base_tier_use_dock6 = True

    use_gnina = bool(gnina_global and base_tier_use_gnina)
    use_ledock = bool(ledock_enabled and base_tier_use_ledock)
    use_dock6 = bool(dock6_enabled and base_tier_use_dock6)

    triggers_str = ",".join(metrics.triggers)
    policy = EnginePolicy(
        tier=metrics.tier,
        use_gnina=use_gnina,
        use_ledock=use_ledock,
        use_dock6=use_dock6,
        reason=f"metrics;triggers={triggers_str}",
    )
    family_triggers = getattr(metrics, "family_prior_triggers", None) or []
    family_triggers_str = ";".join(family_triggers)
    protein_class = getattr(metrics, "protein_class", None) or ""
    logger.info(
        "[druggability.orchestrator.policy] pdb_id=%s variant=%s ph=%s tier_base=%s tier_family=%s tier_effective=%s class=%s druggability=%.3f volume=%.1f open=%.3f polar_frac=%.3f has_metal=%d uniprot=%s family=%s family_triggers=%s use_pd=1 use_gnina=%d use_ledock=%d use_dock6=%d triggers=%s",
        pdb_id,
        variant,
        ph_label,
        base_tier,
        family_prior or "",
        effective_tier,
        protein_class,
        metrics.druggability,
        metrics.volume,
        metrics.openness,
        metrics.polar_fraction,
        int(metrics.has_metal),
        getattr(metrics, "uniprot_id", None) or "",
        getattr(metrics, "protein_family", None) or "",
        family_triggers_str,
        int(policy.use_gnina),
        int(policy.use_ledock),
        int(policy.use_dock6),
        triggers_str,
    )
    return policy

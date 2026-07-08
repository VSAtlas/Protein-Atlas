from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple

from config.normalize import _to_bool


@dataclass
class EnginePolicy:
    tier: str
    use_gnina: bool
    use_ledock: bool
    use_dock6: bool
    reason: str = ""


def _use_gnina_flag(cfg: Dict[str, Any]) -> bool:
    raw = cfg.get("USE_GNINA", cfg.get("use_gnina", False))
    return bool(_to_bool(raw, default=False))


def decide_engine_policy(
    *,
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
    Decide optional follow-up engine usage for this (pdb_id, variant, pH).

    Protein-druggability-based engine routing was removed from the supported
    public surface; follow-up engines now honor their direct enablement flags.
    """
    policy = EnginePolicy(
        tier="legacy",
        use_gnina=_use_gnina_flag(cfg),
        use_ledock=bool(ledock_enabled),
        use_dock6=bool(dock6_enabled),
        reason="direct_engine_flags",
    )
    logger.info(
        "[engine.policy] pdb_id=%s variant=%s ph=%s center=%s use_gnina=%d use_ledock=%d use_dock6=%d reason=%s",
        pdb_id,
        variant,
        ph_label,
        "set" if center is not None else "missing",
        int(policy.use_gnina),
        int(policy.use_ledock),
        int(policy.use_dock6),
        policy.reason,
    )
    return policy

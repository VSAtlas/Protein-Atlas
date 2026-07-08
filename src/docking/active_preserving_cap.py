from __future__ import annotations

import hashlib
from collections.abc import Callable, Sequence
from typing import TypeVar

T = TypeVar("T")

_ACTIVE_PREFIXES = (
    "actives_final",
    "active_final",
    "actives_",
    "active_",
    "actives-",
    "active-",
)


def is_dud_active_key(value: str) -> bool:
    """Return true for Atlas/DUD-style active ligand base names."""

    key = str(value or "").strip().lower()
    return any(key.startswith(prefix) for prefix in _ACTIVE_PREFIXES)


def stable_active_preserving_limit(
    items: Sequence[T],
    *,
    limit: int,
    seed: int,
    scope: str,
    key: Callable[[T], str],
) -> list[T]:
    """Deterministically cap a DUD-style library while preserving actives.

    Water-policy enrichment benchmarks are invalid if a small random cap happens
    to select only decoys. This helper reserves a small active quota whenever
    active-like names are present, then fills the rest of the cap from decoys or
    other compounds using the same deterministic hashing style as the older cap.
    """

    ordered = sorted(items, key=key)
    if limit <= 0 or len(ordered) <= limit:
        return ordered

    active_items = [item for item in ordered if is_dud_active_key(key(item))]
    other_items = [item for item in ordered if not is_dud_active_key(key(item))]
    if not active_items or not other_items:
        return _stable_limit(ordered, limit=limit, seed=seed, scope=scope, key=key)

    active_quota = min(len(active_items), max(1, int(limit) // 20))
    if other_items and active_quota >= limit:
        active_quota = max(1, limit - 1)
    other_quota = max(0, limit - active_quota)

    selected = _stable_limit(
        active_items,
        limit=active_quota,
        seed=seed,
        scope=f"{scope}|actives",
        key=key,
    )
    selected.extend(
        _stable_limit(
            other_items,
            limit=other_quota,
            seed=seed,
            scope=f"{scope}|nonactives",
            key=key,
        )
    )
    if len(selected) < limit:
        selected_keys = {key(item) for item in selected}
        remaining_actives = [
            item for item in active_items if key(item) not in selected_keys
        ]
        selected.extend(
            _stable_limit(
                remaining_actives,
                limit=limit - len(selected),
                seed=seed,
                scope=f"{scope}|active_fill",
                key=key,
            )
        )
    return sorted(selected[:limit], key=key)


def _stable_limit(
    items: Sequence[T],
    *,
    limit: int,
    seed: int,
    scope: str,
    key: Callable[[T], str],
) -> list[T]:
    if limit <= 0 or len(items) <= limit:
        return sorted(items, key=key)
    ranked = sorted(
        items,
        key=lambda item: hashlib.sha1(
            f"{seed}|{scope}|{key(item)}".encode("utf-8")
        ).hexdigest(),
    )
    return sorted(ranked[:limit], key=key)
